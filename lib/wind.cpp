#include "wind.h"


Wind3d_bump::Wind3d_bump(const param& params, std::mt19937& generator) :
generator{generator} {
    N_sec = params.d.at("ep_length");
    N_sec_therm = params.d.at("therm_length");
    dt_s = params.d.at("series_time");
    dt_i = params.d.at("int_steps");

    N_steps_s = int(N_sec/dt_s);
    N_steps_therm_s = int(N_sec_therm/dt_s);
    f = vecd(N_steps_s, 0);
    frac_dt = int (dt_s/dt_i);

    poisson = std::poisson_distribution<int>(lambda*N_sec);
    uniform = std::uniform_real_distribution<double>(0.0, N_sec);
    gauss = std::normal_distribution<double>(0.0, 1.0);
}


double* Wind3d_bump::init(double x, double y, double z) {
    N_events = poisson(generator);

    double t = 0;

    for (size_t i = 0; i < N_events; i++) {
      mu = uniform(generator);
      double a = gauss(generator);
      for (size_t j = 0; j < N_steps_s; j++) {
        t = j*dt_s;
        f[j] += a * exp(-(t-mu)*(t-mu) / (2*sigma*sigma));
      }
    }

    return velocity(x, y, z, 0);
}


double* Wind3d_bump::velocity(double x, double y, double z, double t) {

  t_s = int(t/dt_s);
  m_vel[0] = w*z - f[t_s]*psi(x)*z;
  m_vel[1] = 0.0;
  m_vel[2] = f[t_s]*d_psi(x)*z*z/2;

  return m_vel;
}


double Wind3d_bump::d_psi(double x){
  return exp((-pow(x,2)) / (2*pow(sigma_x,2)));
};


double Wind3d_bump::psi(double x){
  return cost_pi*(erf(x/cost)-1);
};


Wind3d* get_wind3d(const param& params, std::mt19937& generator) {

    std::string wind_type = params.s.at("wind_type");
    if (wind_type == "bump") { // Gaussian bumps on linear profile
        return new Wind3d_bump(params, generator);
    }
    /*
    if (wind_type == "turbo") { // Turbolent flow
        return new Wind3d_turbo(params);
    }*/

    else throw std::invalid_argument( "Invalid wind type" );
}
