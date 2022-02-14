#include "spider.h"


// KITE 2D HAVING THE ATTACK ANGLE AS AGGREGATE STATE


/* Constructor */
Spider3d::Spider3d(const param& params, Wind3d* wind, std::mt19937& generator) :
Spider{params, generator}, wind{wind} {


    // SETTING SPECIFIC PARAMETERS
    init_xspider = params.d.at("init_xspider");
    init_yspider = params.d.at("init_yspider");
    init_zspider = params.d.at("init_zspider");

    vrel_par = vecd(3, 0);
    vrel_perp = vecd(3, 0);
    v_rel = vecd(3, 0);
    u = vecd(3, 0);

    vel_traj = vec2d(0);
    pos_traj = vec2d(0);
    distance_traj = vecd(0);
    duration_traj = vecd(0);
    p = std::vector<std::vector<bead> > (N_lines, std::vector<bead> (N_beads));

    //debug_file.open("debug.txt");
}

/* Initial configuration given the initial position of the spider*/
void Spider3d::reset_spider(){
    std::normal_distribution<double> distribution(1.0, 0.01);
    std::normal_distribution<double> distribution2(0.0, 0.01);

    time_sec = 0;

    spider.q[0] = init_xspider;
    spider.q[1] = init_yspider;
    spider.q[2] = init_zspider;

    for (size_t i = 0; i < 3; i++) {
      spider.vel[i] = 0;
    }

    double* v_wind = (*wind).init(spider.q[0], spider.q[1], spider.q[2]);

    for (size_t i = 0; i < N_lines; i++) {
      for (size_t j = 0; j < N_beads; j++) {
          double a = distribution(m_generator);
          double b = distribution2(m_generator);
          double c = distribution2(m_generator);
          p[i][j].q[0] = spider.q[0]+s0*(j+1);// + s0*a;                       //initialize positions
          p[i][j].q[1] = spider.q[1];// + s0*b;
          p[i][j].q[2] = spider.q[2];// + s0*c;
      }
    }
}


const std::string Spider3d::descr() const {
    return "3d spider composed of body and multiple lines. " + wind->descr();
}

bool Spider3d::integrate_trajectory(bool therm) {

    compute_Fel();
    //compute_FKP();
    compute_Faer();
    update_state(therm);

    time_sec += dt_i;

    //if (int(time_sec/dt_i)%1000==0){
      //std::cout << time_sec << " " << spider.q[0] << " " << spider.F_el[0] << " " <<  spider.F_drag[0] << '\n';
    //}

    if (therm) return false;

    // Check if a terminal state is reached. Spider fallen
    if ((spider.q[2] <= fall_limit) || (time_sec > ep_length*dt_i)){
        std::cout << "The spider has fallen after " << time_sec << " seconds with x = " << spider.q[0] << 'm' << '\n';
        duration_traj.push_back(time_sec);
        distance_traj.push_back(spider.q[0]);
        return true;
    }
    return false;
}


void Spider3d::compute_distances(){
  double d = 0;
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 1; j < N_beads; j++) {
      for (size_t k = 0; k < 3; k++) {
        d += (p[i][j].q[k]-p[i][j-1].q[k])*(p[i][j].q[k]-p[i][j-1].q[k]);
      }
      p[i][j].r = sqrt(d);
      d = 0;
    }
    p[i][0].r = 0;
    for (size_t k = 0; k < 3; k++) {
      p[i][0].r += sqrt((p[i][0].q[k]-spider.q[k])*(p[i][0].q[k]-spider.q[k]));
    }
  }
};


void Spider3d::compute_angles(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads-1; j++) {
      p[i][j].theta = acos((p[i][j+1].q[2] - p[i][j].q[2]) / p[i][j+1].r);
      p[i][j].phi = atan2((p[i][j+1].q[1] - p[i][j].q[1]), (p[i][j+1].q[0] - p[i][j].q[0]));
    }
    spider.theta = acos((p[i][0].q[2] - spider.q[2]) / p[i][0].r);
    spider.phi = atan2((p[i][0].q[1] - spider.q[1]), (p[i][0].q[0] - spider.q[0]));
  }
};


void Spider3d::compute_Fel(){

  compute_distances();
  compute_angles();

  for (size_t i = 0; i < 3; i++) {
    spider.F_el[i] = 0;
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads-1; j++) {
      p[i][j].F_el[0] = k*(p[i][j+1].q[0] - p[i][j].q[0] - s0*cos(p[i][j].phi)*sin(p[i][j].theta));
      p[i][j].F_el[1] = k*(p[i][j+1].q[1] - p[i][j].q[1] - s0*sin(p[i][j].phi)*sin(p[i][j].theta));
      p[i][j].F_el[2] = k*(p[i][j+1].q[2] - p[i][j].q[2] - s0*cos(p[i][j].theta));
    }

    for (size_t k = 0; k < 3; k++) {
      p[i][N_beads-1].F_el[k] = 0;
    }

    for (size_t j = 1; j < N_beads; j++) {
      p[i][j].F_el[0] -= k*(p[i][j].q[0] - p[i][j-1].q[0] - s0*cos(p[i][j-1].phi)*sin(p[i][j-1].theta));
      p[i][j].F_el[1] -= k*(p[i][j].q[1] - p[i][j-1].q[1] - s0*sin(p[i][j-1].phi)*sin(p[i][j-1].theta));
      p[i][j].F_el[2] -= k*(p[i][j].q[2] - p[i][j-1].q[2] - s0*cos(p[i][j-1].theta));
    }

    p[i][0].F_el[0] -= k*(p[i][0].q[0] - spider.q[0] - s0*cos(spider.phi)*sin(spider.theta));
    p[i][0].F_el[1] -= k*(p[i][0].q[1] - spider.q[1] - s0*sin(spider.phi)*sin(spider.theta));
    p[i][0].F_el[2] -= k*(p[i][0].q[2] - spider.q[2] - s0*cos(spider.theta));

    spider.F_el[0] += k*(p[i][0].q[0] - spider.q[0] - s0*cos(spider.phi)*sin(spider.theta));
    spider.F_el[1] += k*(p[i][0].q[1] - spider.q[1] - s0*sin(spider.phi)*sin(spider.theta));
    spider.F_el[2] += k*(p[i][0].q[2] - spider.q[2] - s0*cos(spider.theta));
  }
};


void Spider3d::compute_scalarprod(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 1; j < N_beads-1; j++) {
      p[i][j].sp = 0;
      for (size_t k = 0; k < 3; k++) {
        p[i][j].sp += (p[i][j-1].q[k] - p[i][j].q[k]) * (p[i][j+1].q[k]-p[i][j].q[k]);
      }
    }

    p[i][0].sp = 0;
    for (size_t k = 0; k < 3; k++) {
      p[i][0].sp += (spider.q[k] - p[i][0].q[k]) * (p[i][1].q[k]-p[i][0].q[k]);
    }
  }
};


void Spider3d::compute_FKP(){

  compute_scalarprod();

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 2; j < N_beads; j++) {
      for (size_t k = 0; k < 3; k++) {
        p[i][j].F_KP[k] = -J*((p[i][j-2].q[k] - p[i][j-1].q[k])/(p[i][j-1].r*(p[i][j].r)) - p[i][j-1].sp*(p[i][j].q[k] - p[i][j-1].q[k])/(p[i][j-1].r*pow(p[i][j].r, 3)));
      }
    }

    for (size_t k = 0; k < 3; k++) {
      spider.F_KP[k] = -J*((p[i][1].q[k] - p[i][0].q[k])/(p[i][0].r*p[i][1].r) - p[i][0].sp*(spider.q[k] - p[i][0].q[k])/(pow(p[i][0].r, 3)*(p[i][1].r)));
      p[i][0].F_KP[k] = -J*((2*p[i][0].q[k] - spider.q[k] - p[i][1].q[k])/(p[i][0].r*p[i][1].r) + p[i][0].sp*((spider.q[k] - p[i][0].q[k])/(pow(p[i][0].r, 3)*p[i][1].r) + (p[i][1].q[k] - p[i][0].q[k])/(p[i][0].r*pow(p[i][1].r, 3))));;
      p[i][1].F_KP[k] = -J*((spider.q[k] - p[i][0].q[k])/(p[i][0].r*(p[i][1].r)) - p[i][0].sp*(p[i][1].q[k] - p[i][0].q[k])/(p[i][0].r*pow(p[i][1].r, 3)));
    }

    for (size_t j = 1; j < N_beads-1; j++) {
      for (size_t k = 0; k < 3; k++) {
        p[i][j].F_KP[k] -= J*((2*p[i][j].q[k] - p[i][j-1].q[k] - p[i][j+1].q[k])/(p[i][j].r*p[i][j+1].r) + p[i][j].sp*((p[i][j-1].q[k] - p[i][j].q[k])/(pow(p[i][j].r, 3)*p[i][j+1].r) + (p[i][j+1].q[k] - p[i][j].q[k])/(p[i][j].r*pow(p[i][j+1].r, 3))));
      }
    }

    for (size_t j = 0; j < N_beads-2; j++) {
      for (size_t k = 0; k < 3; k++) {
        p[i][j].F_KP[k] -= J*((p[i][j+2].q[k] - p[i][j+1].q[k])/(p[i][j+1].r*p[i][j+2].r) - p[i][j+1].sp*(p[i][j].q[k] - p[i][j+1].q[k])/(pow(p[i][j+1].r, 3)*(p[i][j+2].r)));
      }
    }
  }
};


void Spider3d::getWindVel(){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      p[i][j].wind_vel = (*wind).velocity(p[i][j].q[0], p[i][j].q[1], p[i][j].q[2], time_sec);
    }
  }

  spider.wind_vel = (*wind).velocity(spider.q[0], spider.q[1], spider.q[2], time_sec);
};


void Spider3d::compute_Faer(){

  double norm;
  double vrel_par_sp;
  getWindVel();

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 1; j < N_beads-1; j++) {
      norm = 0;
      vrel_par_sp = 0;
      for (size_t k = 0; k < 3; k++) {
        u[k] = p[i][j+1].q[k] - p[i][j-1].q[k];
        norm += u[k]*u[k];
      }
      for (size_t k = 0; k < 3; k++) {
        u[k] = u[k]/sqrt(norm);
        v_rel[k] = p[i][j].vel[k] - p[i][j].wind_vel[k];
        vrel_par_sp += v_rel[k]*u[k];
      }
      for (size_t k = 0; k < 3; k++) {
        vrel_par[k] = vrel_par_sp*u[k];
        vrel_perp[k] = v_rel[k] - vrel_par[k];
        p[i][j].F_drag[k] = -b*vrel_par[k] - 2*b*vrel_perp[k];
      }
    }
  }

  for (size_t i = 0; i < 3; i++) {
    v_rel[i] = spider.vel[i] - spider.wind_vel[i];
    spider.F_drag[i] = -b_spider*v_rel[i];
  }
};


void Spider3d::update_state(bool therm){
  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      for (size_t k = 0; k < 3; k++) {
        p[i][j].acc[k] = (p[i][j].F_el[k] + p[i][j].F_KP[k] + p[i][j].F_drag[k])/m_bead;
        p[i][j].vel[k] += p[i][j].acc[k]*dt_i;
        p[i][j].q[k] += p[i][j].vel[k]*dt_i;
      }
    }
  }

  if (!therm){
    for (size_t i = 0; i < 3; i++) {
      spider.acc[i] = (spider.F_el[i] + spider.F_KP[i] + spider.F_weight[i] + spider.F_drag[i])/m_spider;
      spider.vel[i] += spider.acc[i]*dt_i;
      spider.q[i] += spider.vel[i]*dt_i;
    }
  }
};


void Spider3d::build_traj(){
  vecd v = vecd(0);
  for (size_t i = 0; i < 3; i++) {
    v.push_back(spider.wind_vel[i]);
  }
  vel_traj.push_back(v);
};


void Spider3d::debug_traj(){
  //vecd F_aer = vecd(0);
  //vecd F_drag = vecd(0);
  //vecd F_el = vecd(0);
  //vecd F_KP = vecd(0);
  vecd pos = vecd(0);

  for (size_t i = 0; i < 3; i++) {
      pos.push_back(spider.q[i]);
  }

  for (size_t i = 0; i < N_lines; i++) {
    for (size_t j = 0; j < N_beads; j++) {
      for (size_t k = 0; k < 3; k++) {
        pos.push_back(p[i][j].q[k]);
      }
    }
  }

  pos_traj.push_back(pos);
}


void Spider3d::print_traj(std::string out_dir) const{
  std::ofstream out_v, out_d, out_db;
  out_v.open(out_dir + "vel_traj.txt");
  out_d.open(out_dir + "distance_traj.txt");
  out_db.open(out_dir + "debug_traj.txt");

  for (size_t i = 0; i < vel_traj.size(); i++) {
    for (size_t j = 0; j < 3; j++) {
      out_v << vel_traj[i][j] << " ";
    }
    out_v << "\n";
  }

  out_d << "Episode_length\tReached distance\n";

  for (size_t i = 0; i < distance_traj.size(); i++) {
    out_d << duration_traj[i] << "\t" << distance_traj[i] << "\n";
  }

  for (size_t i = 0; i < pos_traj.size(); i++) {
    for (size_t j = 0; j < (N_lines*N_beads+1)*3; j++) {
      out_db << pos_traj[i][j] << "\t";
      if ((j+1) % 3 == 0){
        out_db << "\n";
      }
    }
    out_db << "\n";
  }

  out_v.close();
  out_d.close();
  out_db.close();
}
