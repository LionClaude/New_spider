#ifndef WIND_H
#define WIND_H

#include "utils.h"

/* Abstract class */
class Wind3d {
    protected:
        /* Velocity variable that is referenced to by velocity method*/
        double m_vel[3];

    public:
        virtual double* init(double x0, double y0, double z0) { return velocity(x0,y0,z0,0.0); };
        virtual double* velocity(double x, double y, double z, double t) = 0;
        virtual const std::string descr() const = 0;
};


/* Linear wind with gaussian bumps*/
class Wind3d_bump : public Wind3d {

    private:

        int N_sec;

        int N_sec_therm;

        int N_steps_therm_s;

        int frac_dt;

        int t_s;

        double dt_s;

        double dt_i;

        int N_steps_s;

        vecd f;

        int N_events;

        double mu;

        const double sigma = 1.0;

        const double sigma_x = 1.0;

        const double cost = sigma_x*sqrt(2);

        const double cost_pi = sigma_x*sqrt(PI/2);

        const double w = 0.8;

        const double lambda = 0.4;

        std::mt19937 generator;

        std::poisson_distribution<int> poisson;
        std::uniform_real_distribution<double> uniform;
        std::normal_distribution<double> gauss;

        double psi(double x);
        double d_psi(double x);

    public:
        Wind3d_bump(const param& params, std::mt19937& generator);
        virtual double* init(double x, double y, double z);
        virtual double* velocity(double x, double y, double z, double t);

        const std::string descr() const
        { return "3d linear wind with gaussian bumps."; }
};


Wind3d* get_wind3d(const param& params, std::mt19937& generator);

#endif
