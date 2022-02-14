#ifndef SPIDER_H
#define SPIDER_H

#include "../env.h"
#include "../wind.h"

struct bead {
  vecd q;
  vecd vel;
  vecd acc;
  double* wind_vel;
  vecd F_el;
  vecd F_KP;
  vecd F_drag;
  vecd F_weight;
  double r;                  //distances between beads
  double sp;                 //scalar products
  double theta;
  double phi;
  double h;
  double m;
  bead() {
    q = vecd(3, 0);
    vel = vecd(3, 0);
    acc = vecd(3, 0);
    F_el = vecd(3, 0);
    F_KP = vecd(3, 0);
    F_drag = vecd(3, 0);
    F_weight = vecd(3, 0);
    F_weight[2] = -m_bead*g;
    r = 0;
    theta = 0;
    phi = 0;
    sp = 0;
  }
};


class Spider : public Environment {

	private:

		// INTERNAL VARIABLE
		/* Current step of the episode */
		int curr_ep_step;

		// SET BY CONSTRUCTOR
		/* Episode length */

    int therm_length;

    int steps_btw_series;




	protected:

		// FIXED CONSTANTS
		/* Spider mass, kg */
		const double m_spider = 0.01;

    const int N_beads = 15;

    const int N_lines = 3;

    const double k = 20.0;

    const double J = 0.01;

    const double s0 = 0.1;

    const double b = 0.05;

    const double b_spider = 0.1;

		// SET BY CONSTRUCTOR

		double init_xspider;
		/* Initial x position of the spider*/
		double init_yspider;
		/* Initial y position of the spider*/
    double init_zspider;
		/* Initial z position of the spider*/
		bool xrand_start;
		/* Flag for random start on the x or not*/
		bool yrand_start;
		/* Flag for random start on the y or not*/
    double dt_i;

    double dt_s;

    int ep_length;

		std::ofstream debug_file;

	public:

		/* Constructor */
		Spider(const param& params, std::mt19937& generator);

        void evolve();




		/* Spider description */
		virtual const std::string descr() const = 0;

		/* Abstract. Set the specific kite in the initial state */
        virtual void reset_spider() = 0;

        /* Abstract. One temporal step of the internal dynamics */
        virtual bool integrate_trajectory(bool therm) = 0;

        virtual void build_traj() = 0;

        virtual void print_traj(std::string out_dir) const = 0;

        virtual void debug_traj() = 0;



};

class Spider3d : public Spider {

	protected:

		/* Height below which the spider is considered fallen */
		const double fall_limit = 0.0;

		/* Generic wind type */
		Wind3d* wind;

    std::vector<std::vector<bead> > p;
    bead spider;

    double mu;
    double h;
    double t;

    int t_s;

    vecd vrel_par;
    vecd vrel_perp;
    vecd v_rel;
    vecd u;
    double time_sec;

    vec2d vel_traj;
    vec2d pos_traj;
    vecd distance_traj;
    vecd duration_traj;

		// AUXILIARY FUNCTION FOR THE DYNAMICS
    void compute_distances();
    void compute_angles();
    void compute_Fel();
    void compute_scalarprod();
    void compute_FKP();
    void getWindVel();
    void compute_Faer();
    void update_state(bool therm);
    //void addNoise();

	public:

		Spider3d(const param& params, Wind3d* wind, std::mt19937& generator);
		virtual ~Spider3d() { delete wind; }
		virtual const std::string descr() const;

        virtual void reset_spider();
        virtual bool integrate_trajectory(bool therm);
        virtual void build_traj();
        virtual void print_traj(std::string out_dir) const;
        virtual void debug_traj();

};

#endif
