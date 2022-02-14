#ifndef ENV_H
#define ENV_H

#include "utils.h"



/* Abstract. It contains all the info of the single player MDP to be solved
   by the Natural Actor Critic algorithm with state aggregation */
class Environment {

    protected:

        /* Random number generator */
        std::mt19937 m_generator;

    public:

        // CONSTRUCTOR
        Environment(const param& par, std::mt19937& generator) : m_generator{generator} {};
        virtual ~Environment() {};

        /* Abstract. Get the description of the environment */
        virtual const std::string descr() const = 0;

        virtual void evolve() = 0;

        virtual void build_traj() = 0;

        virtual void print_traj(std::string out_dir) const = 0;

};


#endif
