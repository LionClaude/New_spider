#include "../lib/envs/spider.h"


Environment* get_env(std::string env_name, const param& params, std::mt19937& generator);

void run(Environment* env, const param& params);

void print_output(Environment* env, std::string dir);


int main(int argc, char const *argv[]) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);

  std::string env_name(argv[1]), data_dir = "../data/";
  param params_env = parse_param_file(data_dir+env_name+"/param_env.txt");

  Environment* env = get_env(env_name, params_env, generator); // Def in kite.h
  std::cout << "Environment successfully built:\n" << (*env).descr() << "\n\n";

  run(env, params_env);

  print_output(env, data_dir + env_name + "/");
  std::cout << "Trajectories successfully printed at " << data_dir + env_name + "/" << "\n";

  delete env;

  return 0;
}


Environment* get_env(std::string env_name, const param& params, std::mt19937& generator) {
		Wind3d* wind = get_wind3d(params, generator);
    Environment* env = new Spider3d(params, wind, generator);
    return env;
}


void run(Environment* env, const param& params) {

    int n_episodes;
    try {
        n_episodes = params.d.at("n_episodes");
    } catch (std::exception) {
        throw std::invalid_argument("Invalid temporal parameters of the algorithm.");
    }


    for (size_t i = 0; i < n_episodes; i++) {
        std::cout << "Episode " << i << '\n';
        (*env).evolve();
    }

}


void print_output(Environment* env, std::string dir) {
  (*env).print_traj(dir);
}
