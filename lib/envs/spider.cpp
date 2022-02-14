#include "spider.h"


Spider::Spider(const param& params, std::mt19937& generator) : Environment{params, generator} {

	// Check temporal inconsistencies
	if (params.d.at("ep_length") < params.d.at("series_time") || params.d.at("ep_length") < params.d.at("int_steps")
		|| params.d.at("series_time") < params.d.at("int_steps") || params.d.at("ep_length") < params.d.at("therm_length"))
		throw std::runtime_error ( "Temporal scales of the spider are not consistent\n" );

	try {
		ep_length = int((params.d.at("ep_length")+params.d.at("therm_length"))/params.d.at("int_steps"));
    therm_length = int(params.d.at("therm_length")/params.d.at("int_steps"));
		steps_btw_series = int(params.d.at("series_time")/params.d.at("int_steps"));
		dt_i = params.d.at("int_steps");
    dt_s = params.d.at("series_time");
	} catch (std::exception) {
		throw std::invalid_argument("Invalid temporal parameter");
	}
}


void Spider::evolve(){

  reset_spider();

  bool done;
  for (size_t i = 0; i < therm_length; i++) {
    done = integrate_trajectory(true);
		if (i%steps_btw_series == 0){
				build_traj();
				debug_traj();
		}
  }

  for (size_t i = 0; i < ep_length-therm_length; i++) {
    done = integrate_trajectory(false);
		debug_traj();
		
    if (done) break;
  }
}
