#include <Rcpp.h>

#include "survey_ns.hpp"
#include "survey_sim.hpp"

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;
	
	function("Rcpp_survey_ns", &survey_ns);
	
	function("Rcpp_survey_sim_std", &survey_sim_std);
	function("Rcpp_survey_sim_nstd", &survey_sim_nstd);

}
