#include <Rcpp.h>

#include "survey_sim.hpp"
#include "utilities.hpp"

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;
	
	function("rgamma_cv", &rgamma_cv);
	function("rnbinom_cv", &rnbinom_cv);
	function("rbeta_cv", &rbeta_cv);

	function("Rcpp_survey_sim_std", &survey_sim_std);
	function("Rcpp_survey_sim_nstd", &survey_sim_nstd);

}
