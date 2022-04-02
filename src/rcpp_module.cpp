#include <Rcpp.h>

#include "survey_sim.hpp"
#include "utilities.hpp"

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;
	
	function("rgamma_cv", &rgamma_cv);
	function("rnbinom_cv", &rnbinom_cv);
	function("rbeta_cv", &rbeta_cv);

	function("Rcpp_survey_sim", &survey_sim);

}
