#include <Rcpp.h>

#include "survey_sim.hpp"
#include "utilities.hpp"

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;
	
	function("Rcpp_rgamma_cv", &rgamma_cv);
	function("Rcpp_rnbinom_cv", &rnbinom_cv);
	function("Rcpp_rbeta_cv", &rbeta_cv);
	function("Rcpp_count_time", &count_time);

	function("Rcpp_survey_sim", &survey_sim);

}
