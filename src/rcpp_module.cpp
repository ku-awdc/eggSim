#include <Rcpp.h>

#include "survey_ns.hpp"

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;
	function("Rcpp_survey_ns", &survey_ns);

}
