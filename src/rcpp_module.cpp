#include <Rcpp.h>

#include "survey_ns.hpp"
#include "survey_sim.hpp"

// TODO: templates with compile-time polymorphism for common day/aliquot types

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;
	function("Rcpp_survey_ns", &survey_ns);
	function("Rcpp_survey_sim", &survey_sim);

}
