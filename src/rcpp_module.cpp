#include <Rcpp.h>

#include "survey_sim.h"
#include "utilities.h"
#include "CountSummariseTest.h"

RCPP_MODULE(eggSimModule){

	using namespace Rcpp;

	function("Rcpp_rgamma_cv", &rgamma_cv);
	function("Rcpp_rnbinom_cv", &rnbinom_cv);
	function("Rcpp_rbeta_cv", &rbeta_cv);
	function("Rcpp_results_levels", &results_levels);

	function("Rcpp_survey_sim", &survey_sim);

	class_<CountSummariseTest>("CountSummariseTest")

    .constructor("Default constructor")

		.method("show", &CountSummariseTest::show, "The show method")
		.method("add_counts", &CountSummariseTest::add_counts, "Add counts")
		.property("result", &CountSummariseTest::get_result, "Get result")

	;

}
