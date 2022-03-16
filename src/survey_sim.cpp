#include "survey_sim.hpp"

#include "survey_ns.hpp"
#include <Rcpp.h>

// TODO: compile-time polymorphism for survey design and n_day/n_aliquot


Rcpp::DataFrame survey_sim_nstd(const std::string design, 
								const int n_day_pre, const int n_aliquot_pre, 
								const int n_day_post, const int n_aliquot_post, 
								const Rcpp::DataFrame parameters)
{
	const int n = parameters.nrow();

	const Rcpp::NumericVector n_individ = parameters["n_individ"];
	const Rcpp::NumericVector mu_pre = parameters["mu_pre"];
	const Rcpp::NumericVector reduction = parameters["reduction"];
	const Rcpp::NumericVector weight = parameters["weight"];
	const Rcpp::NumericVector performance = parameters["performance"];
	const Rcpp::NumericVector cost_sample = parameters["cost_sample"];
	const Rcpp::NumericVector cost_aliquot = parameters["cost_aliquot"];
	const Rcpp::NumericVector individ_cv = parameters["individ_cv"];
	const Rcpp::NumericVector day_cv = parameters["day_cv"];
	const Rcpp::NumericVector aliquot_cv = parameters["aliquot_cv"];
	const Rcpp::NumericVector reduction_cv = parameters["reduction_cv"];
	
	Rcpp::NumericVector efficacy(n);
	Rcpp::NumericVector cost(n);
	
	for(int i=0L; i<n; ++i)
	{
		// TODO: return std::array instead
		Rcpp::NumericVector rv = survey_ns(n_individ[i], n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post,
											mu_pre[i], reduction[i], weight[i], performance[i],
											cost_sample[i], cost_aliquot[i], individ_cv[i], 
											day_cv[i], aliquot_cv[i], reduction_cv[i]);
											
		efficacy[i] = rv[0L];
		cost = rv[1L];											
	}
	
	Rcpp::DataFrame df = Rcpp::DataFrame::create( Rcpp::_["efficacy"] = efficacy, Rcpp::_["cost"] = cost );
	return df;
}


Rcpp::DataFrame survey_sim_std(const std::string design, const Rcpp::DataFrame parameters)
{
	Rcpp::DataFrame rv;
	
	if( design == "NS_2x2" )
	{
		rv = survey_sim_nstd("NS", 2L, 2L, 2L, 2L, parameters);
	}
	else
	{
		Rcpp::stop("Unrecognised survey design");
	}
	
	return rv;
}
