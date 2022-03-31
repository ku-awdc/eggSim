#include "survey_sim.hpp"

#include "survey_ns.hpp"
#include "survey_ss.hpp"
#include "survey_ssr.hpp"
#include <Rcpp.h>

// TODO: compile-time polymorphism for survey design and n_day/n_aliquot


Rcpp::DataFrame survey_sim_nstd(const std::string& design,
								const Rcpp::DataFrame& parameters)
{
	const int n = parameters.nrow();

	const Rcpp::IntegerVector n_day_screen = parameters["n_day_screen"];
	const Rcpp::IntegerVector n_aliquot_screen = parameters["n_aliquot_screen"];
	const Rcpp::IntegerVector n_day_pre = parameters["n_day_pre"];
	const Rcpp::IntegerVector n_aliquot_pre = parameters["n_aliquot_pre"];
	const Rcpp::IntegerVector n_day_post = parameters["n_day_post"];
	const Rcpp::IntegerVector n_aliquot_post = parameters["n_aliquot_post"];
	const Rcpp::NumericVector n_individ = parameters["n_individ"];
	const Rcpp::NumericVector mu_pre = parameters["mu_pre"];
	const Rcpp::NumericVector reduction = parameters["reduction"];
	const Rcpp::NumericVector weight = parameters["weight"];
	const Rcpp::NumericVector performance = parameters["performance"];
	const Rcpp::NumericVector individ_cv = parameters["individ_cv"];
	const Rcpp::NumericVector day_cv = parameters["day_cv"];
	const Rcpp::NumericVector aliquot_cv = parameters["aliquot_cv"];
	const Rcpp::NumericVector reduction_cv = parameters["reduction_cv"];
	const Rcpp::NumericVector count_intercept = parameters["count_intercept"];
	const Rcpp::NumericVector count_coefficient = parameters["count_coefficient"];
	const Rcpp::NumericVector count_add = parameters["count_add"];
	const Rcpp::NumericVector count_mult = parameters["count_mult"];

	Rcpp::NumericVector efficacy(n);
	Rcpp::IntegerVector n_screen(n);
	Rcpp::IntegerVector n_pre(n);
	Rcpp::IntegerVector n_post(n);
	Rcpp::NumericVector time_count(n);

	if( design == "NS" )
	{
		for(int i=0L; i<n; ++i)
		{
			survey_ns(n_individ[i], n_day_pre[i], n_aliquot_pre[i], n_day_post[i], n_aliquot_post[i],
												mu_pre[i], reduction[i], weight[i], performance[i],
												individ_cv[i], day_cv[i], aliquot_cv[i], reduction_cv[i],
												count_intercept[i], count_coefficient[i], count_add[i], count_mult[i],
												efficacy[i], n_pre[i], n_post[i], time_count[i]);

			n_screen[i] = 0L;
		}
	}else if( design == "SS" )
	{
		for(int i=0L; i<n; ++i)
		{
			survey_ss(n_individ[i], n_day_pre[i], n_aliquot_pre[i], n_day_post[i], n_aliquot_post[i],
												mu_pre[i], reduction[i], weight[i], performance[i],
												individ_cv[i], day_cv[i], aliquot_cv[i], reduction_cv[i],
												count_intercept[i], count_coefficient[i], count_add[i], count_mult[i],
												efficacy[i], n_pre[i], n_post[i], time_count[i]);

			n_screen[i] = 0L;
		}
	}else if( design == "SSR" )
	{
		for(int i=0L; i<n; ++i)
		{
			survey_ssr(n_individ[i], n_day_screen[i], n_aliquot_screen[i], n_day_pre[i], n_aliquot_pre[i], n_day_post[i], n_aliquot_post[i],
												mu_pre[i], reduction[i], weight[i], performance[i],
												individ_cv[i], day_cv[i], aliquot_cv[i], reduction_cv[i],
												count_intercept[i], count_coefficient[i], count_add[i], count_mult[i],
												efficacy[i], n_screen[i], n_pre[i], n_post[i], time_count[i]);
		}
	}else
	{
		Rcpp::stop("Unrecognised design type");
	}

	Rcpp::DataFrame df = Rcpp::DataFrame::create( Rcpp::_["efficacy"] = efficacy, Rcpp::_["n_screen"] = n_screen, Rcpp::_["n_pre"] = n_pre,
		 											Rcpp::_["n_post"] = n_post, Rcpp::_["time_count"] = time_count);
	return df;
}


Rcpp::DataFrame survey_sim_std(const std::string& design, const Rcpp::DataFrame& parameters)
{
	Rcpp::DataFrame rv;

	if( design == "NS_2x2" )
	{
		rv = survey_sim_nstd("NS", parameters);
	}
	else
	{
		Rcpp::stop("Unrecognised survey design");
	}

	return rv;
}
