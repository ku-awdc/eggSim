#include "survey_sim.hpp"

#include "survey_ns.hpp"
#include "survey_ss.hpp"
#include "survey_ssr.hpp"
#include <Rcpp.h>

// TODO: compile-time polymorphism for survey design and n_day/n_aliquot


Rcpp::DataFrame survey_sim_nstd(const std::string& design,
								const Rcpp::DataFrame& parameters,
								const bool summarise)
{
	const int n = parameters.nrow();

	// TODO: meanepg_weight_recovery should be pre-computed
	
	const Rcpp::IntegerVector n_day_screen = parameters["n_day_screen"];
	const Rcpp::IntegerVector n_aliquot_screen = parameters["n_aliquot_screen"];
	const Rcpp::IntegerVector n_day_pre = parameters["n_day_pre"];
	const Rcpp::IntegerVector n_aliquot_pre = parameters["n_aliquot_pre"];
	const Rcpp::IntegerVector n_day_post = parameters["n_day_post"];
	const Rcpp::IntegerVector n_aliquot_post = parameters["n_aliquot_post"];
	const Rcpp::NumericVector n_individ = parameters["n_individ"];
	const Rcpp::NumericVector mu_pre = parameters["mean_epg"];
	const Rcpp::NumericVector reduction = parameters["reduction"];
	const Rcpp::NumericVector weight = parameters["weight"];
	const Rcpp::NumericVector performance = parameters["recovery"];
	const Rcpp::NumericVector individ_cv = parameters["individ_cv"];
	const Rcpp::NumericVector day_cv = parameters["day_cv"];
	const Rcpp::NumericVector aliquot_cv = parameters["aliquot_cv"];
	const Rcpp::NumericVector reduction_cv = parameters["reduction_cv"];
	const Rcpp::NumericVector count_intercept = parameters["count_intercept"];
	const Rcpp::NumericVector count_coefficient = parameters["count_coefficient"];
	const Rcpp::NumericVector count_add = parameters["count_add"];
	const Rcpp::NumericVector count_mult = parameters["count_mult"];

	const Rcpp::NumericVector cost_consumables_screen = parameters["cost_consumables_screen"];
	const Rcpp::NumericVector cost_consumables_pre = parameters["cost_consumables_pre"];
	const Rcpp::NumericVector cost_consumables_post = parameters["cost_consumables_post"];
	const Rcpp::NumericVector time_consumables_screen = parameters["time_consumables_screen"];
	const Rcpp::NumericVector time_consumables_pre = parameters["time_consumables_pre"];
	const Rcpp::NumericVector time_consumables_post = parameters["time_consumables_post"];
	const Rcpp::NumericVector n_technicians = parameters["n_technicians"];
	const Rcpp::NumericVector n_team = parameters["n_team"];
	const Rcpp::NumericVector cost_salary = parameters["cost_salary"];
	const Rcpp::NumericVector cost_travel = parameters["cost_travel"];
	
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

	Rcpp::NumericVector consumables_cost(n);
	Rcpp::NumericVector salary_cost(n);
	Rcpp::NumericVector travel_cost(n);
	Rcpp::NumericVector total_cost(n);

	// Do cost calculations:
	for(int i=0L; i<n; ++i)
	{
		const double cons_cost = n_screen[i] * cost_consumables_screen[i] + n_pre[i] * cost_consumables_pre[i] + n_post[i] * cost_consumables_post[i];
		
		const double totaltime = n_screen[i] * time_consumables_screen[i] + n_pre[i] * time_consumables_pre[i] + n_post[i] * time_consumables_post[i] + time_count[i];
		const double ndays = totaltime / (n_technicians[i] * 4.0 * 60.0 * 60.0);
		const double sal_cost = ndays * n_team[i] * cost_salary[i];
		
		const int tra_cost = ndays * cost_travel[i];
		
		consumables_cost[i] = cons_cost;
		salary_cost[i] = sal_cost;
		travel_cost[i] = tra_cost;
		total_cost[i] = cons_cost + sal_cost + tra_cost;
	}
	
	Rcpp::DataFrame df = Rcpp::DataFrame::create( Rcpp::_["efficacy"] = efficacy, Rcpp::_["n_screen"] = n_screen, Rcpp::_["n_pre"] = n_pre,
		 											Rcpp::_["n_post"] = n_post, Rcpp::_["time_count"] = time_count, 
													Rcpp::_["consumables_cost"] = consumables_cost, Rcpp::_["salary_cost"] = salary_cost,
													Rcpp::_["travel_cost"] = travel_cost, Rcpp::_["total_cost"] = total_cost);
	return df;
}


Rcpp::DataFrame survey_sim_std(const std::string& design, const Rcpp::DataFrame& parameters, const bool summarise)
{
	Rcpp::DataFrame rv;

	if( design == "NS_2x2" )
	{
		rv = survey_sim_nstd("NS", parameters, summarise);
	}
	else
	{
		Rcpp::stop("Unrecognised survey design");
	}

	return rv;
}
