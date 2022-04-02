#include <Rcpp.h>

#include "enums.hpp"
#include "survey_class.hpp"

template<designs design, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
Rcpp::DataFrame survey_template(const Rcpp::IntegerVector& all_ns, const Rcpp::DataFrame& parameters,
								                const Rcpp::IntegerVector& n_individ, const bool summarise)
{
	const int np = parameters.nrow();
	const int ni = n_individ.length();

	const int n_day_screen = all_ns["n_day_screen"];
	const int n_aliquot_screen = all_ns["n_aliquot_screen"];
	const int n_day_pre = all_ns["n_day_pre"];
	const int n_aliquot_pre = all_ns["n_aliquot_pre"];
	const int n_day_post = all_ns["n_day_post"];
	const int n_aliquot_post = all_ns["n_aliquot_post"];
  
	const Rcpp::IntegerVector replicateID = parameters["replicateID"];

	const Rcpp::NumericVector mu_pre = parameters["mu_pre"];
	const Rcpp::NumericVector reduction = parameters["reduction"];
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

	// Create the output vectors:
	Rcpp::NumericVector efficacy(np*ni);
	Rcpp::NumericVector n_screen(np*ni);
	Rcpp::NumericVector n_pre(np*ni);
	Rcpp::NumericVector n_post(np*ni);
	Rcpp::NumericVector time_count(np*ni);

	Rcpp::NumericVector consumables_cost(np*ni);
	Rcpp::NumericVector salary_cost(np*ni);
	Rcpp::NumericVector travel_cost(np*ni);
	Rcpp::NumericVector total_cost(np*ni);

  // TODO: method can be either specified (in which case count_int/coef/add/mult must be constant and as expected) or custom
  // TODO: the parent function could check to see if e.g. aliquot_cv is always 0 and specify dist_aliquot as rpois
  
  // Create the survey class once:
  survey_class<design, method, dist_individ, dist_day, dist_aliquot, dist_red>
    survey(n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post);
  
	// Memory offset for indexing results for different n_individ:
	ptrdiff_t offset = 1L;

	int ind = 0L;
	for(int p=0L; p<np; ++p)
	{

    survey.run(n_individ, mu_pre[p], reduction[p], individ_cv[p], day_cv[p], aliquot_cv[p], reduction_cv[p],
          count_intercept[p], count_coefficient[p], count_add[p], count_mult[p],
          &efficacy[ind], &n_screen[ind], &n_pre[ind], &n_post[ind], &time_count[ind], offset);
          
		// Do cost calculations:
		for(int i=0; i<ni; ++i)
		{
			const double cons_cost = n_screen[ind] * cost_consumables_screen[p] +
										n_pre[ind] * cost_consumables_pre[p] +
										n_post[ind] * cost_consumables_post[p];

			const double totaltime = n_screen[ind] * time_consumables_screen[p] +
										n_pre[ind] * time_consumables_pre[p] +
										n_post[ind] * time_consumables_post[p] +
										time_count[ind];

			const double ndays = totaltime / (n_technicians[p] * 4.0 * 60.0 * 60.0);
			const double sal_cost = ndays * n_team[p] * cost_salary[p];

			const int tra_cost = ndays * cost_travel[p];

			consumables_cost[ind] = cons_cost;
			salary_cost[ind] = sal_cost;
			travel_cost[ind] = tra_cost;
			total_cost[ind] = cons_cost + sal_cost + tra_cost;

			ind++;
		}

	}

	// Create the final output:
	Rcpp::Function expgrd("expand.grid");
	// NB: deliberately backwards as we use expand.grid and not expand_grid:
	Rcpp::DataFrame repmean = expgrd(n_individ, replicateID);
	Rcpp::IntegerVector nind = repmean[0L];
	Rcpp::IntegerVector repID = repmean[1L];

	Rcpp::DataFrame df = Rcpp::DataFrame::create( 	Rcpp::_["replicateID"] = repID, Rcpp::_["n_individ"] = nind,
													Rcpp::_["efficacy"] = efficacy, Rcpp::_["n_screen"] = n_screen, Rcpp::_["n_pre"] = n_pre,
		 											Rcpp::_["n_post"] = n_post, Rcpp::_["time_count"] = time_count,
													Rcpp::_["consumables_cost"] = consumables_cost, Rcpp::_["salary_cost"] = salary_cost,
													Rcpp::_["travel_cost"] = travel_cost, Rcpp::_["total_cost"] = total_cost);
	return df;
}
