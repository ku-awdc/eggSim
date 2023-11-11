#include <Rcpp.h>
#include <math.h>

#include "enums.h"
#include "survey_class.h"

template<bool t_summarise, designs design, bool t_fixed_n, int nd0, int na0, int nd1, int na1, int nd2, int na2,
          methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
Rcpp::DataFrame survey_template(const Rcpp::IntegerVector& all_ns, const Rcpp::DataFrame& parameters,
								                const Rcpp::IntegerVector& n_individ, const bool summarise)
{
	const Rcpp::IntegerVector scenario_int = parameters["scenario_int"];

	const int np = parameters.nrow();
	const int ni = n_individ.length();
	const int ns = Rcpp::max(scenario_int) + 1L;
  if(Rcpp::min(scenario_int)!=0L) Rcpp::stop("Invalid scenario_int");

  // If summarising then we only need one output per n_individ+scenario:
  const int ol = t_summarise ? (ni*ns) : (np*ni);

	const int n_day_screen = all_ns["n_day_screen"];
	const int n_aliquot_screen = all_ns["n_aliquot_screen"];
	const int n_day_pre = all_ns["n_day_pre"];
	const int n_aliquot_pre = all_ns["n_aliquot_pre"];
	const int n_day_post = all_ns["n_day_post"];
	const int n_aliquot_post = all_ns["n_aliquot_post"];
  const int min_pos_screen = all_ns["min_positive_screen"];
  const int min_pos_pre = all_ns["min_positive_pre"];

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

  // TODO: the parent function could check to see if e.g. aliquot_cv is always 0 and specify dist_aliquot as rpois

  // Create the survey class once:
  survey_class<design, t_fixed_n, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
    survey(n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post, min_pos_screen, min_pos_pre);

	// Memory offset for indexing results for different n_individ:
	ptrdiff_t offset = 1L;

  // Eventual output:
  Rcpp::DataFrame df;

  if constexpr(t_summarise)
  {
  	const Rcpp::NumericVector cutoff = parameters["cutoff"];

  	// Create the output vectors:
    // Rcpp::IntegerVector result(ol);
  	Rcpp::NumericVector n_screen(ol);
  	Rcpp::NumericVector n_pre(ol);
  	Rcpp::NumericVector n_post(ol);
  	Rcpp::NumericVector efficacy(ol);
  	Rcpp::NumericVector efficacy_var(ol);

  	Rcpp::NumericVector mean_pre(ol);
  	Rcpp::NumericVector mean_post(ol);
  	Rcpp::NumericVector imean_pre(ol);
  	Rcpp::NumericVector imean_post(ol);
  	//Rcpp::NumericVector time_screen(ol);
  	//Rcpp::NumericVector time_pre(ol);
  	//Rcpp::NumericVector time_post(ol);

  	Rcpp::NumericVector mean_days(ol);

  	//Rcpp::NumericVector consumables_cost(ol);
  	//Rcpp::NumericVector salary_cost(ol);
  	//Rcpp::NumericVector travel_cost(ol);
  	Rcpp::NumericVector mean_cost(ol);
  	Rcpp::NumericVector var_cost(ol);

    Rcpp::IntegerVector n_result_0(ol);
    Rcpp::IntegerVector n_result_1(ol);
    Rcpp::IntegerVector n_result_2(ol);
    Rcpp::IntegerVector n_result_3(ol);
    Rcpp::IntegerVector n_total(ol);

    Rcpp::IntegerVector n_above_cutoff(ol);
    Rcpp::IntegerVector n_below_cutoff(ol);

    Rcpp::IntegerVector n_individ_out(ol);
    Rcpp::IntegerVector scenario_int_out(ol);

    // Counters for running means with partially missing data:
    std::vector<int> imean_pre_nn(ol, 0L);
    std::vector<int> imean_post_nn(ol, 0L);

    // Create local vectors:
    std::vector<int> result_ll(ni);
    std::vector<double> n_screen_ll(ni);
    std::vector<double> n_pre_ll(ni);
    std::vector<double> n_post_ll(ni);
    std::vector<double> mean_pre_ll(ni);
    std::vector<double> mean_post_ll(ni);
    std::vector<double> imean_pre_ll(ni);
    std::vector<double> imean_post_ll(ni);
    std::vector<double> time_screen_ll(ni);
    std::vector<double> time_pre_ll(ni);
    std::vector<double> time_post_ll(ni);

    /* Not necessary as Rcpp vectors initialise with 0 anyway:
		for(int i=0; i<ol; ++i)
		{
      // consumables_cost[i] = 0.0;
      // salary_cost[i] = 0.0;
      // travel_cost[i] = 0.0;

      mean_cost[i] = 0.0;
      var_cost[i] = 0.0;
      efficacy[i] = 0.0;
      efficacy_var[i] = 0.0;

      n_screen[i] = 0.0;
      n_pre[i] = 0.0;
      n_post[i] = 0.0;
      mean_pre[i] = 0.0;
      mean_post[i] = 0.0;
      imean_pre[i] = 0.0;
      imean_post[i] = 0.0;

      n_result_0[i] = 0L;
      n_result_1[i] = 0L;
      n_result_2[i] = 0L;
      n_result_3[i] = 0L;

      n_above_cutoff[i] = 0L;
      n_below_cutoff[i] = 0L;

      imean_pre_nn[i] = 0L;
      imean_post_nn[i] = 0L;
    }
    */

    // Save output and n_individ:
    int ind=0L;
    for(int s=0L; s<ns; ++s)
    {
      for(int i=0; i<ni; ++i)
      {
        scenario_int_out[ind] = s;
        n_individ_out[ind] = n_individ[i];
        ind++;
      }
    }

    for(int p=0L; p<np; ++p)
	  {
      // If summarising then pass a references to local vectors:
      survey.run(n_individ, mu_pre[p], reduction[p], individ_cv[p], day_cv[p], aliquot_cv[p], reduction_cv[p],
            count_intercept[p], count_coefficient[p], count_add[p], count_mult[p],
            &result_ll[0L], &n_screen_ll[0L], &n_pre_ll[0L], &n_post_ll[0L],
            &mean_pre_ll[0L], &mean_post_ll[0L], &imean_pre_ll[0L], &imean_post_ll[0L],
            &time_screen_ll[0L], &time_pre_ll[0L], &time_post_ll[0L], offset);

  		// Do cost calculations and running means:
  		for(int i=0; i<ni; ++i)
  		{
        // Output index accounts for n_individ and scenario:
        const int ind = scenario_int[p] * ni + i;
        n_total[ind]++;

  			const double cons_cost = n_screen_ll[i] * cost_consumables_screen[p] +
  										n_pre_ll[i] * cost_consumables_pre[p] +
  										n_post_ll[i] * cost_consumables_post[p];

  			const double days_screen = (n_screen_ll[i] * time_consumables_screen[p] + time_screen_ll[i]) / (n_technicians[p] * 4.0 * 60.0 * 60.0);
  			const double days_pre = (n_pre_ll[i] * time_consumables_pre[p] +time_pre_ll[i]) / (n_technicians[p] * 4.0 * 60.0 * 60.0);
  			const double days_post = (n_post_ll[i] * time_consumables_post[p] + time_post_ll[i]) / (n_technicians[p] * 4.0 * 60.0 * 60.0);

        const double ndays =  (n_screen_ll[i]<0.5 ? 0.0 : ceil(days_screen)) +
                              (n_pre_ll[i]<0.5 ? 0.0 : ceil(days_pre)) +
                              (n_post_ll[i]<0.5 ? 0.0 : ceil(days_post));
  			const double sal_cost = ndays * n_team[p] * cost_salary[p];

  			const double tra_cost = ndays * cost_travel[p];

        const double tot_cost = cons_cost + sal_cost + tra_cost;

        /*
  			consumables_cost[ind] -= (consumables_cost[ind] - cons_cost) / n_total[ind];
  			salary_cost[ind] -= (salary_cost[ind] - sal_cost) / n_total[ind];
  			travel_cost[ind] -= (travel_cost[ind] - tra_cost) / n_total[ind];
        */
        const double cdelta = tot_cost - mean_cost[ind];
        mean_cost[ind] += cdelta / n_total[ind];
        var_cost[ind] += cdelta * (tot_cost - mean_cost[ind]);
        
        mean_days[ind] += (static_cast<double>(ndays) - mean_days[ind]) / n_total[ind];

        if(result_ll[i] > 3L || result_ll[i] < 0L) Rcpp::stop("Unhandled result_ll");
        n_result_0[ind] += static_cast<int>(result_ll[i] == 0L);
        n_result_1[ind] += static_cast<int>(result_ll[i] == 1L);
        n_result_2[ind] += static_cast<int>(result_ll[i] == 2L);
        n_result_3[ind] += static_cast<int>(result_ll[i] == 3L);

        // Only use these if the scenario was successful:
        if(result_ll[i] == 0L)
        {
          const int nmi = n_result_0[ind];
          const double eff = 1.0 - mean_post_ll[i] / mean_pre_ll[i];

          const int blc = static_cast<int>(eff < cutoff[p]);
          n_below_cutoff[ind] += blc;
          n_above_cutoff[ind] += (1L-blc);

      		const double edelta = eff - efficacy[ind];
          efficacy[ind] += edelta / nmi;
      		efficacy_var[ind] += edelta * (eff - efficacy[ind]);

          n_screen[ind] -= (n_screen[ind] - n_screen_ll[i]) / nmi;
          n_pre[ind] -= (n_pre[ind] - n_pre_ll[i]) / nmi;
          n_post[ind] -= (n_post[ind] - n_post_ll[i]) / nmi;

          mean_pre[ind] -= (mean_pre[ind] - mean_pre_ll[i]) / nmi;
          mean_post[ind] -= (mean_post[ind] - mean_post_ll[i]) / nmi;
        }

        // Only use these if there was a non-missing imean_pre/imean_post:
        if(!Rcpp::NumericVector::is_na(imean_pre_ll[i]))
        {
          imean_pre_nn[ind]++;
          imean_pre[ind] -= (imean_pre[ind] - imean_pre_ll[i]) / imean_pre_nn[ind];
        }
        if(!Rcpp::NumericVector::is_na(imean_post_ll[i]))
        {
          imean_post_nn[ind]++;
          imean_post[ind] -= (imean_post[ind] - imean_post_ll[i]) / imean_post_nn[ind];
        }
  		}
    }

    // Set means to NA where all elements were NA:
    for(int ind=0L; ind<ol; ++ind)
    {
      if(n_result_0[ind] == 0L)
      {
        efficacy[ind] = NA_REAL;
        efficacy_var[ind] = NA_REAL;
        mean_pre[ind] = NA_REAL;
        mean_post[ind] = NA_REAL;
      }
      else
      {
        // Also calculate variance:
        efficacy_var[ind] /= static_cast<double>(n_result_0[ind] - 1L);
      }
      if(imean_pre_nn[ind] == 0L)
      {
        imean_pre[ind] = NA_REAL;
      }
      if(imean_post_nn[ind] == 0L)
      {
        imean_post_nn[ind] = NA_REAL;
      }
      // And cost variance:
      var_cost[ind] /= static_cast<double>(n_total[ind] - 1L);
    }

  	Rcpp::DataFrame df1 = Rcpp::DataFrame::create( 	Rcpp::_["scenario_int"] = scenario_int_out, Rcpp::_["n_individ"] = n_individ_out,
                            Rcpp::_["mean_efficacy"] = efficacy, Rcpp::_["var_efficacy"] = efficacy_var,
                            Rcpp::_["n_below_cutoff"] = n_below_cutoff, Rcpp::_["n_above_cutoff"] = n_above_cutoff,
                            Rcpp::_["mean_pre"] = mean_pre, Rcpp::_["mean_post"] = mean_post,
                            Rcpp::_["imean_pre"] = imean_pre, Rcpp::_["imean_post"] = imean_post );
    Rcpp::DataFrame df2 = Rcpp::DataFrame::create( 	Rcpp::_["mean_n_screen"] = n_screen, 
                            Rcpp::_["mean_n_pre"] = n_pre,	Rcpp::_["mean_n_post"] = n_post,
                            Rcpp::_["n_result_0"] = n_result_0, Rcpp::_["n_result_1"] = n_result_1,
                            Rcpp::_["n_result_2"] = n_result_2, Rcpp::_["n_result_3"] = n_result_3,
                            Rcpp::_["n_total"] = n_total, Rcpp::_["mean_days"] = mean_days,
  													Rcpp::_["mean_cost"] = mean_cost, Rcpp::_["var_cost"] = var_cost );
                            
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("dplyr");
    Rcpp::Function bcol = pkg["bind_cols"];
    // Rcpp::Function bcol("cbind");
    df = bcol(df1, df2);
  }
  else
  {
  	const Rcpp::IntegerVector replicateID = parameters["replicateID"];

  	// Create the output vectors:
  	Rcpp::IntegerVector result(ol);
  	Rcpp::NumericVector n_screen(ol);
  	Rcpp::NumericVector n_pre(ol);
  	Rcpp::NumericVector n_post(ol);
  	Rcpp::NumericVector efficacy(ol);
  	Rcpp::NumericVector mean_pre(ol);
  	Rcpp::NumericVector mean_post(ol);
  	Rcpp::NumericVector imean_pre(ol);
  	Rcpp::NumericVector imean_post(ol);

  	Rcpp::NumericVector time_screen(ol);
  	Rcpp::NumericVector time_pre(ol);
  	Rcpp::NumericVector time_post(ol);
  	Rcpp::NumericVector time_screen_count(ol);
  	Rcpp::NumericVector time_pre_count(ol);
  	Rcpp::NumericVector time_post_count(ol);

  	Rcpp::NumericVector total_days(ol);
    
  	Rcpp::NumericVector consumables_cost(ol);
  	Rcpp::NumericVector salary_cost(ol);
  	Rcpp::NumericVector travel_cost(ol);
  	Rcpp::NumericVector total_cost(ol);

  	int ind = 0L;
    for(int p=0L; p<np; ++p)
    {
      // Otherwise pass reference to the output vector:
      survey.run(n_individ, mu_pre[p], reduction[p], individ_cv[p], day_cv[p], aliquot_cv[p], reduction_cv[p],
            count_intercept[p], count_coefficient[p], count_add[p], count_mult[p],
            &result[ind], &n_screen[ind], &n_pre[ind], &n_post[ind],
            &mean_pre[ind], &mean_post[ind], &imean_pre[ind], &imean_post[ind],
            &time_screen_count[ind], &time_pre_count[ind], &time_post_count[ind], offset);

  		// Do cost calculations:
  		for(int i=0; i<ni; ++i)
  		{
  			const double cons_cost = n_screen[ind] * cost_consumables_screen[p] +
  										n_pre[ind] * cost_consumables_pre[p] +
  										n_post[ind] * cost_consumables_post[p];

        time_screen[ind] = n_screen[ind] * time_consumables_screen[p];
        time_pre[ind] = n_pre[ind] * time_consumables_pre[p];
        time_post[ind] = n_post[ind] * time_consumables_post[p];

  			const double days_screen = (time_screen[ind] + time_screen_count[ind]) / (n_technicians[p] * 4.0 * 60.0 * 60.0);
  			const double days_pre = (time_pre[ind] + time_pre_count[ind]) / (n_technicians[p] * 4.0 * 60.0 * 60.0);
  			const double days_post = (time_post[ind] + time_post_count[ind]) / (n_technicians[p] * 4.0 * 60.0 * 60.0);

        const double ndays =  (n_screen[ind]<0.5 ? 0.0 : ceil(days_screen)) +
                              (n_pre[ind]<0.5 ? 0.0 : ceil(days_pre)) +
                              (n_post[ind]<0.5 ? 0.0 : ceil(days_post));
  			const double sal_cost = ndays * n_team[p] * cost_salary[p];

  			const double tra_cost = ndays * cost_travel[p];

        total_days[ind] = ndays;
  			consumables_cost[ind] = cons_cost;
  			salary_cost[ind] = sal_cost;
  			travel_cost[ind] = tra_cost;
  			total_cost[ind] = cons_cost + sal_cost + tra_cost;

        efficacy[ind] = 1.0 - mean_post[ind] / mean_pre[ind];

  			ind++;
  		}
    }

  	// Create the final non-summarised output:
  	Rcpp::Function expgrd("expand.grid");
  	// NB: deliberately backwards as we use expand.grid and not expand_grid:
  	Rcpp::DataFrame repmean = expgrd(n_individ, replicateID);
  	Rcpp::IntegerVector nind = repmean[0L];
  	Rcpp::IntegerVector repID = repmean[1L];

    // Note: Rcpp::DataFrame::create is limited to 20 columns
  	Rcpp::DataFrame df1 = Rcpp::DataFrame::create(
                            Rcpp::_["replicateID"] = repID, Rcpp::_["n_individ"] = nind,
  													Rcpp::_["result"] = result, Rcpp::_["efficacy"] = efficacy,
  													Rcpp::_["pre_mean"] = mean_pre, Rcpp::_["post_mean"] = mean_post,
  													Rcpp::_["pre_imean"] = imean_pre, Rcpp::_["post_imean"] = imean_post,
                            Rcpp::_["n_screen"] = n_screen, Rcpp::_["n_pre"] = n_pre,	Rcpp::_["n_post"] = n_post);

  	Rcpp::DataFrame df2 = Rcpp::DataFrame::create(
                            Rcpp::_["total_days"] = total_days, 
                            Rcpp::_["time_screen"] = time_screen, Rcpp::_["time_screen_count"] = time_screen_count,
                            Rcpp::_["time_pre"] = time_pre, Rcpp::_["time_pre_count"] = time_pre_count,
                            Rcpp::_["time_post"] = time_post, Rcpp::_["time_post_count"] = time_post_count,
  													Rcpp::_["consumables_cost"] = consumables_cost, Rcpp::_["salary_cost"] = salary_cost,
  													Rcpp::_["travel_cost"] = travel_cost, Rcpp::_["total_cost"] = total_cost);

    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("dplyr");
    Rcpp::Function bcol = pkg["bind_cols"];
    // Rcpp::Function bcol("cbind");
    df = bcol(df1, df2);
	}

	return df;
}
