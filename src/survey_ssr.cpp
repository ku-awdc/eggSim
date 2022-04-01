#include "survey_ssr.hpp"

#include <Rcpp.h>
#include "utilities.hpp"

// TODO: meanepg_weight_performance should be pre-computed

void survey_ssr(const int N_individ, const int N_day_screen, const int N_aliquot_screen,
                 const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
                 const double reduction, const double weight, const double performance,
                 const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
                 const double count_intercept, const double count_coefficient,
                 const double count_add, const double count_mult,
                 double &efficacy, int &n_screen, int &n_pre, int &n_post, double &t_count)
{
  double pre_mean = 0.0;
  double post_mean = 0.0;
  int screen_n=0L;
  int pre_n=0L;
  int post_n=0L;

  const double wp = weight * performance;
  t_count = 0.0;

  for(int ind=0L; ind<N_individ; ++ind)
  {
    double mu_ind = rgamma_cv(mu_pre, individ_cv);

    bool included = false;
    for(int day=0L; day<N_day_screen; ++day)
    {
      const double mu_day = rgamma_cv(mu_ind, day_cv) * wp;
      for(int aliquot=0L; aliquot<N_aliquot_screen; ++aliquot)
      {
				const int count = rnbinom_cv(mu_day, aliquot_cv);
        t_count += count_time((static_cast<double>(count)+count_add)*count_mult, count_intercept, count_coefficient);
				included = included || count > 0L;
        screen_n++;
      }
    }

		if(included){
	    for(int day=0L; day<N_day_pre; ++day)
	    {
	      const double mu_day = rgamma_cv(mu_ind, day_cv) * wp;
	      for(int aliquot=0L; aliquot<N_aliquot_pre; ++aliquot)
	      {
	        /*
	        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
	        int count = rpois(mu_aliquot);
	        */
	        const double count = static_cast<double>(rnbinom_cv(mu_day, aliquot_cv));
	        pre_mean -= (pre_mean - count) / static_cast<double>(++pre_n);
			    t_count += count_time((count+count_add)*count_mult, count_intercept, count_coefficient);
	      }
	    }

	    mu_ind *= rbeta_cv(reduction, reduction_cv);
	    for(int day=0L; day<N_day_post; ++day)
	    {
	      const double mu_day = rgamma_cv(mu_ind, day_cv) * wp;
	      for(int aliquot=0L; aliquot<N_aliquot_post; ++aliquot)
	      {
	        /*
	        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
	        int count = rpois(mu_aliquot);
	        */
	        const double count = static_cast<double>(rnbinom_cv(mu_day, aliquot_cv));
	        post_mean -= (post_mean - count) / static_cast<double>(++post_n);
	  		  t_count += count_time((count+count_add)*count_mult, count_intercept, count_coefficient);
	      }
	    }
    }
  }

  // If zero eggs observed (safe float comparison: fewer than 0.5 eggs in total):
  if(pre_mean < (0.5/(N_individ*N_day_pre*N_aliquot_pre)))
  {
    efficacy = NA_REAL;
  }
  else
  {
    efficacy = 1.0 - post_mean/pre_mean;
  }

  n_screen = screen_n;
  n_pre = pre_n;
  n_post = post_n;

}
