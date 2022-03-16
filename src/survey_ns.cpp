#include "survey_ns.hpp"

#include <Rcpp.h>
#include "utilities.hpp"

Rcpp::NumericVector survey_ns(const int N_individ, const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
                 const double reduction, const double weight, const double performance,
                 const double cost_sample, const double cost_aliquot,
                 const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv)
{
  double pre_mean = 0.0;
  double post_mean = 0.0;
  int pre_n=0L;
  int post_n=0L;

  const double wp = weight * performance;

  for(int ind=0L; ind<N_individ; ++ind)
  {
    double mu_ind = rgamma_cv(mu_pre, individ_cv);
    for(int day=0L; day<N_day_pre; ++day)
    {
      const double mu_day = rgamma_cv(mu_ind, day_cv) * wp;
      for(int aliquot=0L; aliquot<N_aliquot_pre; ++aliquot)
      {
        /*
        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
        int count = rpois(mu_aliquot);
        */
        const int count = rnbinom_cv(mu_day, aliquot_cv);
        pre_mean -= (pre_mean - static_cast<double>(count)) / static_cast<double>(++pre_n);
      }
    }

    mu_ind *= (1.0 - rbeta_cv(reduction, reduction_cv));
    for(int day=0L; day<N_day_post; ++day)
    {
      const double mu_day = rgamma_cv(mu_ind, day_cv) * wp;
      for(int aliquot=0L; aliquot<N_aliquot_post; ++aliquot)
      {
        /*
        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
        int count = rpois(mu_aliquot);
        */
        const int count = rnbinom_cv(mu_day, aliquot_cv);
        post_mean -= (post_mean - static_cast<double>(count)) / static_cast<double>(++post_n);
      }
    }
  }

  Rcpp::NumericVector rv(2L);

  // If zero eggs observed (safe float comparison: fewer than 0.5 eggs in total):
  if(pre_mean < (0.5/(N_individ*N_day_pre*N_aliquot_pre)))
  {
    rv[0] = NA_REAL;
  }
  else
  {
    rv[0] = 1.0 - post_mean/pre_mean;
  }

  // Total cost:
  const double cost = N_individ * N_day_pre * (cost_sample + cost_aliquot * N_aliquot_pre) +
                N_individ * N_day_post * (cost_sample + cost_aliquot * N_aliquot_post);
  rv[1] = cost;

  return rv;
}
