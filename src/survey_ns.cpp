#include "survey_ns.hpp"

#include <Rcpp.h>
#include "utilities.hpp"

void survey_ns(Rcpp::IntegerVector N_individ, const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
                 const double reduction, const double weight, const double performance,
                 const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
				 const double count_intercept, const double count_coefficient,
				 const double count_add, const double count_mult,
				 double* efficacy, double* n_screen, double* n_pre, double* n_post,
				 double* time_count, ptrdiff_t offset)
{
  const double wp = weight * performance;
  
  double pre_mean = 0.0;
  double post_mean = 0.0;
  int pre_n=0L;
  int post_n=0L;

  double t_count = 0.0;
  
  int outn = 0L;
  ptrdiff_t outoffset = 0L;

  for(int ind=0L; ind<N_individ[N_individ.length()-1L]; ++ind)
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
    
    // Save output:
    if((ind+1L) == N_individ[outn])
    {
      // If zero eggs observed (safe float comparison: fewer than 0.5 eggs in total):
      if(pre_mean < (0.5/(static_cast<double>((ind+1L)*N_day_pre*N_aliquot_pre))))
      {
        *(efficacy+outoffset) = NA_REAL;
      }
      else
      {
        *(efficacy+outoffset) = 1.0 - post_mean/pre_mean;
      }

      *(n_screen+outoffset) = 0.0;
      *(n_pre+outoffset) = static_cast<double>(pre_n);
      *(n_post+outoffset) = static_cast<double>(post_n);

      *(time_count+outoffset) = t_count;

      outn++;
      outoffset += offset;
    }
  }

}
