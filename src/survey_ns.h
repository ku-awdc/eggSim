#ifndef SURVEY_NS_H
#define SURVEY_NS_H

#include <Rcpp.h>

#include "utilities.h"
#include "enums.h"
#include "count_timer.h"
#include "distribution.h"

template<bool t_fixed_n, int nd1, int na1, int nd2, int na2,
        methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
void survey_ns(const int N_day_pre_, const int N_aliquot_pre_,
                 const int N_day_post_, const int N_aliquot_post_,
                 const int min_pos_pre,
				 const Rcpp::IntegerVector& N_individ, const double mu_pre,
                 const double reduction, const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
				 const double count_intercept, const double count_coefficient,
				 const double count_add, const double count_mult,
				 int* result, double* efficacy, double* n_screen, double* n_pre, double* n_post,
				 double* time_count, ptrdiff_t offset)
{

  const int N_day_pre = t_fixed_n ? nd1 : N_day_pre_;
  const int N_aliquot_pre = t_fixed_n ? na1 : N_aliquot_pre_;
  const int N_day_post = t_fixed_n ? nd2 : N_day_post_;
  const int N_aliquot_post = t_fixed_n ? na2 : N_aliquot_post_;

  double pre_mean = 0.0;
  double post_mean = 0.0;
  int pre_n=0L;
  int post_n=0L;
  int npos=0L;

  count_timer<method> counter_pre(count_intercept, count_coefficient, count_add, count_mult);
  count_timer<method> counter_post(count_intercept, count_coefficient, count_add, count_mult);
  distribution<dist_individ> rindivid(individ_cv);
  distribution<dist_day> rday(day_cv);
  distribution<dist_aliquot> raliquot(aliquot_cv);
  // For the beta distribution it is more efficient to know the mean in advance:
  distribution<dist_red> rred(reduction_cv, reduction);

  int outn = 0L;
  ptrdiff_t outoffset = 0L;

  for(int ind=0L; ind<N_individ[N_individ.length()-1L]; ++ind)
  {
    bool ispos=false;
    double mu_ind = rindivid.draw(mu_pre);
    for(int day=0L; day<N_day_pre; ++day)
    {
      const double mu_day = rday.draw(mu_ind);
      for(int aliquot=0L; aliquot<N_aliquot_pre; ++aliquot)
      {
        /*
        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
        int count = rpois(mu_aliquot);
        */
        const int icount = raliquot.draw(mu_day);
        ispos = ispos || (icount > 0L);
        const double count = static_cast<double>(icount);
        pre_mean -= (pre_mean - count) / static_cast<double>(++pre_n);
        counter_pre.add_count(count);
      }
    }
    npos += static_cast<int>(ispos);

    mu_ind *= rred.draw();
    for(int day=0L; day<N_day_post; ++day)
    {
      const double mu_day = rday.draw(mu_ind);
      for(int aliquot=0L; aliquot<N_aliquot_post; ++aliquot)
      {
        /*
        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
        int count = rpois(mu_aliquot);
        */
        const double count = static_cast<double>(raliquot.draw(mu_day));
        post_mean -= (post_mean - count) / static_cast<double>(++post_n);
  		  counter_post.add_count(count);
      }
    }

    // Save output:
    if((ind+1L) == N_individ[outn])
    {
      Rcpp::stop("Handle result:  0=success, 1=failed min_pos_pre, 2=failed min_pos_screen, 3=unexpected error");
      if(npos < min_pos_pre)
      {
        *(efficacy+outoffset) = NA_REAL;

        *(n_screen+outoffset) = 0.0;
        *(n_pre+outoffset) = static_cast<double>(pre_n);
        *(n_post+outoffset) = 0.0;

        *(time_count+outoffset) = counter_pre.get_total();
      }
      else
      {
        // If zero eggs observed (safe float comparison: fewer than 0.5 eggs in total):
        if(pre_n == 0L || pre_mean < (0.5/(static_cast<double>((ind+1L)*N_day_pre*N_aliquot_pre))))
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

        *(time_count+outoffset) = counter_pre.get_total() + counter_post.get_total();
      }
      outn++;
      outoffset += offset;
    }
  }

}


#endif // SURVEY_NS_H
