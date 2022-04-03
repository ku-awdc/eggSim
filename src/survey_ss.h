#ifndef SURVEY_SS_H
#define SURVEY_SS_H

#include <Rcpp.h>

#include "utilities.h"
#include "enums.h"
#include "count_timer.h"
#include "distribution.h"

template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
void survey_ss(const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, 
				 const Rcpp::IntegerVector& N_individ, const double mu_pre,
                 const double reduction, const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
				 const double count_intercept, const double count_coefficient,
				 const double count_add, const double count_mult,
				 double* efficacy, double* n_screen, double* n_pre, double* n_post,
				 double* time_count, ptrdiff_t offset)
{

  #include "survey_ss_body.h"

}

template<int nd1, int na1, int nd2, int na2, 
        methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
void survey_ss_tt(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                 const double reduction, const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
				 const double count_intercept, const double count_coefficient,
				 const double count_add, const double count_mult,
				 double* efficacy, double* n_screen, double* n_pre, double* n_post,
				 double* time_count, ptrdiff_t offset)
{

  const int N_day_pre = nd1;
  const int N_aliquot_pre = na1;
  const int N_day_post = nd2;
  const int N_aliquot_post = na2;

  #include "survey_ss_body.h"

}

#endif // SURVEY_SS_H
