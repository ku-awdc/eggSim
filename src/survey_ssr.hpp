#ifndef SURVEY_SSR_HPP
#define SURVEY_SSR_HPP

#include <Rcpp.h>

#include "utilities.hpp"
#include "enums.hpp"
#include "count_timer.hpp"
#include "distribution.hpp"

template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
void survey_ssr(const int N_day_screen, const int N_aliquot_screen,
                const int N_day_pre, const int N_aliquot_pre,
                const int N_day_post, const int N_aliquot_post, 
				        const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
                const double count_intercept, const double count_coefficient,
                const double count_add, const double count_mult,
                double* efficacy, double* n_screen, double* n_pre, double* n_post,
                double* time_count, ptrdiff_t offset)
{

  #include "survey_ssr_body.hpp"

}

template<int nd0, int na0, int nd1, int na1, int nd2, int na2, 
        methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
void survey_ssr_tt(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
                const double count_intercept, const double count_coefficient,
                const double count_add, const double count_mult,
                double* efficacy, double* n_screen, double* n_pre, double* n_post,
                double* time_count, ptrdiff_t offset)
{

  const int N_day_screen = nd0;
  const int N_aliquot_screen = na0;
  const int N_day_pre = nd1;
  const int N_aliquot_pre = na1;
  const int N_day_post = nd2;
  const int N_aliquot_post = na2;

  #include "survey_ssr_body.hpp"

}

#endif // SURVEY_SSR_HPP
