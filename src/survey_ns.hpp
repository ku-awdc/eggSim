#include <Rcpp.h>

Rcpp::NumericVector survey_ns(const int N_individ, const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
                 const double weight, const double performance,
                 const double cost_sample, const double cost_aliquot,
                 const double individ_k, const double day_k,
                 const double aliquot_k, const double efficacy_a, const double efficacy_b);
