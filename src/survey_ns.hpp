#include <Rcpp.h>

Rcpp::NumericVector survey_ns(const int N_individ, const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
				 const double reduction, const double weight, const double performance,
                 const double cost_sample, const double cost_aliquot,
                 const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv);

