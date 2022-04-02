#include <Rcpp.h>

void survey_ns(Rcpp::IntegerVector N_individ, const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
                 const double reduction, const double weight, const double performance,
                 const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
				 const double count_intercept, const double count_coefficient,
				 const double count_add, const double count_mult,
				 double* efficacy, double* n_screen, double* n_pre, double* n_post,
				 double* time_count, ptrdiff_t offset);

