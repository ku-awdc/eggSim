#include <Rcpp.h>

void survey_ssr(const int N_individ, const int N_day_screen, const int N_aliquot_screen, 
                 const int N_day_pre, const int N_aliquot_pre,
                 const int N_day_post, const int N_aliquot_post, const double mu_pre,
                 const double reduction, const double weight, const double performance,
                 const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv,
				 const double count_intercept, const double count_coefficient,
				 const double count_add, const double count_mult,
				 double &efficacy, int &n_screen, int &n_pre, int &n_post,
				 double &t_count);

