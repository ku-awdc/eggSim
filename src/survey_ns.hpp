#include <Rcpp.h>
using namespace Rcpp;

double survey_ns(int N_individ, int N_day_pre, int N_sample_pre,
                 int N_day_post, int N_sample_post, double mu_pre,
                 double individ_k, double day_k,
                 double sample_k, double efficacy_a, double efficacy_b);
