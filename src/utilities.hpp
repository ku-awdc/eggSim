#include <Rcpp.h>

// TODO: local compilation within e.g. survey_ns for optimisation?

double rgamma_cv(const double mu, const double cv);
int rnbinom_cv(const double mu, const double cv);
double rbeta_cv(const double mu, const double cv);
double count_time(const double count, const double intercept, const double coefficient);
