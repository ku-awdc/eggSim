#include <Rcpp.h>

Rcpp::NumericVector rgamma_cv(const int n, const double mu, const double cv);
Rcpp::IntegerVector rnbinom_cv(const int n, const double mu, const double cv);
Rcpp::NumericVector rbeta_cv(const int n, const double mu, const double cv);
double count_time(const Rcpp::IntegerVector& count, const double intercept, const double coefficient, const double add, const double mult);
