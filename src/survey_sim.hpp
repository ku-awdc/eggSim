#include <Rcpp.h>

Rcpp::DataFrame survey_sim_std(const std::string& design, const Rcpp::DataFrame& parameters);

Rcpp::DataFrame survey_sim_nstd(const std::string& design, const Rcpp::DataFrame& parameters);
