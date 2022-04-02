#include <Rcpp.h>

Rcpp::DataFrame survey_sim_std(const std::string& design, const Rcpp::DataFrame& parameters, const Rcpp::IntegerVector& n_individ, const bool summarise);

Rcpp::DataFrame survey_sim_nstd(const std::string& design, const Rcpp::DataFrame& parameters, const Rcpp::IntegerVector& n_individ, const bool summarise);
