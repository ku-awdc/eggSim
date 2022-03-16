#include <Rcpp.h>

Rcpp::DataFrame survey_sim_std(const std::string design, const Rcpp::DataFrame parameters);

Rcpp::DataFrame survey_sim_nstd(const std::string design, 
								const int n_day_pre, const int n_aliquot_pre, 
								const int n_day_post, const int n_aliquot_post, 
								const Rcpp::DataFrame parameters);
