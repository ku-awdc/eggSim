#include <Rcpp.h>

Rcpp::DataFrame survey_sim(const std::string& design, const std::string& dist_string, const Rcpp::IntegerVector& all_ns, 
                    const Rcpp::DataFrame& parameters, const Rcpp::DataFrame& count_parameters,
                    const Rcpp::IntegerVector& n_individ, const bool summarise);
