#include <Rcpp.h>

Rcpp::DataFrame survey_sim(const std::string& design, const Rcpp::StringVector& all_dists, const Rcpp::IntegerVector& all_ns, 
                    const Rcpp::DataFrame& parameters, const Rcpp::IntegerVector& n_individ, const bool summarise);
