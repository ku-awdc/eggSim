#include "survey_sim.hpp"

#include <Rcpp.h>

#include "enums.hpp"
#include "survey_template.hpp"

Rcpp::DataFrame survey_sim(const std::string& design,
                const Rcpp::IntegerVector& all_ns, const Rcpp::DataFrame& parameters,
								const Rcpp::IntegerVector& n_individ, const bool summarise)
{
	Rcpp::DataFrame rv;
  
  // TODO: specialise on if aliquot_cv is always 0
  
	if( design == "NS_11" )
	{
    rv = survey_template<designs::NS, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rgamma>
      (all_ns, parameters, n_individ, summarise);
	}
	else if( design == "NS_12" )
	{
    rv = survey_template<designs::NS, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rgamma>
      (all_ns, parameters, n_individ, summarise);
	}
	else if( design == "NS" )
	{
    rv = survey_template<designs::NS, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rgamma>
      (all_ns, parameters, n_individ, summarise);
	}
	else
	{
		Rcpp::stop("Unrecognised survey design");
	}

	return rv;
}
