#include "survey_sim.hpp"

#include <Rcpp.h>

#include "enums.hpp"
#include "survey_template.hpp"

Rcpp::DataFrame survey_sim(const std::string& design, const Rcpp::StringVector& all_dists,
                const Rcpp::IntegerVector& all_ns, const Rcpp::DataFrame& parameters,
								const Rcpp::IntegerVector& n_individ, const bool summarise)
{

  // TODO: redo this as a string e.g. "ga_ga_nb_be" and template up available options
  const Rcpp::String aliquot_dist = all_dists[2L];
  
  // TODO: a pre-processor macro to generate the template set for each design
  // or alternatively, how much performance is actually lost by having run-time polymorphism on survey type within survey_template???

	Rcpp::DataFrame rv;

	if( design == "NS_12" )
	{
    if( aliquot_dist == "rpois" )
    {
      rv = survey_template<designs::NS12, methods::custom, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta>
        (all_ns, parameters, n_individ, summarise);
    }
    else if( aliquot_dist == "rnbinom" )
    {
      rv = survey_template<designs::NS12, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta>
        (all_ns, parameters, n_individ, summarise);
    }
    else
    {
      Rcpp::stop("Unrecognised aliquot_dist");
    }
	}

	else if( design == "NS_11" || design == "NS" )
	{
	  if( aliquot_dist == "rpois" )
	  {
	    rv = survey_template<designs::NS, methods::custom, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta>
	    (all_ns, parameters, n_individ, summarise);
	  }
	  else if( aliquot_dist == "rnbinom" )
	  {
	    rv = survey_template<designs::NS, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta>
	    (all_ns, parameters, n_individ, summarise);
	  }
	  else
	  {
	    Rcpp::stop("Unrecognised aliquot_dist");
	  }
	}
	else if( design == "SS" || design == "SS_11" || design == "SS_12" )
	{
	  if( aliquot_dist == "rpois" )
	  {
	    rv = survey_template<designs::SS, methods::custom, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta>
	    (all_ns, parameters, n_individ, summarise);
	  }
	  else if( aliquot_dist == "rnbinom" )
	  {
	    rv = survey_template<designs::SS, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta>
	    (all_ns, parameters, n_individ, summarise);
	  }
	  else
	  {
	    Rcpp::stop("Unrecognised aliquot_dist");
	  }
	}
	else if( design == "SSR" || design == "SSR_11" || design == "SSR_12" )
	{
	  if( aliquot_dist == "rpois" )
	  {
	    rv = survey_template<designs::SSR, methods::custom, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta>
	    (all_ns, parameters, n_individ, summarise);
	  }
	  else if( aliquot_dist == "rnbinom" )
	  {
	    rv = survey_template<designs::SSR, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta>
	    (all_ns, parameters, n_individ, summarise);
	  }
	  else
	  {
	    Rcpp::stop("Unrecognised aliquot_dist");
	  }
	}
	else
	{
		Rcpp::stop("Unrecognised survey design");
	}

	return rv;
}
