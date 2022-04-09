#include "survey_sim.h"

#include <Rcpp.h>

#include "enums.h"
#include "survey_template.h"

// Hacky macro function to avoid lots of copy/paste:

#define EXPAND(DESIGN, VARIANT, FIXN, D0, A0, D1, A1, D2, A2) ({ \
  if( design == (VARIANT) ){ \
    if( dist_string == "cs_ga_ga_po_be" ) \
    { \
      rv = survey_template<designs::DESIGN, FIXN, D0, A0, D1, A1, D2, A2, methods::custom, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta> \
      (all_ns, parameters, n_individ, summarise); \
    } \
    else if( dist_string == "cs_ga_ga_nb_be" ) \
    { \
      rv = survey_template<designs::DESIGN, FIXN, D0, A0, D1, A1, D2, A2, methods::custom, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta> \
        (all_ns, parameters, n_individ, summarise); \
    } \
    else \
    { \
      Rcpp::stop("Unhandled dist_string"); \
    } \
    handled = true; \
  } \
})

Rcpp::DataFrame survey_sim(const std::string& design, const std::string& dist_string,
                const Rcpp::IntegerVector& all_ns, const Rcpp::DataFrame& parameters,
								const Rcpp::IntegerVector& n_individ, const bool summarise)
{

	Rcpp::DataFrame rv;
  bool handled = false;

  EXPAND(NS, "NS", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(NS, "NS_11", true, 0L, 0L, 1L, 1L, 1L, 1L);
  EXPAND(NS, "NS_12", true, 0L, 0L, 1L, 1L, 1L, 2L);
  EXPAND(SS, "SS", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(SS, "SS_11", true, 0L, 0L, 1L, 1L, 1L, 1L);
  EXPAND(SS, "SS_12", true, 0L, 0L, 1L, 1L, 1L, 2L);
  EXPAND(SSR, "SSR", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(SSR, "SSR_11", true, 1L, 1L, 1L, 1L, 1L, 1L);
  EXPAND(SSR, "SSR_12", true, 1L, 1L, 1L, 1L, 1L, 2L);

  if(!handled)
  {
    Rcpp::stop("Unrecognised survey design");
  }

	return rv;
}
