#include "survey_sim.h"

#include <Rcpp.h>

#include "enums.h"
#include "survey_template.h"

// Hacky macro function to avoid lots of copy/paste:

#define EXPAND(SUMMARISE, DESIGN, VARIANT, FIXN, D0, A0, D1, A1, D2, A2) ({ \
  if( summarise == SUMMARISE && design == (VARIANT) ){ \
    if( dist_string == "cs_ga_ga_po_be_mean" ) \
    { \
      rv = survey_template<SUMMARISE, designs::DESIGN, FIXN, D0, A0, D1, A1, D2, A2, methods::mean, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta> \
      (all_ns, parameters, count_parameters, n_individ, summarise); \
    } \
    else if( dist_string == "cs_ga_ga_nb_be_mean" ) \
    { \
      rv = survey_template<SUMMARISE, designs::DESIGN, FIXN, D0, A0, D1, A1, D2, A2, methods::mean, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta> \
        (all_ns, parameters, count_parameters, n_individ, summarise); \
    } \
    else if( dist_string == "cs_ga_ga_po_be_delta" ) \
    { \
      rv = survey_template<SUMMARISE, designs::DESIGN, FIXN, D0, A0, D1, A1, D2, A2, methods::delta, dists::rgamma, dists::rgamma, dists::rpois, dists::rbeta> \
      (all_ns, parameters, count_parameters, n_individ, summarise); \
    } \
    else if( dist_string == "cs_ga_ga_nb_be_delta" ) \
    { \
      rv = survey_template<SUMMARISE, designs::DESIGN, FIXN, D0, A0, D1, A1, D2, A2, methods::delta, dists::rgamma, dists::rgamma, dists::rnbinom, dists::rbeta> \
        (all_ns, parameters, count_parameters, n_individ, summarise); \
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
                const Rcpp::DataFrame& count_parameters,
								const Rcpp::IntegerVector& n_individ, const bool summarise)
{

  Rcpp::warning("Allow methods both delta and means");

	Rcpp::DataFrame rv;
  bool handled = false;

  EXPAND(false, NS, "NS", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(false, NS, "NS_11", true, 0L, 0L, 1L, 1L, 1L, 1L);
  EXPAND(false, NS, "NS_12", true, 0L, 0L, 1L, 1L, 1L, 2L);
  EXPAND(true, NS, "NS", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(true, NS, "NS_11", true, 0L, 0L, 1L, 1L, 1L, 1L);
  EXPAND(true, NS, "NS_12", true, 0L, 0L, 1L, 1L, 1L, 2L);

  // NB: the definition of SS_12 is different to that of NS_12 and SSR_12!
  EXPAND(false, SS, "SS", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(false, SS, "SS_11", true, 0L, 0L, 1L, 1L, 1L, 1L);
  EXPAND(false, SS, "SS_12", true, 0L, 0L, 1L, 2L, 1L, 2L);
  EXPAND(true, SS, "SS", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(true, SS, "SS_11", true, 0L, 0L, 1L, 1L, 1L, 1L);
  EXPAND(true, SS, "SS_12", true, 0L, 0L, 1L, 2L, 1L, 2L);

  EXPAND(false, SSR, "SSR", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(false, SSR, "SSR_11", true, 1L, 1L, 1L, 1L, 1L, 1L);
  EXPAND(false, SSR, "SSR_12", true, 1L, 1L, 1L, 1L, 1L, 2L);
  EXPAND(true, SSR, "SSR", false, -1L, -1L, -1L, -1L, -1L, -1L);
  EXPAND(true, SSR, "SSR_11", true, 1L, 1L, 1L, 1L, 1L, 1L);
  EXPAND(true, SSR, "SSR_12", true, 1L, 1L, 1L, 1L, 1L, 2L);

  if(!handled)
  {
    Rcpp::stop("Unrecognised survey design");
  }

	return rv;
}
