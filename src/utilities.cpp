#include "utilities.h"

#include <Rcpp.h>

#include "distribution.h"
//#include "count_timer.h"
#include "enums.h"

Rcpp::NumericVector rgamma_cv(const int n, const double mu, const double cv)
{
  Rcpp::NumericVector rv(n);
  distribution<dists::rgamma> distn(cv, mu);
  for(int i=0L; i<n; ++i)
  {
    rv[i] = distn.draw();
  }
	return rv;
}

Rcpp::IntegerVector rnbinom_cv(const int n, const double mu, const double cv)
{
  Rcpp::IntegerVector rv(n);
  distribution<dists::rnbinom> distn(cv, mu);
  for(int i=0L; i<n; ++i)
  {
    rv[i] = distn.draw();
  }
	return rv;
}

Rcpp::NumericVector rbeta_cv(const int n, const double mu, const double cv)
{
  Rcpp::NumericVector rv(n);
  distribution<dists::rbeta> distn(cv, mu);
  for(int i=0L; i<n; ++i)
  {
    rv[i] = distn.draw();
  }
	return rv;
}

/*
double count_time(const Rcpp::IntegerVector& count, const double intercept,
                  const double coefficient, const double add, const double mult)
{
  count_timer<methods::custom> counter(intercept, coefficient, add, mult);
  for(int i=0L; i<count.length(); ++i)
  {
    counter.add_count(count[i]);
  }

	return counter.get_total();
}
*/

Rcpp::CharacterVector results_levels()
{
  //const std::vector<const std::string> reslevs = ResultsLevels();
  Rcpp::CharacterVector rv(ResultsLevels.size());
  for (size_t i=0L; i<rv.size(); ++i) {
    rv[i] = ResultsLevels[i];
  }
  return rv;
}
