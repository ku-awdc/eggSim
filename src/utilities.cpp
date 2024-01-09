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

std::vector<std::string> GetResultsLevels()
{
  const std::vector<std::string> ResultsLevels = {
      ResultToString(Results::zero_pre),
      ResultToString(Results::few_screen),
      ResultToString(Results::few_pre),
      ResultToString(Results::efficacy_below),
      ResultToString(Results::efficacy_above),
      ResultToString(Results::class_fail),
      ResultToString(Results::resistant),
      ResultToString(Results::low_resistant),
      ResultToString(Results::inconclusive),
      ResultToString(Results::susceptible)
    };
    return ResultsLevels;
}

Rcpp::CharacterVector results_levels()
{
  //const std::vector<std::string> reslevs = ResultsLevels();
  const std::vector<std::string> ResultsLevels = GetResultsLevels();
  Rcpp::CharacterVector rv(ResultsLevels.size());
  for (uint_fast8_t i=0L; i<rv.size(); ++i) {
    rv[i] = ResultsLevels[i];
  }
  return rv;
}

Rcpp::LogicalVector results_success()
{
  const std::vector<std::string> ResultsLevels = GetResultsLevels();
  Rcpp::LogicalVector rv(ResultsLevels.size());
  for (int i=0L; i<rv.size(); ++i) {
    rv[i] = ResultIsSuccess(i);
  }
  return rv;
}
