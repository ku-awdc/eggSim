#include <Rcpp.h>
#include "utilities.hpp"

double rgamma_cv(const double mu, const double cv)
{
  if(cv <= 0.0)
  {
    return mu;
  }
  const double k = pow(cv, -2.0);
  const double rv = R::rgamma(k, mu/k);
	return rv;
}

int rnbinom_cv(const double mu, const double cv)
{
  if(cv <= 0.0)
  {
    return R::rpois(mu);
  }
  const double k = pow(cv, -2.0);
  const int rv = rnbinom_mu(k, mu);    // get compiler error with rnbinom_mu for some reason
  // TODO: check for int overflow
  
	return rv;
}

double rbeta_cv(const double mu, const double cv)
{
  if(cv <= 0.0)
  {
    return mu;
  }
  
  const double sd = cv * mu;
  const double a = mu * ( (mu*(1.0-mu) / pow(sd,2.0)) - 1.0 );
  const double b = (1.0 - mu) * ( (mu*(1.0-mu) / pow(sd,2.0)) - 1.0);
  const double rv = R::rbeta(a, b);
	return rv;
}

double count_time(const double count, const double intercept, const double coefficient)
{
  // TODO: should this be log10 or loge???
	// log10(time to read in sec) = int + coef*log10(egg counts+1)^2 - these are raw egg counts (not in EPG)
	const double rv = std::pow(10.0, intercept + coefficient*std::pow(std::log10(count+1.0), 2.0));
	return rv;
}

