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
    return mu;
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