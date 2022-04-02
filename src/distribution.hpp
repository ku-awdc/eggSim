#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <Rcpp.h>

#include "enums.hpp"

template<dists dist>
class distribution;

template<>
class distribution<dists::rgamma> {
private:
  const double m_mu = 0.0;
  const double m_cv = 0.0;
  const double m_k = 0.0;
public:
  distribution(const double cv, const double mu) : m_mu(mu), m_cv(cv), m_k(pow(cv, -2.0))
  {

  }

  distribution(const double cv) : m_mu(1.0), m_cv(cv), m_k(pow(cv, -2.0))
  {

  }

  double draw() const
  {
    return draw(m_mu);
  }

  double draw(const double mu) const
  {
    if(m_cv <= 0.0)
    {
      return mu;
    }
    const double rv = R::rgamma(m_k, mu/m_k);
  	return rv;
  }

};

template<>
class distribution<dists::rnbinom> {
private:
  const double m_mu = 0.0;
  const double m_cv = 0.0;
  const double m_k = 0.0;
public:
  distribution(const double cv, const double mu) : m_mu(mu), m_cv(cv), m_k(pow(cv, -2.0))
  {

  }
  
  distribution(const double cv) : m_mu(0.0), m_cv(cv), m_k(pow(cv, -2.0))
  {

  }

  double draw() const
  {
    return draw(m_mu);
  }
  
  double draw(const double mu) const
  {
    if(m_cv <= 0.0)
    {
      return R::rpois(mu);
    }
    
    const double rv = rnbinom_mu(m_k, mu);    // get compiler error with rnbinom_mu for some reason
  	return rv;
  }

};

template<>
class distribution<dists::rpois> {
private:
  const double m_cv = 0.0;
public:
  distribution(const double cv) : m_cv(cv)
  {
    if(m_cv > 0.0) Rcpp::stop("Poisson distribution used where cv > 0");
  }

  // Note: must have a mu argument here
  double draw(const double mu) const
  {
    return R::rpois(mu);
  }

};

template<>
class distribution<dists::rbeta> {
private:
  const double m_mu = 0.0;
  const double m_cv = 0.0;
  double m_a = 0.0;
  double m_b = 0.0;
public:
  
  distribution(const double cv, const double mu) : m_mu(mu), m_cv(cv)
  {
    const double sd = cv * mu;
    m_a = mu * ( (mu*(1.0-mu) / pow(sd,2.0)) - 1.0 );
    m_b = (1.0 - mu) * ( (mu*(1.0-mu) / pow(sd,2.0)) - 1.0);
  }

  double draw() const
  {
    if(m_cv <= 0.0)
    {
      return m_mu;
    }
    
    const double rv = R::rbeta(m_a, m_b);
  	return rv;    
  }

};

template<>
class distribution<dists::identity> {
private:
  const double m_mu = 0.0;
  const double m_cv = 0.0;
public:
  
  distribution(const double cv, const double mu) : m_mu(mu), m_cv(cv)
  {
    if(m_cv > 0.0) Rcpp::stop("Identity distribution used where cv > 0");
  }

  distribution(const double cv) : m_mu(1.0), m_cv(cv)
  {
    if(m_cv > 0.0) Rcpp::stop("Identity distribution used where cv > 0");
  }

  double draw() const
  {
    return draw(m_mu);
  }

  double draw(const double mu) const
  {
    return mu;
  }

};

#endif // DISTRIBUTION_HPP
