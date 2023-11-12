#ifndef COUNT_TIMER_H
#define COUNT_TIMER_H

#include <Rcpp.h>

#include "enums.h"

template<methods method>
class count_timer;

template<>
class count_timer<methods::custom>
{
private:
  double m_total = 0.0;
  const double m_intercept = 0.0;
  const double m_coefficient = 0.0;
  const double m_add = 0.0;
  const double m_mult = 0.0;
public:
  count_timer(const double intercept, const double coefficient, const double add, const double mult)
    : m_intercept(intercept), m_coefficient(coefficient), m_add(add), m_mult(mult)
  {

  }

  void add_count(const double count)
  {
    const double effcount = (count+m_add)*m_mult;
	  // log10(time to read in sec) = int + coef*log10(egg counts+1)^2 - these are raw egg counts (not in EPG)
    m_total += std::pow(10.0, m_intercept + m_coefficient*std::pow(std::log10(effcount+1.0), 2.0));
  }

  double get_total() const
  {
    return m_total;
  }

};

#endif // COUNT_TIMER_H
