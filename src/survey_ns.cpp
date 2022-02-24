#include <Rcpp.h>

#include "survey_ns.hpp"

double survey_ns(int N_individ, int N_day_pre, int N_sample_pre,
                 int N_day_post, int N_sample_post, double mu_pre,
                 double individ_k, double day_k,
                 double sample_k, double efficacy_a, double efficacy_b)
{
  double pre_mean = 0.0;
  int pre_n = 0L;
  for(int ind=0L; ind<N_individ; ++ind)
  {
    double mu_ind = R::rgamma(individ_k, individ_k / mu_pre);
    for(int day=0L; day<N_day_pre; ++day)
    {
      double mu_day = R::rgamma(day_k, day_k / mu_ind);
      for(int sample=0L; sample<N_sample_pre; ++sample)
      {
        double mu_sample = R::rgamma(sample_k, sample_k / mu_day);
        int count = R::rpois(mu_sample);
        pre_mean -= (pre_mean - static_cast<double>(count)) / ++pre_n;
      }
    }
  }

	return pre_mean;
}
