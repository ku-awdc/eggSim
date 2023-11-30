#include <Rcpp.h>

#include "enums.h"
#include "CountSummarise.h"


// A thin wrapper for a testing interface via Rcpp

class CountSummariseTest
{

  const CountParams m_cp;

  // template<methods method, bool t_use_screen, bool t_paired, bool t_testing>
  CountSummarise<methods::delta, false, true, true> m_cs_delta;
  CountSummarise<methods::mean, false, true, true> m_cs_mean;

  CountSummariseTest() = delete;

  Rcpp::NumericVector arr2rcpp(std::array<double, 2L> arr) const
  {
    Rcpp::NumericVector rv(arr.size());
    for(size_t i=0; i<arr.size(); ++i)
    {
      rv[i] = arr[i];
    }
    return rv;
  }

  Rcpp::IntegerVector arr2rcpp(std::array<int, 2L> arr) const
  {
    Rcpp::IntegerVector rv(arr.size());
    for(size_t i=0; i<arr.size(); ++i)
    {
      rv[i] = arr[i];
    }
    return rv;
  }

public:
  CountSummariseTest(const double Teff, const double Tlow, const double tail)
    : m_cp(
        CountParams {
          2.38, // count_intercept,
          0.066, // count_coefficient,
          0.0, // count_add,
          1.0, // count_mult,
          0L, // min_pos_screen,
          0L, // min_pos_pre,
          tail, // double tail;
          Teff, // double Teff;
          Tlow  // double Tlow;
        }
    ), m_cs_delta(m_cp), m_cs_mean(m_cp)
  {

  }

  void add_counts(const Rcpp::IntegerVector& pre, const Rcpp::IntegerVector& post)
  {
    const size_t n = std::max(pre.size(), post.size());

    for(size_t i=0; i<n; ++i)
    {
      if(pre.size() > i)
      {
        m_cs_delta.add_count_pre(pre[i]);
        m_cs_mean.add_count_pre(pre[i]);
      }

      if(post.size() > i)
      {
        m_cs_delta.add_count_post(post[i]);
        m_cs_mean.add_count_post(post[i]);
      }

      m_cs_delta.next_ind();
      m_cs_mean.next_ind();
    }
  }

  Rcpp::List get_result_delta() const
  {
    const CountReturn crt = m_cs_delta.get_result();
    Rcpp::CharacterVector result = ResultToString(crt.result);
    Rcpp::NumericVector means = arr2rcpp(m_cs_delta.get_means());
    Rcpp::NumericVector total_time = arr2rcpp(m_cs_delta.get_total_time());
    Rcpp::NumericVector total_obs = arr2rcpp(m_cs_delta.get_total_obs());
    Rcpp::NumericVector total_ind = arr2rcpp(m_cs_delta.get_total_ind());
    Rcpp::NumericVector total_pos = arr2rcpp(m_cs_delta.get_total_pos());

    Rcpp::List rv = Rcpp::List::create(
      Rcpp::Named("result") = result,
      Rcpp::Named("upper_ci") = crt.target_stat,
      Rcpp::Named("lower_ci") = crt.lower_stat,
      Rcpp::Named("means") = means,
      Rcpp::Named("total_time") = total_time,
      Rcpp::Named("total_observations") = total_obs,
      Rcpp::Named("total_individuals") = total_ind,
      Rcpp::Named("total_positive") = total_pos
    );

    return rv;
  }

  void show() const
  {
    Rcpp::Rcout << "Delta:  " << ResultToString(m_cs_delta.get_result().result) << ",  Mean:  " << ResultToString(m_cs_mean.get_result().result) << std::endl;
  }

};
