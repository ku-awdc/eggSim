#include <Rcpp.h>

#include "enums.h"
#include "CountSummarise.h"


// A thin wrapper for a testing interface via Rcpp

class CountSummariseTest
{

  // template<methods method, bool t_use_screen, bool t_paired, bool t_testing>
  CountSummarise<methods::delta, false, true, true> m_cs;

  // CountSummariseTest() = delete;

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
  CountSummariseTest() // Rcpp::List count_params)
    : m_cs(
        CountParams {
          2.38, // count_intercept,
          0.066, // count_coefficient,
          0.0, // count_add,
          1.0, // count_mult,
          50L, // min_pos_screen,
          50L, // min_pos_pre,
          0.025, // double tail;
          0.7, // double Teff;
          0.5  // double Tlow;
        }
      )
  {

  }

  void add_counts(const Rcpp::IntegerVector& pre, const Rcpp::IntegerVector& post)
  {
    const size_t n = std::max(pre.size(), post.size());

    for(size_t i=0; i<n; ++i)
    {
      if(pre.size() > i)
      {
        m_cs.add_count_pre(pre[i]);
      }

      if(post.size() > i)
      {
        m_cs.add_count_post(post[i]);
      }

      m_cs.next_ind();
    }
  }

  Rcpp::List get_result() const
  {
    Rcpp::CharacterVector result = ResultToString(m_cs.get_result());
    Rcpp::NumericVector means = arr2rcpp(m_cs.get_means());
    Rcpp::NumericVector total_time = arr2rcpp(m_cs.get_total_time());
    Rcpp::NumericVector total_obs = arr2rcpp(m_cs.get_total_obs());
    Rcpp::NumericVector total_ind = arr2rcpp(m_cs.get_total_ind());
    Rcpp::NumericVector total_pos = arr2rcpp(m_cs.get_total_pos());

    Rcpp::List rv = Rcpp::List::create(
      Rcpp::Named("result") = result,
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
    Rcpp::Rcout << ResultToString(m_cs.get_result()) << std::endl;
  }

};
