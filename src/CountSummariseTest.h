#include <Rcpp.h>

#include "enums.h"
#include "CountSummarise.h"


// A thin wrapper for a testing interface via Rcpp

class CountSummariseTest
{

  static const size_t m_tp = 3L;

  // template<methods method, bool t_use_screen, bool t_testing>
  CountSummarise<methods::delta, m_tp==3L ? true : false, true> m_cs;

  // CountSummariseTest() = delete;
  
  Rcpp::NumericVector arr2rcpp(std::array<double, m_tp> arr) const
  {
    Rcpp::NumericVector rv(m_tp);
    for(size_t i=0; i<m_tp; ++i)
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
          50L,
          50L,
          0.025, // double tail;
          0.7, // double Te;
          0.5  // double Tl;
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

  Rcpp::NumericVector get_means() const
  {
    return arr2rcpp(m_cs.get_means());
  }
  
  Rcpp::CharacterVector get_result() const
  {
    return ResultToString(m_cs.get_result());
  }

  void show() const
  {
    Rcpp::Rcout << "hello" << std::endl;
  }

};
