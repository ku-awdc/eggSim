#ifndef COUNT_SUMMARISE_H
#define COUNT_SUMMARISE_H

// [[Rcpp::depends(bayescount)]]
//#include "bayescount/FecrtClassify.h"
#include "bayescount/bnb_pval.h"
#include "bayescount/waavp_ci.h"

#include <Rcpp.h>
#include <array>

#include "enums.h"

/*
  CountSummarise offloads some of the work from the underlying survey methods, including:
  - Tracking if an individual has a positive screen/pre/post total count
  - Calculating pre- and post-treatment (running) means and variance/covariance/whatever else is needed
  - Calculating the observed reduction
  - Generating a result classification for the desired method
  - Tracking time spent in the lab
  - Counting number of pre-treatment positive individuals
  - Counting number of included individuals/samples

  But not;
  - Tracking imeans i.e. mean of child means
*/

struct CountReturn {
  // Results enum is defined in enums.h
  Results result;
  double target_stat = NA_REAL;
  double lower_stat = NA_REAL;
};

template<methods t_method, bool t_use_screen, bool t_paired, bool t_testing, int t_post_mult>
class CountSummarise;

/*

template<bool t_use_screen, bool t_paired, bool t_testing>
class CountSummarise<methods::delta, t_use_screen, t_paired, t_testing>
{
*/

template<methods t_method, bool t_use_screen, bool t_paired, bool t_testing, int t_post_mult>
class CountSummarise
{
private:
  // All relevant input parameters:
  const CountParams m_count_params;

  static constexpr size_t m_tp = (t_use_screen ? 3L : 2L);
  static constexpr size_t m_tscreen = (t_use_screen ? 0L : -1L);
  static constexpr size_t m_tpre = (t_use_screen ? 1L : 0L);
  static constexpr size_t m_tpost = (t_use_screen ? 2L : 1L);

  // Stats for separate screen/pre/post:
  // Total times:
  std::array<double, m_tp> m_total_time = {};  // Zero-initialise
  // Number of individuals with one or more count (may differ screen/pre/post):
  std::array<int, m_tp> m_total_ind = {};  // Zero-initialise
  // Number of counts (may differ screen/pre/post):
  std::array<int, m_tp> m_total_count = {};  // Zero-initialise
  // Total number of non-zero individuals:
  std::array<int, m_tp> m_total_pos = {};  // Zero-initialise

  // Stats for pre&post (i.e. only individuals with both):
  // Number of individuals with pre- AND post-treatment counts:
  int m_num_pp = 0L;
  // Running stats for individuals with pre- AND post-treatment counts:
  std::array<double, 2L> m_means_pp = {};  // Zero-initialise
  std::array<double, 2L> m_varn_pp = {};  // Zero-initialise
  double m_covn_pp = 0.0;

  // Stats for use within individual (calculating individual means):
  // Is this individual positive?
  std::array<bool, m_tp> m_is_pos = {};  // Zero-initialise
  // Number of counts within individual:
  std::array<int, m_tp> m_num_count = {};  // Zero-initialise
  // Mean counts within individual:
  std::array<double, m_tp> m_mean_count = {};  // Zero-initialise

  // Allow extraction of int/double stats when we haven't yet hit the thresholds:
  template <typename T, size_t N>
  const std::array<double, N> apply_minimums(const std::array<T, N> input, const double replacement) const noexcept
  {
    std::array<double, N> arr;
    for (size_t i=0L; i<N; ++i) {
      arr[i] = static_cast<double>(input[i]);
    }
    if constexpr (t_use_screen && N == 3L) {   // Must be using screening
      // Remove pre and post if screening failed:
      if (m_total_pos[m_tscreen] < m_count_params.min_pos_screen)
      {
        arr[m_tpre] = replacement;
        arr[m_tpost] = replacement;
      }
      // Remove post only if pre failed:
      if (m_total_pos[m_tpre] < m_count_params.min_pos_pre)
      {
        arr[m_tpost] = replacement;
      }
    } else if constexpr (t_use_screen && N == 2L) {  // Using screening but output vector is pre/post only
      // Remove pre and post if screening failed:
      if (m_total_pos[m_tscreen] < m_count_params.min_pos_screen)
      {
        arr[0L] = replacement;
        arr[1L] = replacement;
      }
      // Remove post only if pre failed:
      if (m_total_pos[m_tpre] < m_count_params.min_pos_pre)
      {
        arr[1L] = replacement;
      }
    } else if constexpr (!t_use_screen && N == 2L) {  // Not using screening and output vector is pre/post only
      // Only remove post:
      if (m_total_pos[m_tpre] < m_count_params.min_pos_pre) {
        arr[1L] = replacement;
      }
    } else {
      Rcpp::stop("Unhandled use of template apply_minimums");
    }
    return arr;
  }

  void add_time(const double count, const size_t time_point) noexcept
  {
    if constexpr (t_testing) {
      if (time_point >= m_tp) Rcpp::stop("Invalid time_point for add_time");
    }

    const double effcount = (count + m_count_params.add) * m_count_params.mult;
	  // log10(time to read in sec) = int + coef*log10(egg counts+1)^2 - these are raw egg counts (not in EPG)
    m_total_time[time_point] += std::pow(10.0, m_count_params.intercept + m_count_params.coefficient*std::pow(std::log10(effcount+1.0), 2.0));
  }

  // Utility function for squaring (faster than pow with some compilers):
  double square(const double x) const
  {
    return x*x;
  }

  // Underlying functions for the delta methods of obtaining CI:
  std::array<double, 2L> levecke_ci() const
  {
    const double mu1 = m_means_pp[0];
    const double mu2 = m_means_pp[1];
    const double n = static_cast<double>(m_num_pp);
    const double var1 = m_varn_pp[0] / static_cast<double>(m_num_pp-1);
    const double var2 = m_varn_pp[1] / static_cast<double>(m_num_pp-1);
    const double cov12 = std::min(0.99, m_covn_pp / static_cast<double>(m_num_pp-1)); // stop perfect correlations breaking things

    std::array<double, 2L> rv;
    if constexpr (t_paired)
    {
      //  From Levecke et al:  Assessment of Anthelmintic Efficacy of Mebendazole inSchool Children in Six Countries Where Soil-TransmittedHelminths Are Endemic

      // const double cor12 = cov12 / std::sqrt(var1 * var2);
      // double varred = std::pow(mu2/mu1, 2) * (var1/std::pow(mu1,2) + var2/std::pow(mu2,2) - 2.0 * cor12 * std::sqrt(var1 * var2) / (mu1 * mu2));
    	// Equivalent:
    	const double varred = square(mu2/mu1) * (var1/square(mu1) + var2/square(mu2) - 2.0 * cov12 / (mu1 * mu2));
    	const double r = mu2/mu1;
    	const double shape = (square(r) * n) / varred;
    	const double scale = varred / (r * n);

    	// Signature of qgamma is:  qgamma(double p, double alpha, double scale, int lower_tail, int log_p)
      rv[0L] = 1.0 - R::qgamma(1.0-m_count_params.tail, shape, scale, 1L, 0L);
    	rv[1L]= 1.0 - R::qgamma(m_count_params.tail, shape, scale, 1L, 0L);

    }else{

      const double varred = square(mu2/mu1) * (var1/square(mu1) + var2/square(mu2));
    	const double r = mu2/mu1;
    	const double shape = (square(r) * n) / varred;
    	const double scale = varred / (r * n);

    	// Signature of qgamma is:  qgamma(double p, double alpha, double scale, int lower_tail, int log_p)
      rv[0L] = 1.0 - R::qgamma(1.0-m_count_params.tail, shape, scale, 1L, 0L);
    	rv[1L]= 1.0 - R::qgamma(m_count_params.tail, shape, scale, 1L, 0L);

    }

    return rv;
  }

  std::array<double, 2L> waavp_ci() const
  {
    const double mu1 = m_means_pp[0];
    const double mu2 = m_means_pp[1];
    const double n = static_cast<double>(m_num_pp);
    const double var1 = m_varn_pp[0] / static_cast<double>(m_num_pp-1);
    const double var2 = m_varn_pp[1] / static_cast<double>(m_num_pp-1);
    const double cov12 = std::min(0.99, m_covn_pp / static_cast<double>(m_num_pp-1)); // stop perfect correlations breaking things

    std::array<double, 2L> rv;
    if constexpr (t_paired)
    {
      //  std::array<double, 2> waavp_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail)
      rv = bayescount::waavp_p_ci(mu1, mu2, var1, var2, cov12, n, m_count_params.tail);

    }else{
      // std::array<double, 2> waavp_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail)
      rv = bayescount::waavp_u_ci(mu1, mu2, var1, var2, n, n, m_count_params.tail);
    }

    return rv;
  }


public:
  CountSummarise(const CountParams& count_params) noexcept
    : m_count_params(count_params)
  {
    // TODO: argument checks for debug mode
    for(auto &is_pos : m_is_pos)
    {
      is_pos = false;
    }
  }

  void next_ind() noexcept
  {
    for(size_t i=0L; i<m_tp; ++i)
    {
      if(m_num_count[i] > 0L)
      {
        // Only increment number of individuals and update means/vars if any counts were registered:
        m_total_ind[i]++;
        m_total_count[i] += m_num_count[i];
        if(m_is_pos[i])
        {
          m_total_pos[i]++;
        }

        /* Not needed unless we want to track mean of screening or mean of all pre (for SS):
  			const double delta = m_mean_count[i] - m_means[i];
  			m_means[i] += delta / static_cast<double>(m_num_ind[i]);
  			m_varn[i] += delta * (m_mean_count[i] - m_means[i]);
        */
      }
    }

    // Only update mean_pp/var_pp/cov_pp if any counts were registered BOTH pre- and post-tx:
    if(m_num_count[m_tpre] > 0L && m_num_count[m_tpost] > 0L)
    {
      const double num = static_cast<double>(++m_num_pp);

      const double dpre = m_mean_count[m_tpre] - m_means_pp[0L];
      const double dpost = m_mean_count[m_tpost] - m_means_pp[1L];
      m_means_pp[0L] += dpre / num;
      m_means_pp[1L] += dpost / num;

      const double dpost_after = m_mean_count[m_tpost] - m_means_pp[1L];
      m_varn_pp[0L] += dpre * (m_mean_count[m_tpre] - m_means_pp[0L]);
      m_varn_pp[1L] += dpost * dpost_after;

      m_covn_pp += dpre * dpost_after;
    }

    // Reset within-individual stuff:
    for(size_t i=0L; i<m_tp; ++i)
    {
      m_num_count[i] = 0L;
      m_is_pos[i] = false;
      m_mean_count[i] = 0.0;
    }
  }

  void add_count_screen(const int count)
  {
    if constexpr (!t_use_screen)
    {
      Rcpp::stop("Logic error:  CountSummarise not set to use screening");
    }

    const double dcount = static_cast<const double>(count);
    add_time(dcount, m_tscreen);
    m_is_pos[m_tscreen] = m_is_pos[m_tscreen] || (count > 0L);
    m_mean_count[m_tscreen] -= (m_mean_count[m_tscreen] - dcount) / static_cast<double>(++m_num_count[m_tscreen]);
  }

  void add_count_pre(const int count) noexcept
  {
    const double dcount = static_cast<const double>(count);
    add_time(dcount, m_tpre);
    m_is_pos[m_tpre] = m_is_pos[m_tpre] || (count > 0L);
    m_mean_count[m_tpre] -= (m_mean_count[m_tpre] - dcount) / static_cast<double>(++m_num_count[m_tpre]);
  }

  void add_count_post(const int count) noexcept
  {
    const double dcount = static_cast<const double>(count);
    add_time(dcount, m_tpost);
    m_is_pos[m_tpost] = m_is_pos[m_tpost] || (count > 0L);
    m_mean_count[m_tpost] -= (m_mean_count[m_tpost] - dcount) / static_cast<double>(++m_num_count[m_tpost]);
  }

  bool is_pos_screen() const noexcept
  {
    return m_is_pos[m_tscreen];
  }

  bool is_pos_pre() const noexcept
  {
    return m_is_pos[m_tpre];
  }

  bool is_pos_post() const noexcept
  {
    return m_is_pos[m_tpost];
  }

  CountReturn get_result() const
  {

    CountReturn rv;

    if constexpr (t_method == methods::delta) {

      if (m_total_pos[m_tpre] == 0L) {

        rv.result = Results::zero_pre;

      } else if (t_use_screen && (m_total_pos[0L] < m_count_params.min_pos_screen)) {

        rv.result = Results::few_screen;

      } else if (m_total_pos[m_tpre] < m_count_params.min_pos_pre) {

        rv.result = Results::few_pre;

      } else if (m_num_pp < 2) {

        // Variance isn't defined for either 0 or 1 individual (can happen with dropouts):
        rv.result = Results::no_post;

      } else if (m_varn_pp[0] == 0.0 || (m_total_pos[m_tpost] > 0 && m_varn_pp[1] == 0.0)) {
        
        // Pre-treatment k isn't defined and variance (both pre/post) is unreliable:
        rv.result = Results::class_fail;

      /*
      } else if (m_total_pos[m_tpost] == 0L) {

        rv.result = Results::class_fail;
      */

      } else {

//        const std::array<double, 2L> ci = levecke_ci();
        const std::array<double, 2L> ci = waavp_ci();

        rv.target_stat = ci[1];
        rv.lower_stat = ci[0];

        // If zero-mean post-treatment use BNB:
        if ( m_total_pos[m_tpost] == 0 ) {

          // inline double bnb_pval_100(double const sum1, double const N1, double const N2, double const K, double const mean_ratio, double const H0){
          const double N = static_cast<double>(m_num_pp);
          const double mu1 = m_means_pp[0];
          const double var1 = m_varn_pp[0] / static_cast<double>(m_num_pp-1);
          const double sum1 = mu1*N;

          // Cheat a bit and make sure k is not more than 10:
          const double vareff = std::max(var1-mu1, mu1*mu1*0.1);
          // Calculate pre-treatment K:
          const double K = (mu1 * mu1) / (vareff);

          // This only works for the specifically templated designs:
          if constexpr (t_post_mult==0)
          {
            Rcpp::stop("BNB method for non-templated designs (unknown mean ratio) needs fixing");
          }
          constexpr double mean_ratio = static_cast<double>(t_post_mult);
          double const H0 = m_count_params.Tlow;

          double const pval = bayescount::bnb_pval_100(sum1, N, N, K, mean_ratio, H0);

          if( !R_finite(pval) || std::abs(sum1 - (mu1*N)) > 0.001){
            std::string msg = std::format("Error with following parameters: sum1={0}, m_num_pp={1}, m_means_pp[0]={2}, m_varn_pp[0]={3}, vareff={4} -> pval={5}", sum1, m_num_pp, m_means_pp[0], m_varn_pp[0], vareff, pval);
            Rcpp::Rcout << msg << "\n";
            Rcpp::warning(msg);
          }

          // NB: only the lower p-value is relevant as we are only using this for 100% observed reduction!
          rv.target_stat = NA_REAL;
          rv.lower_stat = pval;

          if ( Rcpp::NumericVector::is_na(pval) ) {
            rv.result = Results::class_fail;
          } else if (pval < m_count_params.tail) {
            rv.result = Results::susceptible;
          } else {
            rv.result = Results::inconclusive;
          }

        } else if ( Rcpp::NumericVector::is_na(ci[0L]) || Rcpp::NumericVector::is_na(ci[1L]) ) {
          rv.result = Results::class_fail;  // Should now only happen with zero pre-treatment mean??
        } else if (ci[1L] < m_count_params.Teff && ci[0L] < m_count_params.Tlow) {
          rv.result = Results::resistant;
        } else if (ci[1L] < m_count_params.Teff && ci[0L] >= m_count_params.Tlow) {
          rv.result = Results::low_resistant;
        } else if (ci[1L] >= m_count_params.Teff && ci[0L] < m_count_params.Tlow) {
          rv.result = Results::inconclusive;
        } else if (ci[1L] >= m_count_params.Teff && ci[0L] >= m_count_params.Tlow) {
          rv.result = Results::susceptible;
        } else {
          Rcpp::Rcout << ci[0L] << "-" << ci[1L] << "(" << m_count_params.Teff << ", " << m_count_params.Tlow << ")" << std::endl;
          Rcpp::stop("Unhandled result in CountSummarise");
        }
      }

    } else if constexpr (t_method == methods::mean) {

      rv.target_stat = NA_REAL;
      rv.lower_stat = NA_REAL;

      if (m_total_pos[m_tpre] == 0L) {

        rv.result = Results::zero_pre;

      } else if (t_use_screen && (m_total_pos[0L] < m_count_params.min_pos_screen)) {

        rv.result = Results::few_screen;

      } else if (m_total_pos[m_tpre] < m_count_params.min_pos_pre) {

        rv.result = Results::few_pre;

      } else {

        const std::array<double, 2L> mus = get_means();
        const double eff = 1.0 - (mus[1L]/mus[0L]);

        if (eff < m_count_params.Tlow) {
          rv.result = Results::efficacy_below;
        } else {
          rv.result = Results::efficacy_above;
        }
      }

    } else {
      Rcpp::stop("Unhandled method in CountSummarise");
    }

    return rv;
  }

  // total time and total obs are used for costs, so replacement is 0
  std::array<double, m_tp> get_total_time() const noexcept
  {
    const std::array<double, m_tp> rv = apply_minimums<double, m_tp>(m_total_time, 0.0);
    return rv;
  }

  std::array<double, m_tp> get_total_obs() const noexcept
  {
    const std::array<double, m_tp> rv = apply_minimums<int, m_tp>(m_total_count, 0.0);
    return rv;
  }

  // means, total ind and total pos are not used for costs, so make replacement NA
  // so that summaries are for successful surveys only
  std::array<double, m_tp> get_total_ind() const noexcept
  {
    const std::array<double, m_tp> rv = apply_minimums<int, m_tp>(m_total_ind, NA_REAL);
    return rv;
  }

  std::array<double, m_tp> get_total_pos() const noexcept
  {
    const std::array<double, m_tp> rv = apply_minimums<int, m_tp>(m_total_pos, NA_REAL);
    return rv;
  }

  std::array<double, 2L> get_means() const noexcept
  {
    if (m_num_pp > 0L) {
      const std::array<double, 2L> rv = apply_minimums<double, 2L>(m_means_pp, NA_REAL);
      return rv;
    } else {
      // Note: can't be static or constexpr because NA_REAL isn't
      const std::array<double, 2L> rv = { NA_REAL, NA_REAL };
      return rv;
    }
  }


};

#endif // COUNT_SUMMARISE_H
