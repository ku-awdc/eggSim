#ifndef COUNT_SUMMARISE_H
#define COUNT_SUMMARISE_H

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

  It does not do some things including:
  - Tracking imeans i.e. mean of child means
  - Counting number of included individuals/samples

*/


struct CountParams
{
  double intercept;
  double coefficient;
  double add;
  double mult;
  int min_pos_screen;
  int min_pos_pre;
  double tail;
  double Te;
  double Tl;
};


template<methods method, bool t_use_screen, bool t_testing>
class CountSummarise;

template<bool t_use_screen, bool t_testing>
class CountSummarise<methods::delta, t_use_screen, t_testing>
{
private:
  // All relevant input parameters:
  const CountParams m_count_params;

  static const size_t m_tp = (t_use_screen ? 3L : 2L);
  static const size_t m_tscreen = (t_use_screen ? 0L : NA_REAL);
  static const size_t m_tpre = (t_use_screen ? 1L : 0L);
  static const size_t m_tpost = (t_use_screen ? 2L : 1L);

  // Total times:
  std::array<double, m_tp> m_total_time = {};  // Zero-initialise

  // Number of individuals with one or more count:
  std::array<int, m_tp> m_num_ind = {};  // Zero-initialise

  // Running stats for ALL individuals with one or more count:
  std::array<double, m_tp> m_means = {};  // Zero-initialise
  std::array<double, m_tp> m_varn = {};  // Zero-initialise

  // Number of individuals with pre- AND post-treatment counts:
  int m_num_pp = 0L;

  // Running stats for individuals with pre- AND post-treatment counts:
  std::array<double, 2L> m_means_pp = {};  // Zero-initialise
  std::array<double, 2L> m_varn_pp = {};  // Zero-initialise
  double m_covn_pp = 0.0;

  // Total number of positive individuals:
  std::array<int, m_tp> m_num_pos = {};  // Zero-initialise

  // Is this individual positive?
  std::array<bool, m_tp> m_is_pos = {};  // Zero-initialise

  // Number of counts within individual:
  std::array<int, m_tp> m_num_count = {};  // Zero-initialise

  // Mean counts within individual:
  std::array<double, m_tp> m_mean_count = {};  // Zero-initialise

  void add_time(const double count, const int time_point) noexcept
  {
    const double effcount = (count + m_count_params.add) * m_count_params.mult;
	  // log10(time to read in sec) = int + coef*log10(egg counts+1)^2 - these are raw egg counts (not in EPG)
    m_total_time[time_point] += std::pow(10.0, m_count_params.intercept + m_count_params.coefficient*std::pow(std::log10(effcount+1.0), 2.0));
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
      // Only increment number of individuals and update means/vars if any counts were registered:
      if(m_num_count[i] > 0L)
      {
  			const double delta = m_mean_count[i] - m_means[i];
  			m_means[i] += delta / static_cast<double>(++m_num_ind[i]);
  			m_varn[i] += delta * (m_mean_count[i] - m_means[i]);
      }
    }

    // Only update mean_pp/var_pp/cov_pp if any counts were registered BOTH pre- and post-tx:
    if(m_num_count[m_tpre] > 0L && m_num_count[m_tpost] > 0L)
    {
      const double num = static_cast<double>(++m_num_pp);

      const double dpre = m_mean_count[m_tpre] - m_means[m_tpre];
      const double dpost = m_mean_count[m_tpost] - m_means[m_tpost];
      m_means[m_tpre] += dpre / num;
      m_means[m_tpost] += dpost / num;

      const double dpost_after = m_mean_count[m_tpost] - m_means[m_tpost];
      m_varn[m_tpre] += dpre * (m_mean_count[m_tpre] - m_means[m_tpre]);
      m_varn[m_tpost] += dpost * dpost_after;

      m_covn_pp += dpre * dpost_after;
    }

    // Reset within-individual stuff:
    for(size_t i=0L; i<3L; ++i)
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
    m_mean_count[m_tscreen] -= (m_mean_count[m_tscreen] - dcount) / static_cast<double>(++m_num_count[m_tscreen]);
  }

  void add_count_pre(const int count) noexcept
  {
    const double dcount = static_cast<const double>(count);
    add_time(dcount, m_tpre);
    m_mean_count[m_tpre] -= (m_mean_count[m_tpre] - dcount) / static_cast<double>(++m_num_count[m_tpre]);
  }

  void add_count_post(const int count) noexcept
  {
    const double dcount = static_cast<const double>(count);
    add_time(dcount, m_tpost);
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

  Results get_result() const noexcept
  {
    // Note:  divide m_varm_pp and m_covm_pp by (num-1l)
		// rvar[i] = rmm[i] / (double) (_n - 1);


    Rcpp::warning("FIXME");
    // results enum is defined in enums.h
    return Results::success;
  }

  std::array<double, m_tp> get_means() const noexcept
  {
    return m_means;
  }

  std::array<double, m_tp> get_total_time() const noexcept
  {
    return m_total_time;
  }

  std::array<int, m_tp> get_num_ind() const noexcept
  {
    return m_num_ind;
  }

  std::array<int, m_tp> get_num_pos() const noexcept
  {
    return m_num_pos;
  }

};

#endif // COUNT_SUMMARISE_H
