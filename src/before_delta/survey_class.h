#include <Rcpp.h>

#include "enums.h"
#include "survey_ns.h"
#include "survey_ss.h"
#include "survey_ssr.h"

template<designs design, bool t_fixed_n, int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class;


// Specialisation for NS:
template<bool t_fixed_n, int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::NS, t_fixed_n, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  static constexpr int m_day_screen = 0L;
  static constexpr int m_aliquot_screen = 0L;
  const int m_day_pre;
  const int m_aliquot_pre;
  const int m_day_post;
  const int m_aliquot_post;
  static constexpr int m_min_pos_screen = 0L;
  const int m_min_pos_pre;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post,
                  const int min_pos_screen, const int min_pos_pre) :
                  m_day_pre(day_pre), m_aliquot_pre(aliquot_pre),
                  m_day_post(day_post), m_aliquot_post(aliquot_post),
                  m_min_pos_pre(min_pos_pre)

  {
    if(day_screen != m_day_screen) Rcpp::stop("Invalid N_day_screen");
    if(aliquot_screen != m_aliquot_screen) Rcpp::stop("Invalid N_aliquot_screen");
    if(min_pos_screen != m_min_pos_screen) Rcpp::stop("Invalid min_positive_screen");
    if(min_pos_pre < 0L) Rcpp::stop("Invalid min_positive_pre < 0L");

    // Note: if constexpr() requires C++17
    if constexpr(t_fixed_n)
    {
      if(day_screen != nd0) Rcpp::stop("Invalid N_day_screen (does not match template)");
      if(aliquot_screen != na0) Rcpp::stop("Invalid N_aliquot_screen (does not match template)");
      if(day_pre != nd1) Rcpp::stop("Invalid N_day_pre (does not match template)");
      if(aliquot_pre != na1) Rcpp::stop("Invalid N_aliquot_pre (does not match template)");
      if(day_post != nd2) Rcpp::stop("Invalid N_day_post (does not match template)");
      if(aliquot_post != na2) Rcpp::stop("Invalid N_aliquot_post (does not match template)");
    }
    else
    {
      if(day_pre <= 0L) Rcpp::stop("Invalid N_day_pre <= 0L");
      if(aliquot_pre <= 0L) Rcpp::stop("Invalid N_aliquot_pre <= 0L");
      if(day_post <= 0L) Rcpp::stop("Invalid N_day_post <= 0L");
      if(aliquot_post <= 0L) Rcpp::stop("Invalid N_aliquot_post <= 0L");
    }
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
           const double reduction, const double individ_cv, const double day_cv,
           const double aliquot_cv, const double reduction_cv,
           const double count_intercept, const double count_coefficient,
           const double count_add, const double count_mult,
           int* result, double* n_screen, double* n_pre, double* n_post,
  				 double* mean_pre, double* mean_post, double* imean_pre, double* imean_post,
  				 double* time_screen, double* time_pre, double* time_post, ptrdiff_t offset) const
  {
    survey_ns<t_fixed_n, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>(
      m_day_pre, m_aliquot_pre, m_day_post, m_aliquot_post,
      m_min_pos_pre,
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			result, n_screen, n_pre, n_post,
      mean_pre, mean_post, imean_pre, imean_post,
      time_screen, time_pre, time_post, offset);
  }

};

// Specialisation for SS:
template<bool t_fixed_n, int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::SS, t_fixed_n, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  static constexpr int m_day_screen = 0L;
  static constexpr int m_aliquot_screen = 0L;
  const int m_day_pre;
  const int m_aliquot_pre;
  const int m_day_post;
  const int m_aliquot_post;
  static constexpr int m_min_pos_screen = 0L;
  const int m_min_pos_pre;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post,
                  const int min_pos_screen, const int min_pos_pre) :
                  m_day_pre(day_pre), m_aliquot_pre(aliquot_pre),
                  m_day_post(day_post), m_aliquot_post(aliquot_post),
                  m_min_pos_pre(min_pos_pre)

  {
    if(day_screen != m_day_screen) Rcpp::stop("Invalid N_day_screen");
    if(aliquot_screen != m_aliquot_screen) Rcpp::stop("Invalid N_aliquot_screen");
    if(min_pos_screen != m_min_pos_screen) Rcpp::stop("Invalid min_positive_screen");
    if(min_pos_pre < 0L) Rcpp::stop("Invalid min_positive_pre < 0L");

    // Note: if constexpr() requires C++17
    if constexpr(t_fixed_n)
    {
      if(day_screen != nd0) Rcpp::stop("Invalid N_day_screen (does not match template)");
      if(aliquot_screen != na0) Rcpp::stop("Invalid N_aliquot_screen (does not match template)");
      if(day_pre != nd1) Rcpp::stop("Invalid N_day_pre (does not match template)");
      if(aliquot_pre != na1) Rcpp::stop("Invalid N_aliquot_pre (does not match template)");
      if(day_post != nd2) Rcpp::stop("Invalid N_day_post (does not match template)");
      if(aliquot_post != na2) Rcpp::stop("Invalid N_aliquot_post (does not match template)");
    }
    else
    {
      if(day_pre <= 0L) Rcpp::stop("Invalid N_day_pre <= 0L");
      if(aliquot_pre <= 0L) Rcpp::stop("Invalid N_aliquot_pre <= 0L");
      if(day_post <= 0L) Rcpp::stop("Invalid N_day_post <= 0L");
      if(aliquot_post <= 0L) Rcpp::stop("Invalid N_aliquot_post <= 0L");
    }
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
           const double reduction, const double individ_cv, const double day_cv,
           const double aliquot_cv, const double reduction_cv,
           const double count_intercept, const double count_coefficient,
           const double count_add, const double count_mult,
           int* result, double* n_screen, double* n_pre, double* n_post,
  				 double* mean_pre, double* mean_post, double* imean_pre, double* imean_post,
  				 double* time_screen, double* time_pre, double* time_post, ptrdiff_t offset) const
  {
    survey_ss<t_fixed_n, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>(
      m_day_pre, m_aliquot_pre, m_day_post, m_aliquot_post,
      m_min_pos_pre,
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			result, n_screen, n_pre, n_post,
      mean_pre, mean_post, imean_pre, imean_post,
      time_screen, time_pre, time_post, offset);
  }

};


// Specialisation for SSR:
template<bool t_fixed_n, int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::SSR, t_fixed_n, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen;
  const int m_aliquot_screen;
  const int m_day_pre;
  const int m_aliquot_pre;
  const int m_day_post;
  const int m_aliquot_post;
  const int m_min_pos_screen;
  const int m_min_pos_pre;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post,
                  const int min_pos_screen, const int min_pos_pre) :
                  m_day_screen(day_screen), m_aliquot_screen(aliquot_screen),
                  m_day_pre(day_pre), m_aliquot_pre(aliquot_pre),
                  m_day_post(day_post), m_aliquot_post(aliquot_post),
                  m_min_pos_screen(min_pos_screen), m_min_pos_pre(min_pos_pre)

  {
    if(min_pos_screen < 0L) Rcpp::stop("Invalid min_positive_screen < 0L");
    if(min_pos_pre < 0L) Rcpp::stop("Invalid min_positive_pre < 0L");

    // Note: if constexpr() requires C++17
    if constexpr(t_fixed_n)
    {
      if(day_screen != nd0) Rcpp::stop("Invalid N_day_screen (does not match template)");
      if(aliquot_screen != na0) Rcpp::stop("Invalid N_aliquot_screen (does not match template)");
      if(day_pre != nd1) Rcpp::stop("Invalid N_day_pre (does not match template)");
      if(aliquot_pre != na1) Rcpp::stop("Invalid N_aliquot_pre (does not match template)");
      if(day_post != nd2) Rcpp::stop("Invalid N_day_post (does not match template)");
      if(aliquot_post != na2) Rcpp::stop("Invalid N_aliquot_post (does not match template)");
    }
    else
    {
      if(day_screen <= 0L) Rcpp::stop("Invalid N_day_screen <= 0L");
      if(aliquot_screen <= 0L) Rcpp::stop("Invalid N_aliquot_screen <= 0L");
      if(day_pre <= 0L) Rcpp::stop("Invalid N_day_pre <= 0L");
      if(aliquot_pre <= 0L) Rcpp::stop("Invalid N_aliquot_pre <= 0L");
      if(day_post <= 0L) Rcpp::stop("Invalid N_day_post <= 0L");
      if(aliquot_post <= 0L) Rcpp::stop("Invalid N_aliquot_post <= 0L");
    }
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
           const double reduction, const double individ_cv, const double day_cv,
           const double aliquot_cv, const double reduction_cv,
           const double count_intercept, const double count_coefficient,
           const double count_add, const double count_mult,
           int* result, double* n_screen, double* n_pre, double* n_post,
  				 double* mean_pre, double* mean_post, double* imean_pre, double* imean_post,
  				 double* time_screen, double* time_pre, double* time_post, ptrdiff_t offset) const
  {
    survey_ssr<t_fixed_n, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>(
      m_day_screen, m_aliquot_screen, m_day_pre, m_aliquot_pre, m_day_post, m_aliquot_post,
      m_min_pos_screen, m_min_pos_pre,
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			result, n_screen, n_pre, n_post,
      mean_pre, mean_post, imean_pre, imean_post,
      time_screen, time_pre, time_post, offset);
  }

};
