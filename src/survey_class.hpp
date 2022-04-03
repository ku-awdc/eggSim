#include <Rcpp.h>

#include "enums.hpp"
#include "survey_ns.hpp"
#include "survey_ss.hpp"
#include "survey_ssr.hpp"

template<designs design, int nd0, int na0, int nd1, int na1, int nd2, int na2,
          methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class;


// Specialisation for general NS case:
template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::NS, -1L, -1L, -1L, -1L, -1L, -1L, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen = 0L;
  const int m_aliquot_screen = 0L;
  const int m_day_pre = 0L;
  const int m_aliquot_pre = 0L;
  const int m_day_post = 0L;
  const int m_aliquot_post = 0L;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post) :
                  m_day_screen(day_screen), m_aliquot_screen(aliquot_screen),
                  m_day_pre(day_pre), m_aliquot_pre(aliquot_pre),
                  m_day_post(day_post), m_aliquot_post(aliquot_post)

  {
    if(m_day_screen != 0L) Rcpp::stop("Invalid N_day_screen");
    if(m_aliquot_screen != 0L) Rcpp::stop("Invalid N_aliquot_screen");
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
			 const double count_intercept, const double count_coefficient,
			 const double count_add, const double count_mult,
			 double* efficacy, double* n_screen, double* n_pre, double* n_post,
			 double* time_count, ptrdiff_t offset) const
  {
    survey_ns<method, dist_individ, dist_day, dist_aliquot, dist_red>(
      m_day_pre, m_aliquot_pre, m_day_post, m_aliquot_post,
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			efficacy, n_screen, n_pre, n_post, time_count, offset);
  }

};

// Specialisation for specific NS cases:
template<int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::NS, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen = nd0;
  const int m_aliquot_screen = na0;
  const int m_day_pre = nd1;
  const int m_aliquot_pre = na1;
  const int m_day_post = nd2;
  const int m_aliquot_post = na2;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post)

  {
    if(day_screen != m_day_screen) Rcpp::stop("Invalid N_day_screen");
    if(aliquot_screen != m_aliquot_screen) Rcpp::stop("Invalid N_aliquot_screen");
    if(day_pre != m_day_pre) Rcpp::stop("Invalid N_day_pre");
    if(aliquot_pre != m_aliquot_pre) Rcpp::stop("Invalid N_aliquot_pre");
    if(day_post != m_day_post) Rcpp::stop("Invalid N_day_post");
    if(aliquot_post != m_aliquot_post) Rcpp::stop("Invalid N_aliquot_post");
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
			 const double count_intercept, const double count_coefficient,
			 const double count_add, const double count_mult,
			 double* efficacy, double* n_screen, double* n_pre, double* n_post,
			 double* time_count, ptrdiff_t offset) const
  {
    survey_ns_tt<nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>(
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			efficacy, n_screen, n_pre, n_post, time_count, offset);
  }

};

// Specialisation for general SS case:
template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::SS, -1L, -1L, -1L, -1L, -1L, -1L, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen = 0L;
  const int m_aliquot_screen = 0L;
  const int m_day_pre = 0L;
  const int m_aliquot_pre = 0L;
  const int m_day_post = 0L;
  const int m_aliquot_post = 0L;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post) :
                  m_day_screen(day_screen), m_aliquot_screen(aliquot_screen),
                  m_day_pre(day_pre), m_aliquot_pre(aliquot_pre),
                  m_day_post(day_post), m_aliquot_post(aliquot_post)

  {
    if(m_day_screen != 0L) Rcpp::stop("Invalid N_day_screen");
    if(m_aliquot_screen != 0L) Rcpp::stop("Invalid N_aliquot_screen");
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
			 const double count_intercept, const double count_coefficient,
			 const double count_add, const double count_mult,
			 double* efficacy, double* n_screen, double* n_pre, double* n_post,
			 double* time_count, ptrdiff_t offset) const
  {
    survey_ss<method, dist_individ, dist_day, dist_aliquot, dist_red>(
      m_day_pre, m_aliquot_pre, m_day_post, m_aliquot_post,
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			efficacy, n_screen, n_pre, n_post, time_count, offset);
  }

};

// Specialisation for specific SS cases:
template<int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::SS, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen = nd0;
  const int m_aliquot_screen = na0;
  const int m_day_pre = nd1;
  const int m_aliquot_pre = na1;
  const int m_day_post = nd2;
  const int m_aliquot_post = na2;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post)

  {
    if(day_screen != m_day_screen) Rcpp::stop("Invalid N_day_screen");
    if(aliquot_screen != m_aliquot_screen) Rcpp::stop("Invalid N_aliquot_screen");
    if(day_pre != m_day_pre) Rcpp::stop("Invalid N_day_pre");
    if(aliquot_pre != m_aliquot_pre) Rcpp::stop("Invalid N_aliquot_pre");
    if(day_post != m_day_post) Rcpp::stop("Invalid N_day_post");
    if(aliquot_post != m_aliquot_post) Rcpp::stop("Invalid N_aliquot_post");
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
			 const double count_intercept, const double count_coefficient,
			 const double count_add, const double count_mult,
			 double* efficacy, double* n_screen, double* n_pre, double* n_post,
			 double* time_count, ptrdiff_t offset) const
  {
    survey_ss_tt<nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>(
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			efficacy, n_screen, n_pre, n_post, time_count, offset);
  }

};

// Specialisation for general SSR case:
template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::SSR, -1L, -1L, -1L, -1L, -1L, -1L, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen = 0L;
  const int m_aliquot_screen = 0L;
  const int m_day_pre = 0L;
  const int m_aliquot_pre = 0L;
  const int m_day_post = 0L;
  const int m_aliquot_post = 0L;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post) :
                  m_day_screen(day_screen), m_aliquot_screen(aliquot_screen),
                  m_day_pre(day_pre), m_aliquot_pre(aliquot_pre),
                  m_day_post(day_post), m_aliquot_post(aliquot_post)

  {
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
			 const double count_intercept, const double count_coefficient,
			 const double count_add, const double count_mult,
			 double* efficacy, double* n_screen, double* n_pre, double* n_post,
			 double* time_count, ptrdiff_t offset) const
  {
    survey_ssr<method, dist_individ, dist_day, dist_aliquot, dist_red>(
      m_day_screen, m_aliquot_screen,
      m_day_pre, m_aliquot_pre, m_day_post, m_aliquot_post,
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			efficacy, n_screen, n_pre, n_post, time_count, offset);
  }

};

// Specialisation for specific SSR cases:
template<int nd0, int na0, int nd1, int na1, int nd2, int na2, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::SSR, nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>
{
private:
  const int m_day_screen = nd0;
  const int m_aliquot_screen = na0;
  const int m_day_pre = nd1;
  const int m_aliquot_pre = na1;
  const int m_day_post = nd2;
  const int m_aliquot_post = na2;

public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post)

  {
    if(day_screen != m_day_screen) Rcpp::stop("Invalid N_day_screen");
    if(aliquot_screen != m_aliquot_screen) Rcpp::stop("Invalid N_aliquot_screen");
    if(day_pre != m_day_pre) Rcpp::stop("Invalid N_day_pre");
    if(aliquot_pre != m_aliquot_pre) Rcpp::stop("Invalid N_aliquot_pre");
    if(day_post != m_day_post) Rcpp::stop("Invalid N_day_post");
    if(aliquot_post != m_aliquot_post) Rcpp::stop("Invalid N_aliquot_post");
  }

  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
                const double reduction, const double individ_cv, const double day_cv,
                const double aliquot_cv, const double reduction_cv,
			 const double count_intercept, const double count_coefficient,
			 const double count_add, const double count_mult,
			 double* efficacy, double* n_screen, double* n_pre, double* n_post,
			 double* time_count, ptrdiff_t offset) const
  {
    survey_ssr_tt<nd0, na0, nd1, na1, nd2, na2, method, dist_individ, dist_day, dist_aliquot, dist_red>(
      N_individ, mu_pre, reduction,
      individ_cv, day_cv, aliquot_cv, reduction_cv,
      count_intercept, count_coefficient, count_add, count_mult,
			efficacy, n_screen, n_pre, n_post, time_count, offset);
  }

};
