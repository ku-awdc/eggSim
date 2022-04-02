#include <Rcpp.h>

#include "enums.hpp"
#include "survey_ns.hpp"

template<designs design, methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class;

/*
template<designs design>
class survey_class
{
public:
  survey_class(const int day_screen, const int aliquot_screen,
                  const int day_pre, const int aliquot_pre,
                  const int day_post, const int aliquot_post);

  template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
  void run(const Rcpp::IntegerVector& N_individ, const double mu_pre,
              const double reduction, const double individ_cv, const double day_cv,
              const double aliquot_cv, const double reduction_cv,
		 const double count_intercept, const double count_coefficient,
		 const double count_add, const double count_mult,
		 double* efficacy, double* n_screen, double* n_pre, double* n_post,
		 double* time_count, ptrdiff_t offset) const;
};
*/

//template<>
//class survey_class<designs::NS>
template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
class survey_class<designs::NS, method, dist_individ, dist_day, dist_aliquot, dist_red>
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

//  template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
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

/*
template<>
class survey_class<designs::NS11>
{
private:
  const int m_day_screen = 0L;
  const int m_aliquot_screen = 0L;
  const int m_day_pre = 1L;
  const int m_aliquot_pre = 1L;
  const int m_day_post = 1L;
  const int m_aliquot_post = 1L;

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

  template<methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
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

*/