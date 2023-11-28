#ifndef SURVEY_NS_H
#define SURVEY_NS_H

#include <Rcpp.h>

#include "utilities.h"
#include "enums.h"
#include "CountSummarise.h"
#include "distribution.h"

template<bool t_fixed_n, int nd1, int na1, int nd2, int na2,
        methods method, dists dist_individ, dists dist_day, dists dist_aliquot, dists dist_red>
void survey_ns(const int N_day_pre_, const int N_aliquot_pre_,
                 const int N_day_post_, const int N_aliquot_post_,
                 const Rcpp::IntegerVector& N_individ, const double mu_pre,
                 const double reduction, const double individ_cv, const double day_cv,
                 const double aliquot_cv, const double reduction_cv, const CountParams& count_params,
                 int* result, double* n_screen, double* n_pre, double* n_post,
                 double* n_pos_screen, double* n_pos_pre, double* n_pos_post,
                 double* mean_pre, double* mean_post, double* imean_pre, double* imean_post,
                 double* time_screen, double* time_pre, double* time_post, const ptrdiff_t offset)
{
  // Defined in enums.h:
  TESTING();

  // template<methods method, bool t_use_screen, bool t_paired, bool t_testing>
  CountSummarise<method, false, true, t_testing> count_summarise(count_params);
  
  const int N_day_pre = t_fixed_n ? nd1 : N_day_pre_;
  const int N_aliquot_pre = t_fixed_n ? na1 : N_aliquot_pre_;
  const int N_day_post = t_fixed_n ? nd2 : N_day_post_;
  const int N_aliquot_post = t_fixed_n ? na2 : N_aliquot_post_;

  distribution<dist_individ> rindivid(individ_cv);
  distribution<dist_day> rday(day_cv);
  distribution<dist_aliquot> raliquot(aliquot_cv);
  // For the beta distribution it is more efficient to know the mean in advance:
  distribution<dist_red> rred(reduction_cv, reduction);

  double pre_imean = 0.0;
  double post_imean = 0.0;

  static constexpr size_t total_tp = 2L;
  // static constexpr size_t screen_tp = 0L;
  static constexpr size_t pre_tp = 0L;
  static constexpr size_t post_tp = 1L;

  int outn = 0L;
  ptrdiff_t outoffset = 0L;

  for(int ind=0L; ind<N_individ[N_individ.length()-1L]; ++ind)
  {
    double mu_ind = rindivid.draw(mu_pre);
    pre_imean -= (pre_imean - mu_ind) / static_cast<double>(ind+1L);

    for(int day=0L; day<N_day_pre; ++day)
    {
      const double mu_day = rday.draw(mu_ind);
      for(int aliquot=0L; aliquot<N_aliquot_pre; ++aliquot)
      {
        // Allow test scenarios:
        if constexpr (s_testing==0L) {
          count_summarise.add_count_pre(raliquot.draw(mu_day));
        } else if constexpr (s_testing==1L) {
          count_summarise.add_count_pre(20L);
        } else if constexpr (s_testing==2L) {
          count_summarise.add_count_pre(ind+51L);
        } else {
          Rcpp::stop("Unrecognised testing setting");
        }
      }
    }

    mu_ind *= rred.draw();
    post_imean -= (post_imean - mu_ind) / static_cast<double>(ind+1L);

    for(int day=0L; day<N_day_post; ++day)
    {
      const double mu_day = rday.draw(mu_ind);
      for(int aliquot=0L; aliquot<N_aliquot_post; ++aliquot)
      {
        // Allow test scenarios:
        if constexpr (s_testing==0L) {
          count_summarise.add_count_post(raliquot.draw(mu_day));
        } else if constexpr (s_testing==1L) {

          int icount;
          if(N_aliquot_post == 1L)
          {
            icount = 1L;
          }
          // If this is SS_12 then the counts come from scenarios 1 & 2, respectively:
          else if(N_aliquot_post == 2L && aliquot == 0L)
          {
            icount = 1L;
          }
          else if(N_aliquot_post == 2L && aliquot == 1L)
          {
            icount = ind % 5L;
          }
          else
          {
            Rcpp::stop("Unsupported  test design");
          }

          count_summarise.add_count_post(icount);

        } else if constexpr (s_testing==2L) {

          int icount;
          if(N_aliquot_post == 1L)
          {
            icount = ind % 5L;
          }
          // If this is SS_12 then the counts come from scenarios 1 & 2, respectively:
          else if(N_aliquot_post == 2L && aliquot == 0L)
          {
            icount = 1L;
          }
          else if(N_aliquot_post == 2L && aliquot == 1L)
          {
            icount = ind % 5L;
          }
          else
          {
            Rcpp::stop("Unsupported  test design");
          }
          count_summarise.add_count_post(icount);

        } else {
          Rcpp::stop("Unrecognised testing setting");
        }
      }
    }

    count_summarise.next_ind();
    
    // Save output:
    if((ind+1L) == N_individ[outn])
    {
      {
        const Results out_result = count_summarise.get_result();
        *(result+outoffset) = static_cast<int>(out_result);
      }

      {
        const std::array<double, total_tp> out_n = count_summarise.get_total_obs();
        *(n_screen+outoffset) = 0.0;
        *(n_pre+outoffset) = out_n[pre_tp];
        *(n_post+outoffset) = out_n[post_tp];
      }

      {
        const std::array<double, total_tp> out_n_pos = count_summarise.get_total_pos();
        *(n_pos_screen+outoffset) = 0.0;
        *(n_pos_pre+outoffset) = out_n_pos[pre_tp];
        *(n_pos_post+outoffset) = out_n_pos[post_tp];
      }

      {
        // Note: we never have screening here
        const std::array<double, 2L> out_mean = count_summarise.get_means();
        *(mean_pre+outoffset) = out_mean[0L];
        *(mean_post+outoffset) = out_mean[1L];
      }

      {
        *(imean_pre+outoffset) = pre_imean;
        *(imean_post+outoffset) = post_imean;
      }

      {
        const std::array<double, total_tp> out_time = count_summarise.get_total_time();
        *(time_screen+outoffset) = 0.0;
        *(time_pre+outoffset) = out_time[pre_tp];
        *(time_post+outoffset) = out_time[post_tp];
      }

      outn++;
      outoffset += offset;
    }
  }

}


#endif // SURVEY_NS_H
