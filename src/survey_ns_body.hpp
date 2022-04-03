// This file is a dirty macro-based hack to avoid writing the body of this function twice in survey_ns.hpp

  double pre_mean = 0.0;
  double post_mean = 0.0;
  int pre_n=0L;
  int post_n=0L;

  count_timer<method> counter(count_intercept, count_coefficient, count_add, count_mult);
  distribution<dist_individ> rindivid(individ_cv);
  distribution<dist_day> rday(day_cv);
  distribution<dist_aliquot> raliquot(aliquot_cv);
  // For the beta distribution it is more efficient to know the mean in advance:
  distribution<dist_red> rred(reduction_cv, reduction);

  int outn = 0L;
  ptrdiff_t outoffset = 0L;

  for(int ind=0L; ind<N_individ[N_individ.length()-1L]; ++ind)
  {
    double mu_ind = rindivid.draw(mu_pre);
    for(int day=0L; day<N_day_pre; ++day)
    {
      const double mu_day = rday.draw(mu_ind);
      for(int aliquot=0L; aliquot<N_aliquot_pre; ++aliquot)
      {
        /*
        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
        int count = rpois(mu_aliquot);
        */
        const double count = static_cast<double>(raliquot.draw(mu_day));
        pre_mean -= (pre_mean - count) / static_cast<double>(++pre_n);
        counter.add_count(count);
      }
    }

    mu_ind *= rred.draw();
    for(int day=0L; day<N_day_post; ++day)
    {
      const double mu_day = rday.draw(mu_ind);
      for(int aliquot=0L; aliquot<N_aliquot_post; ++aliquot)
      {
        /*
        double mu_aliquot = rgamma_cv(mu_day, aliquot_cv);
        int count = rpois(mu_aliquot);
        */
        const double count = static_cast<double>(raliquot.draw(mu_day));
        post_mean -= (post_mean - count) / static_cast<double>(++post_n);
  		counter.add_count(count);
      }
    }

    // Save output:
    if((ind+1L) == N_individ[outn])
    {
      // If zero eggs observed (safe float comparison: fewer than 0.5 eggs in total):
      if(pre_n == 0L || pre_mean < (0.5/(static_cast<double>((ind+1L)*N_day_pre*N_aliquot_pre))))
      {
        *(efficacy+outoffset) = NA_REAL;
      }
      else
      {
        *(efficacy+outoffset) = 1.0 - post_mean/pre_mean;
      }

      *(n_screen+outoffset) = 0.0;
      *(n_pre+outoffset) = static_cast<double>(pre_n);
      *(n_post+outoffset) = static_cast<double>(post_n);

      *(time_count+outoffset) = counter.get_total();

      outn++;
      outoffset += offset;
    }
  }

