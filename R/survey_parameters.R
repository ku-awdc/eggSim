survey_samples <- function(n_label = "example",
               n_day_screen=1, n_aliquot_screen=1,
               n_day_pre=1, n_aliquot_pre=1,
               n_day_post=1, n_aliquot_post=1){

  data.frame(n_label=n_label, n_day_screen=n_day_screen, n_aliquot_screen=n_aliquot_screen, n_day_pre=n_day_pre, n_aliquot_pre=n_aliquot_pre, n_day_post=n_day_post, n_aliquot_post=n_aliquot_post)

}

