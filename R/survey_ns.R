#' Title
#'
#' @param N_individ
#' @param N_day_pre
#' @param N_sample_pre
#' @param N_day_post
#' @param N_sample_post
#' @param mu_pre
#' @param individ_k
#' @param day_k
#' @param sample_k
#' @param efficacy_a
#' @param efficacy_b
#'
#' @export
survey_ns <- function(N_individ=10, N_day_pre=1, N_sample_pre=1, N_day_post=1,
                      N_sample_post=1, mu_pre=10, individ_k=1, day_k=1,
                      sample_k=1, efficacy_a=1, efficacy_b=1){


  Rcpp_survey_ns(as.integer(N_individ), as.integer(N_day_pre), as.integer(N_sample_pre), as.integer(N_day_post), as.integer(N_sample_post), as.double(mu_pre), as.double(individ_k), as.double(day_k), as.double(sample_k), as.double(efficacy_a), as.double(efficacy_b))

}
