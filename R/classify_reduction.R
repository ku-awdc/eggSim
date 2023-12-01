#' Extract clasification and confidence intervals, using the same underlying code as for survey_sim
#'
#' @param pre integer vector of pre-treatment data
#' @param post integer vector of post-treatment data
#' @param alpha desired confidence level (0.025 = 95% CI)
#' @param efficacy_expected desired expected efficacy target (proportion, not %)
#' @param efficacy_lower_target desired lower efficacy target (proportion, not %)
#'
#' @export
classify_reduction <- function(pre, post, efficacy_expected=0.95, efficacy_lower_target=0.85, alpha=0.025){

  stopifnot(length(efficacy_expected)==1, length(efficacy_lower_target)==1, length(alpha)==1, alpha>0, alpha<1)
  stopifnot(is.numeric(pre), is.numeric(post), !is.na(pre), !is.na(post), length(pre)==length(post), (pre %% 1)==0L, (post %% 1)==0L)

  cs <- eggSim:::CountSummariseTest$new(as.double(efficacy_expected), as.double(efficacy_lower_target), as.double(alpha))
  cs$add_counts(as.integer(pre), as.integer(post))
  print(cs)

  cs$result_delta

}
