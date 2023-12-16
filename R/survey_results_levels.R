#' Explanation of the meaning of the results from survey_sim
#'
#' @export
survey_results_levels <- function(){

  reslevs <- Rcpp_results_levels()
  tibble(
    #ColName = str_c("n_results_", seq_along(reslevs)-1L),
    LevelName = factor(reslevs, levels=reslevs),
    SurveySuccess = Rcpp_results_success()
  ) |>
    full_join(
      tribble(
        ~LevelName, ~Explanation,
        "FailZeroPre", "The survey failed due to zero pre-treatment counts. Note that this classification takes precedence over FailPositiveScreen and FailPositivePre.",
        "FailPositiveScreen", "The survey failed because the number of positive individuals at screening was lower than the specified min_positive_screen parameter. This only applies to the SSR design.",
        "FailPositivePre", "The survey failed because the number of positive individuals at pre-treatment was lower than the specified min_positive_pre parameter.",
        "EfficacyBelowLT", "The observed efficacy was below the specified efficacy_lower_target parameter. This only applies to the mean analysis method.",
        "EfficacyAboveLT", "The observed efficacy was equal to or above the specified efficacy_lower_target parameter. This only applies to the mean analysis method (and the remaining classifications do not apply for the mean analysis method).",
        "ClassifyFail", "The survey succeeded but the delta method failed to produce valid confidence intervals for the data, either because of a 100% observed reduction or perfect correlation between pre- and post-treatment data.",
        "Resistant", "The upper confidence interval was less than efficacy_expected and the lower confidence interval was less than the efficacy_lower_target.",
        "LowResistant", "The upper confidence interval was less than efficacy_expected and the lower confidence interval was greater than or equal to the efficacy_lower_target.",
        "Inconclusive", "The upper confidence interval was greater than or equal to the efficacy_expected and the lower confidence interval was less than the efficacy_lower_target.",
        "Susceptible", "The upper confidence interval was greater than or equal to the efficacy_expected and the lower confidence interval was greater than or equal to the efficacy_lower_target."
      ),
      by="LevelName"
    )
}

results_to_factor <- function(x){
  if(any(is.na(x))) warning("One or more missing values in results")
  reslevs <- Rcpp_results_levels()
  x <- factor(x, levels=0:(length(reslevs)-1L), labels=reslevs)
  if(any(is.na(x))) warning("One or more missing values in results after factorising")
  return(x)
}


