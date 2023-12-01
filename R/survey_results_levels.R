#' Explanation of the meaning of the results from survey_sim
#'
#' @export
survey_results_levels <- function(){

  reslevs <- Rcpp_results_levels()
  tibble(
    ColName = str_c("n_results_", seq_along(reslevs)-1L),
    LevelName = factor(reslevs, levels=reslevs),
    SurveySuccess = Rcpp_results_success(),
    Explanation = "TODO"
  )

}

results_to_factor <- function(x){
  if(any(is.na(x))) warning("One or more missing values in results")
  reslevs <- Rcpp_results_levels()
  x <- factor(x, levels=0:(length(reslevs)-1L), labels=reslevs)
  if(any(is.na(x))) warning("One or more missing values in results after factorising")
  return(x)
}


