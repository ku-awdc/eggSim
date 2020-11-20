#' Calculate the effective overall k based on coefficients of variation
#'
#' @param cv_between A vector of CV reflecting variation in EPG between individuals in each community
#' @param cv_within Day-to-day variation in EPG within an individual
#' @param cv_slide Variation between faecal samples from the same individual and day
#' @param cv_reduction Variation in efficacy between individuals
#'
#' @return a vector of pre_k and post_k
#'
#' @examples
#' combined_k(1.5, 0.75, 0.25, 0)
#' combined_k(1, 0.75, 0, 0)
#' (1/1^2 * 1/0.75^2) / (1/1^2 + 1/0.75^2 + 1)
#'
#' @export
combined_k <- function(cv_between, cv_within, cv_slide, cv_reduction){

	cveff <- sqrt(cv_between^2 + cv_within^2 + (cv_between^2*cv_within^2))
	cvpre <- sqrt(cveff^2 + cv_slide^2 + (cveff^2*cv_slide^2))
	cvpost <- sqrt(cvpre^2 + cv_reduction^2 + (cvpre^2*cv_reduction^2))

	return(c(pre_k=1/cvpre^2, post_k=1/cvpost^2))
}
