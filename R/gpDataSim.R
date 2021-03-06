#' Simulate ERR data from a correlated gamma distribution
#'
#' @param N The number of (paired) observations to simulate
#' @param mu_pre The pre-treatment arithmetic mean eggs per gram
#' @param mu_post The post-treatment arithmetic mean eggs per gram
#' @param cv_pre The extra-Poisson coefficient of variation (ratio of standard deviation to mean) in the pre-treatment data
#' @param cv_post The extra-Poisson coefficient of variation (ratio of standard deviation to mean) in the post-treatment data
#' @param correlation The extra-Poisson correlation between pre- and post-treatment lambda parameters
#' @param parameters Alternative way to specify the parameters as output of the ERRparameters function
#'
#' @return A data frame containing the simulated data
#'
#' @examples
#' gpDataSim(10, 10, 1, 1, 1, 0.4)
#'
#' @export
gpDataSim <- function(N, mu_pre, mu_post, cv_pre, cv_post, correlation, parameters=NULL){
  stopifnot(length(N) == 1 && N > 0)
  stopifnot(length(mu_pre) == 1 && mu_pre > 0)
  stopifnot(length(mu_post) == 1 && mu_post >= 0)
  stopifnot(length(cv_pre) == 1 && cv_pre > 0)
  stopifnot(length(cv_post) == 1 && cv_post > 0)
  stopifnot(length(correlation) == 1 && correlation >= 0 && correlation <= 1)

  if(!is.null(parameters)){
    stopifnot(!is.null(parameters$gammas))
  }else{
    parameters <- ERRparameters(mu_pre, mu_post, cv_pre, cv_post, correlation, FALSE)
  }


  corrcomp <- rgamma(N, parameters$gammas$beta[1], parameters$gammas$beta[2])
  pre <- rbeta(N, parameters$gammas$g1[1], parameters$gammas$g1[2]) * corrcomp * parameters$gammas$adj[1]
  post <- rbeta(N, parameters$gammas$g2[1], parameters$gammas$g2[2]) * corrcomp * parameters$gammas$adj[2]

  return(data.frame(ID = 1:N, preLambda = pre, postLambda = post))

}
