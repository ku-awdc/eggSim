#' Calculate simulation parameters for a given set of (artithmetic) mean counts, coefficients of variation and correlation
#'
#' @param mu_pre The pre-treatment arithmetic mean count
#' @param mu_post The post-treatment arithmetic mean count
#' @param cv_pre The extra-Poisson coefficient of variation (ratio of standard deviation to mean) in the pre-treatment data
#' @param cv_post The extra-Poisson coefficient of variation (ratio of standard deviation to mean) in the post-treatment data
#' @param correlation The extra-Poisson correlation between pre- and post-treatment lambda parameters (GP) or log-lambda parameters (LP)
#' @param geomgamma Should the geometric mean reduction for the gamma-Poisson data be estimated by simulation?
#'
#' @return a list with components reflecting the necessary simulation parameters for correlated gamma-Poisson and correlatedlognormal-Poisson distributions
#'
#' @examples
#' ERRparameters(100, 10, 2, 2.5, 0.4, TRUE)
#'
#' @export
ERRparameters <- function(mu_pre, mu_post, cv_pre, cv_post, correlation, geomgamma=FALSE)
{
  stopifnot(length(mu_pre) == 1 && mu_pre > 0)
  stopifnot(length(mu_post) == 1 && mu_post >= 0)
  stopifnot(length(cv_pre) == 1 && cv_pre > 0)
  stopifnot(length(cv_post) == 1 && cv_post > 0)
  stopifnot(length(correlation) == 1 && correlation >= 0 && correlation <= 1)

  # Lognormal parameters:
  sds <- sqrt(log(c(cv_pre, cv_post)^2 +1))
  lmus <- log(c(mu_pre, mu_post)) - ((sds^2)/2)
  cov <- correlation * sds[1] * sds[2]
  sigma <- matrix(c(sds[1]^2, cov, cov, sds[2]^2), ncol=2)
  lnormals <- list(lmus=lmus, sigma=sigma)

  # Gamma parameters:
  k <- 1/c(cv_pre, cv_post)^2
  kc <- sqrt(k[1]*k[2]) / correlation
  gammas <- list(beta=c(kc, 1.0), g1=c(k[1], kc-k[1]), g2=c(k[2], kc-k[2]), adj=c(mu_pre/k[1], mu_post/k[2]))
  if(geomgamma){
    N <- 10^6
    corrcomp <- rgamma(N, kc, 1.0)
    pre <- rbeta(N, k[1], kc-k[1]) * corrcomp * mu_pre/k[1]
    post <- rbeta(N, k[2], kc-k[2]) * corrcomp * mu_post/k[2]
    gmg <- 100*(1-exp(mean(log(post))-mean(log(pre))))
  }else{
    gmg <- NA
  }

  # Implied parameters:
  efficacies <- matrix(c(rep(100*(1-mu_post/mu_pre), 2), 100*(1-exp(lmus[2]-lmus[1])), gmg), ncol=2, dimnames=list(c('lognormal','gamma'), c('arithmetic','geometric')))

  return(list(lnormals=lnormals, gammas=gammas, efficacies=efficacies))

}

