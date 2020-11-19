#' Simulate ERR data from a compound gamma Poisson distribution
#'
#' @param R The number of replicate datasets
#' @param N The maximum number of individuals in total
#' @param reduction The arithmetic mean reduction (may be vectorised)
#' @param community_mean A vector of arithmetic mean pre-treatment EPG in each community
#' @param cv_between A vector of CV reflecting variation in EPG between individuals in each community
#' @param cv_within Day-to-day variation in EPG within an individual
#' @param cv_slide Variation between faecal samples from the same individual and day
#' @param cv_reduction Variation in efficacy between individuals
#' @param overall_mean The overall mean (i.e. the mean of the distribution reflecting community means)
#' @param grams The grams of faeces examined (0.0417g i.e. 1/24 EPG is standard for Kato-Katz)
#'
#' @return A data frame containing the simulated data
#'
#' @examples
#' cgpDataSim(10, 10, c(0.05, 0.1, 0.15), c(800,1000), c(1.5, 1), 1, 1, 1)
#'
#' @importFrom tidyr expand_grid gather spread
#' @importFrom magrittr %>%
#'
#' @export
cgpDataSim <- function(R, N, reduction, community_mean, cv_between, cv_within, cv_slide, cv_reduction, overall_mean = mean(community_mean), grams=1/24){

  stopifnot(length(N) == 1 && N > 0)

  # community_mean and cv_between can be vectorised for different communities:
  C <- length(community_mean)
  stopifnot(C > 0)
  stopifnot(length(community_mean) == C && all(community_mean > 0.0))
  if(length(cv_between)==1) cv_between <- rep(cv_between, C)
  stopifnot(length(cv_between) == C && all(cv_between >= 0.0))
  N <- ceiling(N/C)

  # reduction can be vectorised:
  stopifnot(length(reduction) > 0 && all(reduction >= 0.0))

  # other parameters must be scalar:
  stopifnot(length(cv_within) == 1 && cv_within >= 0.0)
  stopifnot(length(cv_slide) == 1 && cv_slide >= 0.0)
  stopifnot(length(cv_reduction) == 1 && cv_reduction >= 0.0)
  stopifnot(length(overall_mean) == 1 && overall_mean >= 0.0)

  # Helper function:
  rgamcv <- function(n, mu, cv){
    rv <- rgamma(n, 1/cv^2, scale=mu*cv^2)
	## Allow cv of zero to remove this distribution:
    rv[cv<=0] <- mu
    return(rv)
  }

  simdata <- expand_grid(Replicate = 1:R, Community = 1:C, Individual = 1:N, Reduction = reduction) %>%
    mutate(TrueGeometric = NA, TrueArithmetic = 100*(1-Reduction), OverallMean = overall_mean, CommunityMean = community_mean[Community]) %>%
    mutate(CVbetween = cv_between[Community], CVwithin = cv_within, CVslide = cv_slide, CVreduction = cv_reduction) %>%
    mutate(PreMean = rgamcv(n(), CommunityMean, CVbetween)) %>%
    mutate(PostMean = rgamcv(n(), PreMean*Reduction, CVreduction)) %>%
    mutate(ScreenDayMean = rgamcv(n(), PreMean, CVwithin), PreDayMean = rgamcv(n(), PreMean, CVwithin)) %>%
    mutate(PostDay1Mean = rgamcv(n(), PostMean, CVwithin), PostDay2Mean = rgamcv(n(), PostMean, CVwithin)) %>%
    mutate(ScreenSlide = rgamcv(n(), ScreenDayMean, CVslide), PreSlide = rgamcv(n(), PreDayMean, CVslide)) %>%
    mutate(PostSlide1a = rgamcv(n(), PostDay1Mean, CVslide), PostSlide1b = rgamcv(n(), PostDay1Mean, CVslide), PostSlide2 = rgamcv(n(), PostDay2Mean, CVslide)) %>%
    mutate(ScreenCount = rpois(n(), ScreenSlide*grams), PreCount = rpois(n(), PreSlide*grams)) %>%
    mutate(PostCount1a = rpois(n(), PostSlide1a*grams), PostCount1b = rpois(n(), PostSlide1b*grams), PostCount2 = rpois(n(), PostSlide2*grams))

  return(simdata)

}


# # First calculate the arithmetic mean efficacy by brute force:
# effdata <- expand_grid(Individual=1:10^6/C, Community=1:C, Reduction=reduction) %>%
#   mutate(CommunityMean = log(community_mean)[Community], CVbetween = cv_between[Community], CVwithin = cv_within, CVslide = cv_slide, CVreduction = cv_reduction) %>%
#   mutate(PreMean = rnorm(n(), CommunityMean, cv_lsd(CVbetween))) %>%
#   mutate(IndRed = log(Reduction) + rnorm(n(), 0.0, cv_lsd(CVreduction))) %>%
#   mutate(PreDayMean = rnorm(n(), PreMean, cv_lsd(CVwithin))) %>%
#   mutate(PostDayMean = rnorm(n(), PreMean+IndRed, cv_lsd(CVwithin))) %>%
#   mutate(PreDaySlide = exp(rnorm(n(), PreDayMean, cv_lsd(CVslide)))) %>%
#   mutate(PostDaySlide = exp(rnorm(n(), PostDayMean, cv_lsd(CVslide)))) %>%
#   group_by(Reduction) %>%
#   summarise(GeometricEfficacy = 100*(1-exp(mean(log(PostDaySlide)) - mean(log(PreDaySlide)))),
#             ArithmeticEfficacy = 100*(1-mean(PostDaySlide)/mean(PreDaySlide)))

