#' Simulate ERR data from a compound lognormal Poisson distribution
#'
#' @param R The number of replicate datasets
#' @param N The maximum number of individuals in total
#' @param reduction The geometric mean reduction (may be vectorised)
#' @param community_mean A vector of geometric mean pre-treatment EPG in each community
#' @param cv_between A vector of CV reflecting variation in EPG between individuals in each community
#' @param cv_within Day-to-day variation in EPG within an individual
#' @param cv_slide Variation between faecal samples from the same individual and day
#' @param cv_reduction Variation in efficacy between individuals
#' @param overall_mean The overall mean (i.e. the mean of the distribution reflecting community means)
#' @param edt The egg detection threshold (24 EPG is standard for Kato-Katz)
#'
#' @return A data frame containing the simulated data
#'
#' @examples
#' clpDataSim(10, 10, c(0.05, 0.1, 0.15), c(800,1000), c(1.5, 1), 1, 1, 1)
#'
#' @importFrom tidyr expand_grid
#' @importFrom magrittr %>%
#'
#' @export
clpDataSim <- function(R, N, reduction, community_mean, cv_between, cv_within, cv_slide, cv_reduction, overall_mean = exp(mean(log(community_mean))), edt=24){

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
  cv_lsd <- function(cv) sqrt(log(cv^2 +1))

  simdata <- expand_grid(Replicate = 1:R, Community = 1:C, Individual = 1:N, Reduction = reduction) %>%
    mutate(TrueGeometric = 100*(1-Reduction), TrueArithmetic = NA, OverallMean = overall_mean, CommunityMean = log(community_mean)[Community]) %>%
    mutate(CVbetween = cv_between[Community], CVwithin = cv_within, CVslide = cv_slide, CVreduction = cv_reduction) %>%
    mutate(PreMean = rnorm(n(), CommunityMean, cv_lsd(CVbetween))) %>%
    mutate(PostMean = PreMean + log(Reduction) + rnorm(n(), 0.0, cv_lsd(CVreduction))) %>%
    mutate(ScreenDayMean = rnorm(n(), PreMean, cv_lsd(CVwithin)), PreDayMean = rnorm(n(), PreMean, cv_lsd(CVwithin))) %>%
    mutate(PostDay1Mean = rnorm(n(), PostMean, cv_lsd(CVwithin)), PostDay2Mean = rnorm(n(), PostMean, cv_lsd(CVwithin))) %>%
    mutate(ScreenSlide = exp(rnorm(n(), ScreenDayMean, cv_lsd(CVslide))), PreSlide = exp(rnorm(n(), PreDayMean, cv_lsd(CVslide)))) %>%
    mutate(PostSlide1a = exp(rnorm(n(), PostDay1Mean, cv_lsd(CVslide))), PostSlide1b = exp(rnorm(n(), PostDay1Mean, cv_lsd(CVslide))), PostSlide2 = exp(rnorm(n(), PostDay2Mean, cv_lsd(CVslide)))) %>%
    mutate(ScreenCount = rpois(n(), ScreenSlide/edt), PreCount = rpois(n(), PreSlide/edt)) %>%
    mutate(PostCount1a = rpois(n(), PostSlide1a/edt), PostCount1b = rpois(n(), PostSlide1b/edt), PostCount2 = rpois(n(), PostSlide2/edt))

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

