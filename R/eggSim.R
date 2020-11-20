#' Simulate and summarise egg reduction rate data
#'
#' @param R The number of replicate datasets
#' @param N The maximum number of individuals in total
#' @param reduction The arithmetic mean reduction (may be vectorised)
#' @param community_mean A vector of arithmetic mean pre-treatment EPG in each community
#' @param cv_between A vector of CV reflecting variation in EPG between individuals in each community
#' @param cv_within Day-to-day variation in EPG within an individual
#' @param cv_slide Variation between faecal samples from the same individual and day
#' @param cv_reduction Variation in efficacy between individuals
#' @param edt The egg detection threshold (24 EPG is standard for Kato-Katz)
#' @param N A named vector of number of included individuals depending on design type
#' @param count Logical flag to base means on count data or the underlying rates
#' @param log_constant A constant to add to the count data before calculating geometric means (ignored if count==FALSE)
#' @param screen_threshold The threshold count on which to screen individuals
#' @param maxN The maximum number of simulated individuals (on which subsampling occurs)
#'
#' @return A data frame containing the simulated data
#'
#' @examples
#'
#'
#' @export
eggSim <- function(reduction, budget=600, second_slide_cost = 0.621, max_screen = 0.9, community_mean=c(800,1000,1200), cv_between=c(0.8,1,1.2), cv_within=1, cv_slide=1, cv_reduction=1, count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0, grams=1/24, R=10^3, type="gamma")
{

  stopifnot(length(type) == 1 && type %in% "gamma")

  # TODO: argument checks
  # TODO: allow community_mean and cv_between to be a matrix so we can alter the population mean values rather than reduction (i.e. make figure 2 as planned)

  # First simulate gamma and lognormal data:
  cat('Simulating gamma data...\n')
  simdatagp <- cgpDataSim(R, maxN, reduction, community_mean, cv_between, cv_within, cv_slide, cv_reduction, edt=edt)
  cat('Simulating lognormal data...\n')
  simdatalp <- clpDataSim(R, maxN, reduction, community_mean, cv_between, cv_within, cv_slide, cv_reduction, edt=edt)

  # Then summarise:
  cat('Summarising gamma data...\n')
  meansgp <- design_means(simdatagp, N=N, count=count) %>% mutate(ComparisonArithmetic = TrueArithmetic, ComparisonGeometric = BestGeometric, IncludeNS = "Arithmetic", Data = if(count) 'Gamma-Poisson' else 'Gamma')
  cat('Summarising lognormal data...\n')
  meanslp <- design_means(simdatalp, N=N, count=count) %>% mutate(ComparisonArithmetic = BestArithmetic, ComparisonGeometric = TrueGeometric, IncludeNS = "Geometric", Data = if(count) 'Lognormal-Poisson' else 'Lognormal')

  # Then join and make plots:
  cat('Creating output...\n')
  means <- bind_rows(meansgp, meanslp) %>%
    select(Design, ComparisonArithmetic, ComparisonGeometric, IncludeNS, Data, OverallMean, Count, ArithmeticEfficacy, GeometricEfficacy) %>%
    gather(Type, Efficacy, -Design, -ComparisonArithmetic, -ComparisonGeometric, -IncludeNS, -Data, -OverallMean, -Count) %>%
    mutate(Type = gsub('Efficacy','',Type), Target = ifelse(Type=='Arithmetic', ComparisonArithmetic, ComparisonGeometric)) %>%
    mutate(Set = paste0(Data, ' - ', Type), Cheating = Design=='NS' & IncludeNS!=Type) %>%
    group_by(Design, Set, Data, Type, Cheating, Target) %>%
    summarise(Bias = mean(Efficacy - Target), LCI = Bias - 1.96*sd(Efficacy - Target)/sqrt(n()), UCI = Bias + 1.96*sd(Efficacy - Target)/sqrt(n()), Variance = var(Efficacy)) %>%
    ungroup()

  # Remove bias estimates where it is cheating:
  means$Bias[means$Cheating] <- NA
  means$LCI[means$Cheating] <- NA
  means$UCI[means$Cheating] <- NA
  means$Cheating <- NULL

  return(means)

}

