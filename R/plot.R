#' Title
#'
#' @param object
#' @param ...
#'
#' @return
#'
#' @import ggplot2
#'
#' @examples
#' means <- eggSim(c(0.2,0.1,0.05), cv_reduction=0, R=10^2, parallelise=FALSE)
#' ggplot2::autoplot(means)
#'
#' @export
autoplot.eggSim <- function(object, ...){

  if(!inherets("means", object)){
    object <- getmeans(object, attr(object, "type", TRUE))
  }

  fig2out$Prevalence <- factor(str_c(round(fig2out$ObsPrev*100), "%"))

  bias <- fig2out %>%
    select(Design, Prevalence, Target, Estimate=MedianBias) %>%
    mutate(Type = "Median Bias")
  variance <- fig2out %>%
    select(Design, Prevalence, Target, Estimate=Variance) %>%
    mutate(Estimate = sqrt(Estimate), Type = "Std. Dev.")
  relvariance <- variance %>%
    filter(Design=="NS") %>%
    select(Prevalence, Target, Reference=Estimate) %>%
    full_join(variance, by = c("Prevalence", "Target")) %>%
    mutate(Estimate = (Estimate / Reference), Type = "Relative SD")

  plotdata <- bind_rows(bias, variance, relvariance) %>%
    mutate(Type = factor(Type, levels=c("Median Bias", "Std. Dev.", "Relative SD")))

  yintdat <- plotdata %>%
    distinct(Design, Type) %>%
    mutate(yint = case_when(
      Type == "Relative Variance" ~ 1,
      Type == "Relative SD" ~ 1,
      TRUE ~ 0
    ))


}
