#' @import ggplot2
autoplot.eggSim <- function(object, ...){

  stop("The autoplot method is a work in process - sorry!")

  if(!inherits("means", object)){
    object <- getsummary(object, attr(object, "type", TRUE))
  }

  fig2out$Prevalence <- factor(paste0(round(fig2out$ObsPrev*100), "%"))

  bias <- fig2out %>%
    select(.data$Design, .data$Prevalence, .data$Target, Estimate=.data$MedianBias) %>%
    mutate(Type = "Median Bias")
  variance <- fig2out %>%
    select(.data$Design, .data$Prevalence, .data$Target, Estimate=.data$Variance) %>%
    mutate(Estimate = sqrt(.data$Estimate), Type = "Std. Dev.")
  relvariance <- variance %>%
    filter(.data$Design=="NS") %>%
    select(.data$Prevalence, .data$Target, Reference=.data$Estimate) %>%
    full_join(variance, by = c("Prevalence", "Target")) %>%
    mutate(Estimate = (.data$Estimate / .data$Reference), Type = "Relative SD")

  plotdata <- bind_rows(bias, variance, relvariance) %>%
    mutate(Type = factor(.data$Type, levels=c("Median Bias", "Std. Dev.", "Relative SD")))

  yintdat <- plotdata %>%
    distinct(.data$Design, .data$Type) %>%
    mutate(yint = case_when(
      .data$Type == "Relative Variance" ~ 1,
      .data$Type == "Relative SD" ~ 1,
      TRUE ~ 0
    ))


}
