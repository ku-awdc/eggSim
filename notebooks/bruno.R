## Example for Bruno

# From CRAN:
library("tidyverse")

# To install eggSim:
remotes::install_github("ku-awdc/eggSim")
library("eggSim")

mean <- 2200
k <- 0.2

## To be adjusted based on estimates from data:
cv_between <- 1.5
cv_within <- 0.75
cv_slide <- 0.25
cv_reduction <- 1
combined_k(cv_between, cv_within, cv_slide, cv_reduction)


mean <- 10
k <- 0.002

## To be adjusted based on estimates from data:
cv_between <- 5
cv_within <- 3
cv_slide <- 1
cv_reduction <- 3
combined_k(cv_between, cv_within, cv_slide, cv_reduction)


# Replicates
rep <- 1
# Individuals
ind <- 100
# Efficacy
eff <- 0.9
# True prevalence:
true_prev <- 0.5


data <- cgpDataSim(rep, ind, 1-eff, mean, cv_between, cv_within, cv_slide, cv_reduction, true_prevalence=true_prev)
# The observed proportion of non-zero counts
data$ObsPrev[1]

data <- data %>%
  select(Individual, Infected, PreMean, PostMean, Day1Mean=ScreenDayMean, Day2Mean=PreDayMean, Day3Mean=PostDay1Mean, Day4Mean=PostDay2Mean) %>%
  mutate(Slide1a = rnbinom(n(), size=1/cv_slide^2, mu=Day1Mean)) %>%
  mutate(Slide1b = rnbinom(n(), size=1/cv_slide^2, mu=Day1Mean)) %>%
  mutate(Slide2a = rnbinom(n(), size=1/cv_slide^2, mu=Day2Mean))

sd(data$Slide1a)/mean(data$Slide1a)
1/k^0.5
