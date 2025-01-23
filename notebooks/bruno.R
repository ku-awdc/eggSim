## Example for Bruno

# From CRAN:
library("tidyverse")

# To install eggSim:
# remotes::install_github("ku-awdc/eggSim")
library("eggSim")

full <- survey_sim(output="full", analysis="delta")
agg <- survey_sim(analysis="delta")


## Adjust for your needs:
mean <- 2200

## Estimates based on data (reduction=0 is an assumption - we did not estimate that):
hw_cv <- c(between=3.0, within=1.25, slide=1.0, reduction=0)
tri_cv <- c(between=1.1, within=0.75, slide=0.7, reduction=0)

## To switch between parasites:
parasite_cv <- hw_cv
parasite_cv <- tri_cv

## For interest:
(parasite_k <- 1/parasite_cv^2)
(total_k <- combined_k(parasite_cv["between"], parasite_cv["within"], parasite_cv["slide"], parasite_cv["reduction"])[1])
# etc


# Replicates
rep <- 1
# Individuals
ind <- 100
# Efficacy
eff <- 0.9
# True prevalence:
true_prev <- 1


data <- cgpDataSim(rep, ind, 1-eff, mean, parasite_cv["between"], parasite_cv["within"], parasite_cv["slide"], parasite_cv["reduction"], true_prevalence=true_prev)

# The observed proportion of non-zero counts
data$ObsPrev[1]

data <- data %>%
  select(Individual, Infected, PreMean, PostMean, Day1Mean=ScreenDayMean, Day2Mean=PreDayMean, Day3Mean=PostDay1Mean, Day4Mean=PostDay2Mean) %>%
  mutate(Slide1a = rnbinom(n(), size=1/parasite_cv["slide"]^2, mu=Day1Mean)) %>%
  mutate(Slide1b = rnbinom(n(), size=1/parasite_cv["slide"]^2, mu=Day1Mean)) %>%
  mutate(Slide2a = rnbinom(n(), size=1/parasite_cv["slide"]^2, mu=Day2Mean))

# Observed total CV based on slide 1a:
sd(data$Slide1a)/mean(data$Slide1a)
# True total CV:
1/total_k^0.5
