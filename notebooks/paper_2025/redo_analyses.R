############################################
##
## Script to re-generate analysis/results
## Matt Denwood, 2025-03-05
## This file is distributed as part of eggSim
## License:  GPL-3
##
############################################

## The tidyverse and remotes packages are available from CRAN
library("tidyverse")
theme_set(theme_light())

## The eggSim package currently must be installed from github
if(!requireNamespace("eggSim")){
  remotes::install_github("ku-awdc/eggSim")
}
library("eggSim")


############################################
## Parameter values
############################################

## Survey and simulation parameters:
iterations <- 1e3
sample_size <- seq(100,2000,by=5)
true_efficacy <- seq(50,100,by=0.25)

expand_grid(
  parasite = c("ascaris","hookworm","trichuris"),
  endemicity = c(5,15,35,65)
) ->
  parameters_scenario


## CV-related parameters:
tribble(~parasite, ~intercept, ~slope, ~day_cv, ~reduction_cv,
  "ascaris", 0.0158, 0.0019, 1.40, 0.063,
  "hookworm", 0.0162, 0.0222, 1.07, 0.146,
  "trichuris", 0.0098, 0.0444, 0.84, 0.068,
) ->
  parameters_cv

## Parameters for dropout and assessing other parasites:
tibble(
  dropout = c("baseline", "with dropouts"),
  dropout_bl = c(0,10),
  dropout_fu = c(0,20),
) |>
  expand_grid(
    addition = c(0, 10, 20)
  ) |>
  filter(dropout=="baseline" | addition==0) ->
  parameters_dropadd

## Parameters for drug efficacy
tribble(~parasite, ~drug, ~WHO.efficacy_lower_target, ~WHO.efficacy_expected, ~FHT.efficacy_lower_target, ~FHT.efficacy_expected,
  "ascaris", "ALB", 85.0, 95.0, 89.9, 99.9,
  "ascaris", "MEB", 85.0, 95.0, 88.0, 98.0,
  "trichuris", "ALB", 40.0, 50.0, 54.5, 64.5,
  "trichuris", "MEB", 40.0, 50.0, 52.7, 62.7,
  "hookworm", "ALB", 80.0, 90.0, 86.2, 96.2,
  "hookworm", "MEB", 60.0, 70.0, 70.6, 80.6
) |>
  pivot_longer(cols=c(-parasite, -drug)) |>
  separate_wider_delim(name, delim=".", names=c("framework", "name")) |>
  pivot_wider(names_from=name, values_from=value) ->
  parameters_thresholds

## Cost parameters:
bind_rows(
  tibble(setting = "Ethiopia") |>
    mutate(cost_sample = 0.57, cost_aliquot_screen = 1.37,
      cost_aliquot_pre = 1.37, cost_aliquot_post_11 = 1.37,
      cost_aliquot_post_12 = 1.51, cost_salary = 22.5,
      cost_travel = 90),
  tibble(setting = "Tanzania") |>
    mutate(cost_sample = 0.62, cost_aliquot_screen = 0.84,
      cost_aliquot_pre = 0.84, cost_aliquot_post_11 = 0.84,
      cost_aliquot_post_12 = 0.90, cost_salary = 42.73,
      cost_travel = 242.3),
) ->
  parameters_cost


## Fixed parameters:
tibble(
  design = c("NS_11","NS_12","SSR_11","SSR_12")
) |>
  mutate(
    method = "kk",
    n_day_screen = if_else(str_detect(design, "NS"), 0, 1),
    n_aliquot_screen = if_else(str_detect(design, "NS"), 0, 1),
    n_day_pre = 1,
    n_aliquot_pre = 1,
    n_day_post = 1,
    n_aliquot_post = if_else(str_detect(design, "11"), 1, 2),
    min_positive_screen = if_else(str_detect(design, "NS"), 0, 1),
    min_positive_pre = 50,
    aliquot_cv = 0,
    weight = 1/24,
    recovery = 1,
    count_add = 0, count_mult = 1, count_intercept = 2.38, count_coefficient = 0.066, n_technicians = 3, n_team = 4,
    time_demography = 15,
    time_prep_screen = if_else(str_detect(design, "SSR"), 67, 0),
    time_prep_pre = 67,
    time_prep_post = if_else(str_detect(design, "11"), 67, 135),
    time_record = 9,
    alpha = 0.05,
  ) ->
  parameters_fixed

stop("INCREASE SAMPLE SIZE WHERE NEEDED!!!")
warning("Improve mechanism of specifying parameter values")


############################################
## Estimating mean_epg and individual_cv
############################################

add_mean_and_cv <- function(x, mu_max=1e4){
  x |>
    distinct(.data$parasite, .data$endemicity, .data$weight) |>
    left_join(
      parameters_cv,
      by = "parasite"
    ) |>
    # NOTE: dividing by 24 IS necessary here:
    mutate(slope = slope / 24) |>
    rowwise() |>
    group_split() |>
    map(function(y){
      int <- y |> pull("intercept")
      slope <- y |> pull("slope")
      zeros <- 1 - (y |> pull("endemicity"))/100
      cv_d <- y |> pull("day_cv")
      wt <-  y |> pull("weight")
      mean <- optimise(function(mu){
        k_i <- int + slope*mu
        abs(integrate(function(x) dnbinom(0, 1/cv_d^2, mu=x) * dgamma(x, k_i, rate=k_i/(wt*mu)), 0, Inf)$value-zeros)
      }, c(0,mu_max))$minimum
      if(mean >= (mu_max*0.99) || mean <= 0.01) stop("mu_max needs adjustment!")
      y |>
        mutate(mean_epg = mean) |>
        mutate(individ_cv = 1 / sqrt(.data$intercept + .data$slope*.data$mean_epg)) |>
        select(-"intercept", -"slope")
    }, .progress=TRUE) |>
    list_rbind() |>
    # mutate(total_cv = sqrt(day_cv^2 + individ_cv^2 + day_cv^2*individ_cv^2)) |>
    identity() ->
    new_df
  left_join(x, new_df, by = join_by(parasite, endemicity, weight))
}


############################################
## Generate full parameter sets
############################################

expand_grid(
  parameters_scenario,
  parameters_fixed,
  parameters_cost,
  parameters_dropadd,
) |>
  mutate(cost_aliquot_post = if_else(str_detect(design, "11"), cost_aliquot_post_11, cost_aliquot_post_12)) |>
  add_mean_and_cv() |>
  full_join(parameters_thresholds, by = "parasite", relationship="many-to-many") |>
  group_by(framework, parasite, endemicity, mean_epg) |>
  group_split() ->
  all_parameters

pp <- all_parameters[[1]]

pp |>
  distinct(parasite, mean_epg) |>
  mutate(true_efficacy = 80) ->
  scenario

cl <- NULL
survey_sim(
  n_individ = sample_size,
  scenario = scenario,
  parameters = pp[1:10,] |> rowwise() |> group_split(),
  iterations = iterations,
  cl = cl,
  output = "summarised",
  analysis = "mean"
)

tibble(
  parasite = "hookworm",
  day_cv = 1,
  endemicity = c(5,15,35,65),
  weight = 1/24
) |>
  add_mean_cv()

100 * (1-dnbinom(0, 1/1.9^2, mu=283/24))

integrate(function(x) dnbinom(0, 1, mu=x) * dgamma(x, 1, 1), 0, Inf)$value

ii <- 1e5
rgamma_mu <- function(n, k, mu) rgamma(n, k, rate=k/mu)
dd <- rpois(ii, rgamma_mu(ii, 1/1^2, rgamma_mu(ii, 1/1.71^2, 354))/24)
sum(dd!=0)/ii *100

tibble(
  EPG = 1:1000,
) |>
  mutate(
    Weight1 = find_cvi(EPG, 0.01589, 0.00193, 1),
    Weight24 = find_cvi(EPG, 0.01589, 0.00193, 1/24),
  ) |>
  pivot_longer(starts_with("Weight"), names_to="Weight", values_to="cv_i") |>
  ggplot(aes(x=EPG, y=cv_i, col=Weight)) +
  geom_line() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  ylim(0,NA)
ggsave("calculating_cv_ind.pdf", width=6, height=5)

find_cvi(1:250, 0.01589, 0.00193, 1)
find_cvi(1:250, 0.01589, 0.00193, 1/24)

find_mu <- function(prev, k){

}
find_mu(0.05, 0.023)*c(34,24,10)
find_mu(0.15, 0.053)*c(34,24,10)


alpha <- c(49.5,30.2,53.8)
beta <- c(12.4, 56.1, 17.9)
mean <- alpha / (alpha+beta)
var <- (alpha*beta) / ((alpha+beta)^2 * (alpha+beta+1))
mean
qbeta(0.025,alpha,beta)
qbeta(0.975,alpha,beta)
ss <- rbeta(1e4,alpha[1],beta[1])
sd(ss)/mean(ss)
sqrt(var)/mean

iters <- 1e5
