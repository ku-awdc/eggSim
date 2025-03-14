############################################
##
## Script to re-generate analysis/results
## Matt Denwood, 2025-03-05
## This file is distributed as part of eggSim
## License:  GPL-3
##
############################################

## The tidyverse, qs and remotes packages are available from CRAN
library("tidyverse")
theme_set(theme_light())
library("qs")

## The eggSim package currently must be installed from github
if(!requireNamespace("eggSim")){
  remotes::install_github("ku-awdc/eggSim")
}
library("eggSim")


############################################
## Parameter values
############################################

## General simulation parameters:
iterations <- 1e4
cl <- 10
individ_min <- 100
individ_increment <- 5
performance_max <- 0.999
individ_fig1 <- 380

expand_grid(
  parasite = c("ascaris","hookworm","trichuris"),
  endemicity = c(2,5,15,35,65)
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
  dropout_screen = c(0,10),
  dropout_pre = c(0,20),
) |>
  expand_grid(
    force_inclusion_prob = c(0, 10, 20)
  ) |>
  filter(dropout=="baseline" | force_inclusion_prob==0) ->
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
  pivot_wider(names_from=name, values_from=value) |>
  mutate(efficacy_lower_target = efficacy_lower_target / 100) |>
  mutate(efficacy_expected = efficacy_expected / 100) ->
  parameters_thresholds

## Parameters for analysis type
parameters_analysis <- tibble(analysis_type = c("mean","delta"))

## Parameters for simulated drug efficacy
parameters_efficacy <- tibble(true_efficacy = seq(50,100,by=0.25)/100)

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
expand_grid(
  variant = c("NS_11","NS_12","SSR_11A","SSR_12A","SSR_11B","SSR_12B"),
  min_positive = c(1, 10, 25, 50, 100)
) |>
  mutate(
    design = str_sub(variant, 1, if_else(str_detect(variant, "NS"), 5, 6)),
    method = "kk",
    min_positive_screen = if_else(str_detect(design, "NS"), 0, min_positive),
    min_positive_pre = case_when(
      str_detect(design, "NS") ~ min_positive,
      str_detect(variant, "A") ~ min_positive,
      str_detect(variant, "B") ~ 1,
    ),
    n_day_screen = if_else(str_detect(design, "NS"), 0, 1),
    n_aliquot_screen = if_else(str_detect(design, "NS"), 0, 1),
    n_day_pre = 1,
    n_aliquot_pre = 1,
    n_day_post = 1,
    n_aliquot_post = if_else(str_detect(design, "11"), 1, 2),
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
  left_join(x, new_df, by = join_by(parasite, endemicity, weight)) |>
    mutate(cost_aliquot_post = if_else(str_detect(design, "11"), cost_aliquot_post_11, cost_aliquot_post_12))
}


############################################
## Utility functions
############################################

fix_n_analysis <- function(parameters, iters=iterations, cl=NULL){

  parameters |>
    group_by(framework, analysis_type, parasite, endemicity, mean_epg, min_positive, variant) |>
    group_split() ->
    pars

  seq_along(pars) |>
    lapply(function(i){
      pars[[i]] |>
        distinct(parasite, mean_epg, true_efficacy, endemicity, framework, analysis_type, individ_min, individ_max, min_positive, variant) |>
        mutate(scenario = row_number()) ->
        scenario

      stopifnot(nrow(scenario |> distinct(individ_min, individ_max))==1L)

      pars[[i]] |>
        select(parasite, !any_of(names(scenario))) |>
        unique() |>
        mutate(parameter_set = row_number()) |>
        rowwise() |>
        group_split() ->
        pp

      stopifnot(nrow(scenario)==1L || length(pp)==1L)

      survey_sim(
        n_individ = seq(scenario$individ_min, scenario$individ_max, by=individ_increment),
        scenario = scenario,
        parameters = pp,
        iterations = iters,
        cl = cl,
        output = "summarised",
        analysis = scenario$analysis_type[1]
      ) |>
        left_join(
          scenario |> select(scenario, framework, analysis_type, endemicity, min_positive, variant),
          by="scenario"
        ) |>
        mutate(Positive = (n_Susceptible+n_LowResistant+n_ClassifyFail), Negative = iters-Positive, Performance = Positive/iters)
    }) |>
    bind_rows() |>
    ungroup()
}


vary_n_analysis <- function(parameters, iters=iterations, cl=NULL){

  parameters |>
    rowwise() |>
    group_split() ->
    pars

  seq_along(pars) |>
    pbapply::pblapply(function(i){
    #lapply(function(i){

      if(is.null(cl)) cat("Parameter cluster ", i, " of ", length(pars), "...\n", sep="")

      capture.output({
      pars[[i]] |>
        mutate(individ_min = individ_min, individ_max = individ_min*200) |>
        fix_n_analysis(iters = 100, cl=NULL) ->
        pilot
      })

      suppressWarnings(mod <- mgcv::gam(cbind(Positive, Negative) ~ s(n_individ), family="binomial", data=pilot))
      pilot$predict <- plogis(predict(mod))

      pilot |>
        filter(predict > performance_max) ->
        perf_ok
      if(nrow(perf_ok)==0L) browser()
      perf_ok |>
        arrange(n_individ) |>
        slice(1) |>
        pull(n_individ) ->
        individ_max

      cat(individ_max, "->")
      individ_max <- ceiling(individ_max*0.11)*10
      cat(individ_max, "\n")

      capture.output({
        pars[[i]] |>
        mutate(individ_min = individ_min, individ_max = individ_max) |>
        fix_n_analysis(iters=iters, cl=NULL) ->
          res
      })

      res
    #}) |>
    }, cl=cl) |>
  bind_rows() |>
    ungroup()
}


plot_data <- function(res){
  res |>
    mutate(aborted = (n_FailZeroPre+n_FailPositiveScreen+n_FailPositivePre)) |>
    mutate(Completion = 1 - aborted / (Positive+Negative)) |>
    mutate(MeanCost = cost_mean, StdvCost = sqrt(cost_variance)) |>
    mutate(Power = Positive / (Positive+Negative-aborted)) |>
    mutate(SampleSize = n_individ)
}

plot_data_cost <- function(res){
  res |>
    plot_data() |>
    pivot_longer(c("Performance","Completion","StdvCost","Power","SampleSize")) |>
    mutate(name = factor(name, levels=c("Completion","Power","Performance","SampleSize","StdvCost")))
}

plot_data_ss <- function(res){
  res |>
    plot_data() |>
    pivot_longer(c("Performance","Completion","StdvCost","Power","MeanCost")) |>
    mutate(name = factor(name, levels=c("Completion","Power","Performance","MeanCost","StdvCost")))
}


############################################
## Recreate figure 1
############################################

set.seed(2025-03-05)

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm", endemicity==15),
  parameters_fixed |> filter(design == "NS_11"),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "baseline", force_inclusion_prob == 0),
  parameters_analysis,
  parameters_efficacy
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB"),
    by = "parasite", relationship="many-to-many"
  ) |>
  mutate(individ_min = individ_fig1, individ_max = individ_fig1) |>
  vary_n_analysis(iterations=iterations, cl=cl) ->
  fig_1_data
#qsave(fig_1_data, "notebooks/paper_2025/fig_1_data.rqs")

fig_1_data |>
  mutate(
    Failed = n_failure + n_FailZeroPre,
    Adequate = n_above_cutoffs + n_Susceptible + n_ClassifyFail,
    Reduced = n_below_cutoffs + n_Resistant + n_LowResistant,
    Inconclusive = n_between_cutoffs + n_Inconclusive
  ) |>
  mutate(Total = Failed + Adequate + Reduced + Inconclusive) |>
  select(true_efficacy, efficacy_expected, analysis, Failed, Adequate, Reduced, Inconclusive) |>
  pivot_longer(Failed:Inconclusive, names_to="classification", values_to="tally") ->
  plotdata

## To insert into paper:
summary((plotdata |> filter(classification=="Failed") |> pull(tally)) / iterations * 100)

plotdata |>
  filter(classification != "Failed") |>
  mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
  group_by(efficacy_expected, analysis, true_efficacy) |>
  arrange(classification) |>
  mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
  ungroup() |>
  ggplot(aes(x=true_efficacy, ymin=ymin, ymax=ymax, fill=classification)) +
  geom_ribbon() +
  facet_grid(efficacy_expected ~ analysis)


############################################
## Re-create figure 2
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm"),
  parameters_fixed,
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "baseline", force_inclusion_prob == 0),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

parameters |>
  vary_n_analysis(cl=10) ->
  res
#qsave(res, "notebooks/paper_2025/temp_res.rqs")

res |>
  plot_data_cost() |>
  filter(name=="Performance", value>0.5, value<0.95) |>
  ggplot(aes(x=MeanCost, y=value, col=variant)) +
  geom_line() +
  facet_grid(min_positive ~ endemicity, scales="free") +
  geom_hline(yintercept = 0.8, lty="dashed")
ggsave("fig_performance.pdf", width=10, height=10)

pdf("fig_cost.pdf", width=10, height=10)
res |>
  group_by(min_positive) |>
  group_split() |>
  lapply(function(x){
    plot_data_cost(x) |>
      ggplot(aes(x=MeanCost, y=value, col=variant)) +
      geom_line() +
      facet_grid(name ~ endemicity, scales="free") +
      ylab(NULL) +
      labs(title=str_c("MinPos: ", x$min_positive[1]))
  }) |>
  print()
dev.off()

pdf("fig_sample.pdf", width=10, height=10)
res |>
  group_by(min_positive) |>
  group_split() |>
  lapply(function(x){
    plot_data_ss(x) |>
      ggplot(aes(x=SampleSize, y=value, col=variant)) +
      geom_line() +
      facet_grid(name ~ endemicity, scales="free") +
      ylab(NULL) +
      labs(title=str_c("MinPos: ", x$min_positive[1]))
  }) |>
  print()
dev.off()



ggplot(res, aes(x=cost_mean, y=Performance, col=design)) +
  geom_line() +
  facet_wrap(~endemicity, scales="free_x") +
  geom_hline(yintercept=0.8, lty="dashed")

ggplot(res, aes(x=cost_mean, y=sqrt(cost_variance), col=design)) +
  geom_line() +
  facet_wrap(~endemicity, scales="free") +
  geom_hline(yintercept=0.8, lty="dashed")

ggplot(res, aes(x=cost_mean, y=n_individ, col=design)) +
  geom_line() +
  facet_wrap(~endemicity, scales="free") +
  geom_hline(yintercept=0.8, lty="dashed")

res |>
  mutate(Completion = 1 - (n_FailZeroPre+n_FailPositiveScreen+n_FailPositivePre) / (Positive+Negative)) |>
  ggplot(aes(x=cost_mean, y=Completion, col=design)) +
  geom_line() +
  facet_wrap(~endemicity, scales="free") +
  geom_hline(yintercept=0.8, lty="dashed")

res |>
  mutate(Power = Positive / (Positive+Negative- (n_FailZeroPre+n_FailPositiveScreen+n_FailPositivePre))) |>
  ggplot(aes(x=cost_mean, y=Power, col=design)) +
  geom_line() +
  facet_wrap(~endemicity, scales="free") +
  geom_hline(yintercept=0.8, lty="dashed")

ggplot(res, aes(x=n_individ, y=sqrt(cost_variance), col=design)) +
  geom_line() +
  facet_wrap(~endemicity, scales="free") +
  geom_hline(yintercept=0.8, lty="dashed")

#ggsave("fig_3b.pdf", height=8, width=10)

res |>
  mutate(Performance = Positive / (Positive+Negative-n_FailPositivePre)) |>
  ggplot(aes(x=cost_mean, y=Performance, col=design)) +
  geom_line()

ggplot(res, aes(x=n_individ, y=Performance, col=design)) +
  geom_line()

ggplot(res, aes(x=cost_mean, y=Performance, col=design)) +
  geom_line()

ggplot(res, aes(x=n_individ, y=(n_FailZeroPre+n_FailPositiveScreen+n_FailPositivePre), col=design)) +
  geom_line()

ggplot(res, aes(x=cost_mean, y=(n_FailZeroPre+n_FailPositiveScreen+n_FailPositivePre), col=design)) +
  geom_line()

expand_grid(
  parameters_scenario,
  parameters_fixed,
  parameters_cost,
  parameters_dropadd,
) |>
  mutate(parameter_set = row_number()) |>
  mutate(cost_aliquot_post = if_else(str_detect(design, "11"), cost_aliquot_post_11, cost_aliquot_post_12)) |>
  add_mean_and_cv() |>
  full_join(parameters_thresholds, by = "parasite", relationship="many-to-many") ->
parameters




  |>
  group_by(framework, parasite, endemicity, mean_epg) |>
  group_split() ->
  all_parameters

pp <- all_parameters[[1]]

pp |>
  distinct(parasite, mean_epg) |>
  mutate(scenario = row_number(), true_efficacy = 80) ->
  scenario

pp |>
  select(parasite, !any_of(names(scenario))) ->
  pp

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


stop("INCREASE SAMPLE SIZE WHERE NEEDED!!!")
warning("Improve mechanism of specifying parameter values")




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
