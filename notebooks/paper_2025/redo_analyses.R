############################################
##
## Script to re-generate analysis/results
## Matt Denwood, 2025-06-24
## This file is distributed as part of eggSim
## License:  GPL-3
##
############################################

## The tidyverse, qs and remotes packages are available from CRAN
library("tidyverse")
theme_set(theme_light())
library("qs")

## The eggSim package currently must be installed from github
## (bayescount-link branch, which requires an in-development version of bayescount)
# if(!requireNamespace("eggSim")) remotes::install_github("ku-awdc/eggSim")
library("eggSim")
if(packageVersion("eggSim") < "0.9.7") stop("You need to install the bayescount-link branch of eggSim")

library("pbapply")
pboptions(use_lb=FALSE)

############################################
## Parameter values
############################################

## General simulation parameters:
iterations <- 1e4
cl <- 10
individ_increment <- 1

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
  dropout_screen = c(0,0.1),
  dropout_pre = c(0,0.2),
) |>
  expand_grid(
    force_inclusion_prob = c(0, 0.05, 0.1, 0.15, 0.2)
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
  min_positive = c(1:10, 25, 50, 100)
) |>
  filter(!str_detect(variant, "A")) |>  # We agreed this makes no sense
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

fix_n_analysis <- function(parameters, iters=iterations, min = 10L, max = 500L, increment=individ_increment, cl=NULL){

  if(!"ParSet" %in% names(parameters)){
    parameters <- parameters |> mutate(ParSet = row_number())
  }

  parameters |>
    rowwise() |>
    group_split() ->
    pars

  if(is.null(cl)){
    lafun <- lapply
  }else{
    if(cl==1){
      lafun <- function(x, ff) pblapply(x, ff, cl=NULL)
    }else{
      lafun <- function(x, ff) pblapply(x, ff, cl=cl)
    }
  }

  pars |>
    lafun(function(pp){

      survey_sim(
        n_individ = seq(min, max, by=increment),
        scenario = pp |> mutate(scenario=1) |> select(parasite, mean_epg, true_efficacy, scenario),
        parameters = pp |> mutate(parameter_set = ParSet) |> select(-mean_epg, -true_efficacy),
        iterations = iters,
        cl = NULL,
        output = "summarised",
        analysis = pp$analysis_type,
        quiet=TRUE
      ) ->
        res

      bind_cols(
        res,
        pp[!names(pp) %in% names(res)]
      ) |>
        mutate(Positive = (n_Susceptible+n_LowResistant), Negative = iters-Positive, Performance = Positive/iters) |>
        select(ParSet, Performance, everything())

    }) |>
    bind_rows() |>
    ungroup()
}


vary_n_analysis <- function(parameters, iters=iterations, performance=NULL, performance_max=0.999, min = 10L, max = 500L, increment=individ_increment, cl=NULL){

  if(!"ParSet" %in% names(parameters)){
    parameters <- parameters |> mutate(ParSet = row_number())
  }

  parameters |>
    slice_sample(prop=1) |>
    rowwise() |>
    group_split() ->
    pars

  stopifnot(nrow(parameters)==length(pars))

  if(is.null(cl)){
    lafun <- lapply
  }else{
    if(cl==1){
      lafun <- function(x, ff) pblapply(x, ff, cl=NULL)
    }else{
      lafun <- function(x, ff) pblapply(x, ff, cl=cl)
    }
  }

  pars |>
    lafun(function(pp){

      try({
        if(is.null(cl)) cat("Parameter cluster ", i, " of ", length(pars), "...\n", sep="")

        ok <- FALSE

        emin <- min
        emax <- max

        while(!ok){
          pp |>
            fix_n_analysis(iters = 100, min=emin, max=emax, increment=10, cl=NULL) ->
            pilot

          suppressWarnings(mod <- mgcv::gam(cbind(Positive, Negative) ~ s(n_individ), family="binomial", data=pilot))
          pilot$predict <- plogis(predict(mod))

          ok <- any(pilot$predict >=  performance_max)

          if(!ok){
            emin <- emax / 2
            emax <- emax * 2
          }
        }

        pilot |>
          filter(predict > performance_max) |>
          arrange(n_individ) |>
          slice(1) |>
          pull(n_individ) ->
          individ_max

        emax <- ceiling(individ_max*0.11)*10

        pp |>
          fix_n_analysis(iters=iters, min = min, max = emax, increment=increment, cl=NULL) |>
          mutate(Status = "OK") |>
          arrange(Performance) |>
          select(ParSet, Status, Performance, everything()) ->
            res

        if(!is.null(performance)){
          lapply(performance, function(pf){
            res |>
              filter(Performance >= pf) |>
              slice(1) |>
              mutate(Target = pf) |>
              select(Target, Low=n_individ)
          }) |>
            bind_rows()  ->
            lr

          lapply(performance, function(pf){
            res |>
              arrange(desc(Performance)) |>
              filter(cummin(Performance) >= pf) |>
              slice(n()) |>
              mutate(Target = pf) |>
              select(Target, High=n_individ)
          }) |>
            bind_rows() ->
            hr

          full_join(lr, hr, by="Target") |>
            mutate(Mean = (Low+High)/2) |>
            select(Target, Mean) |>
            rowwise() |>
            group_split() |>
            lapply(function(x){
              bind_cols(
                x, res
              ) |>
                filter(n_individ >= Mean) |>
                slice(1) |>
                select(-Mean)
            }) |>
            bind_rows() ->
            or

          stopifnot(nrow(or)==length(performance))

          or |>
            select(ParSet, Status, Target, Performance, everything()) ->
            res
        }

        return(res)

      }) -> ss

      if(inherits(ss, "try-error")) return(pp |> mutate(Status = "Failed", MSG = as.character(ss)))

      return(out)

    }) ->
    out

  if(any(sapply(out, \(x) !"Status" %in% names(x) || any(x$Status!="OK")))){
    warning("One or more error encountered")
    return(out)
  }

  ss <- try({
    out |>
      bind_rows() |>
      ungroup() |>
      arrange(ParSet, Performance) ->
      out
  })

  return(out)
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

cols <- c(gg_colour_hue(3),"grey50")
names(cols) <- c("Reduced","Adequate","Inconclusive","Failed")
cols <- c("Reduced" = "#F8766D", "Inconclusive" = "#00BFC4", "Adequate" = "#7CAE00")

expand_grid(
  n_individ = c(20, 50, 100, 250, 300, 500),
  min_positive = c(1),
  endemicity = c(15),
) ->
  all

st <- Sys.time()
all |>
  filter(min_positive <= n_individ/2) |>
  mutate(Row = row_number()) |>
  rowwise() |>
  group_split() |>
  lapply(function(x){

    cat(x$Row, "of", nrow(all), "-", as.numeric(Sys.time()-st, "mins"), "\n")
    print(x)

    expand_grid(
      parameters_scenario |> filter(parasite=="hookworm", endemicity==x$endemicity),
      parameters_fixed |> filter(design == "NS_12", min_positive == 1),
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
      mutate(min_positive_pre = x$min_positive) |>
      fix_n_analysis(iters=iterations, min = x$n_individ, max=x$n_individ, cl=cl) ->
      fig_1_data

    fig_1_data |>
      mutate(
        Adequate = if_else(analysis_type=="delta", n_Susceptible, n_above_cutoffs),
        Reduced = if_else(analysis_type=="delta", n_Resistant + n_LowResistant, n_below_cutoffs),
        Inconclusive = if_else(analysis_type=="delta", n_Inconclusive, n_between_cutoffs),
        Failed = n_FailZeroPre + n_FailPositiveScreen + n_FailPositivePre + if_else(analysis_type=="delta", n_ClassifyFail, 0)
      ) |>
      mutate(Total = Failed + Adequate + Reduced + Inconclusive) |>
      select(true_efficacy, efficacy_expected, analysis, Adequate, Reduced, Inconclusive, Failed) |>
      pivot_longer(Adequate:Failed, names_to="classification", values_to="tally") |>
      mutate(analysis = if_else(analysis=="delta", "hypothesis", analysis)) |>
      bind_cols(x) ->
      plotdata

    plotdata |>
      mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Failed","Reduced"))) |>
      group_by(efficacy_expected, analysis, true_efficacy) |>
      arrange(classification) |>
      mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
      ungroup() |>
      ggplot(aes(x=true_efficacy, ymin=ymin, ymax=ymax, fill=classification)) +
      geom_ribbon() +
      facet_grid(efficacy_expected ~ analysis) +
      theme_bw() +
      geom_vline(aes(xintercept=efficacy_expected)) +
      geom_vline(aes(xintercept=efficacy_expected-0.1)) +
      geom_hline(yintercept=c(0.05,0.95)) +
      scale_fill_manual(values=cols) +
      labs(title = str_c("N = ", x$n_individ, ", MP = ", x$min_positive, ", End = ", x$endemicity, "%")) ->
      plot1

    plotdata |>
      filter(classification!="Failed") |>
      mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
      group_by(efficacy_expected, analysis, true_efficacy) |>
      arrange(classification) |>
      mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
      ungroup() |>
      ggplot(aes(x=true_efficacy, ymin=ymin, ymax=ymax, fill=classification)) +
      geom_ribbon() +
      facet_grid(efficacy_expected ~ analysis) +
      theme_bw() +
      geom_vline(aes(xintercept=efficacy_expected)) +
      geom_vline(aes(xintercept=efficacy_expected-0.1)) +
      geom_hline(yintercept=c(0.05,0.95)) +
      scale_fill_manual(values=cols) +
      labs(title = str_c("N = ", x$n_individ, ", MP = ", x$min_positive, ", End = ", x$endemicity, "% (no failed)")) ->
      plot2

    list(data=plotdata, p1=plot1, p2=plot2)

  }) ->
  plots
#qsave(plots, "notebooks/paper_2025/fig1_res.rqs")

## Figure 1:
plots |>
  lapply(\(x) x$data) |>
  bind_rows() |>
  filter(n_individ==300) |>
  filter(classification != "Failed") |>
  mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
  group_by(efficacy_expected, analysis, true_efficacy, n_individ) |>
  arrange(classification) |>
  mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
  mutate(type = str_c(analysis, " - ", efficacy_expected) |> fct()) |>
  mutate(type = factor(type,
                       levels=c("mean - 0.9", "hypothesis - 0.9", "mean - 0.962", "hypothesis - 0.962"),
                       labels=c(
                         expression("A: Point estimate with T"[l] *"="* " 80% T"[u] *"="* " 90%"),
                         expression("B: Hypothesis testing with T"[l] *"="* " 80% T"[u] *"="* " 90%"),
                         expression("C: Point estimate with T"[l] *"="* " 86.2% T"[u] *"="* " 96.2%"),
                         expression("D: Hypothesis testing with T"[l] *"="* " 86.2% T"[u] *"="* " 96.2%")
                        )
  )) |>
  ungroup() |>
  ggplot(aes(x=true_efficacy*100, ymin=ymin*100, ymax=ymax*100, fill=classification)) +
  geom_ribbon() +
  facet_wrap( ~ type, labeller = label_parsed) +
  # geom_hline(yintercept=c(0.05,0.95)) +
  geom_segment(aes(y=5, x=50, xend=efficacy_expected*100-10), lty="dotted") +
  geom_segment(aes(y=95, x=efficacy_expected*100, xend=100), lty="dotted") +
  geom_vline(aes(xintercept=efficacy_expected*100)) +
  geom_vline(aes(xintercept=efficacy_expected*100-10), lty="dashed") +
  theme_minimal() +
  scale_fill_manual(values=cols, guide = guide_legend(reverse = TRUE)) +
  labs(x = "True efficacy (%)",
       y = "Proportion of iterations (%)",
       fill = "Efficacy classification") +
  theme_minimal() +
  # theme(strip.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = "bottom")
  theme(strip.text = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom")
ggsave("notebooks/paper_2025/fig1.pdf", height=8, width=10)

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm", endemicity==15),
  parameters_fixed,
  parameters_cost |> filter(setting == "Ethiopia")
) |>
  add_mean_and_cv() |>
  distinct(endemicity, mean_epg)

## New figure S1:
plots |>
  lapply(\(x) x$data) |>
  bind_rows() |>
  filter(n_individ!=300, n_individ>=100) |>
  filter(classification != "Failed") |>
  mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
  group_by(efficacy_expected, analysis, true_efficacy, n_individ) |>
  arrange(classification) |>
  mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
  mutate(type = str_c(analysis, " - ", efficacy_expected) |> fct()) |>
  mutate(type = factor(type,
                       levels=c("mean - 0.9", "hypothesis - 0.9", "mean - 0.962", "hypothesis - 0.962"),
                       labels=c(
                         expression(atop("A: Point estimate with", "T"[l] *"="* " 80% T"[u] *"="* " 90%")),
                         expression(atop("B: Hypothesis testing with", "T"[l] *"="* " 80% T"[u] *"="* " 90%")),
                         expression(atop("C: Point estimate with", "T"[l] *"="* " 86.2% T"[u] *"="* " 96.2%")),
                         expression(atop("D: Hypothesis testing with", "T"[l] *"="* " 86.2% T"[u] *"="* " 96.2%"))
                       )
  )) |>
  mutate(n_individ = factor(n_individ, levels=c(20,50,100,250,500), labels=c(
    expression("N "*"="*" 20"), expression("N "*"="*" 50"), expression("N "*"="*" 100"), expression("N "*"="*" 250"), expression("N "*"="*" 500")))
  ) |>
  ungroup() |>
  ggplot(aes(x=true_efficacy*100, ymin=ymin*100, ymax=ymax*100, fill=classification)) +
  geom_ribbon() +
  facet_grid(n_individ ~ type, labeller = label_parsed) +
  geom_segment(aes(y=5, x=50, xend=efficacy_expected*100-10), lty="dotted") +
  geom_segment(aes(y=95, x=efficacy_expected*100, xend=100), lty="dotted") +
  geom_vline(aes(xintercept=efficacy_expected*100)) +
  geom_vline(aes(xintercept=efficacy_expected*100-10), lty="dashed") +
  theme_minimal() +
  scale_fill_manual(values=cols, guide = guide_legend(reverse = TRUE)) +
  labs(x = "True efficacy (%)",
       y = "Proportion of iterations (%)",
       fill = "Efficacy classification") +
  theme_minimal() +
  #theme(strip.text = element_text(size = 12), legend.title = element_text(size = 12))
  theme(strip.text = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom")
ggsave("notebooks/paper_2025/figS1.pdf", height=8, width=12)


############################################
## New figure S2
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm", endemicity==15),
  parameters_fixed |> filter(min_positive == 1),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "baseline", force_inclusion_prob == 0),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) |>
  mutate(min_positive_pre = 1) |>
  fix_n_analysis(iters=iterations, min = 10, max=1000, cl=cl) ->
  fig_S2_data

fig_S2_data |>
  mutate(design = fct(design, levels=c("NS_11","NS_12","SSR_11","SSR_12"))) |>
  mutate(cost_sd = sqrt(cost_variance) / 1e3, cost_mean = cost_mean/1e3, Performance=Performance*1e2) |>
  pivot_longer(c(Performance, cost_mean, cost_sd)) |>
  mutate(Panel = factor(
    name,
    levels = c("Performance","cost_mean","cost_sd"),
    labels = c(
      expression("A: Performance (" * "%" * ")"),
      expression("B: Mean cost"[total]~ "(x1000 US$)"),
      expression("C: Standard deviation of cost"[total]~ "(x1000 US$)")
    )
  )) ->
  ps2d

ps2d |>
  ggplot(aes(x=n_individ, y=value, col=design)) +
  geom_hline(data=tibble(Panel=factor(levels(ps2d$Panel)[1], levels=levels(ps2d$Panel)), yi=80), aes(yintercept=yi), lty="dashed") +
  geom_hline(data=tibble(Panel=factor(levels(ps2d$Panel)[1], levels=levels(ps2d$Panel)), yi=90), aes(yintercept=yi), lty="dotted") +
  geom_line() +
  facet_wrap(~Panel, scales="free_y", ncol=1, labeller = label_parsed) +
  coord_cartesian(xlim=c(50,500)) +
  labs(x="Number of Children", y=NULL) +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1))
ggsave("notebooks/paper_2025/figS2.pdf", width=7, height=6)

## TODO: other ggplot additions package for different y axis breaks


############################################
## Re-create figure 2
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm", endemicity==2),
  parameters_fixed |> filter(min_positive==1),
  parameters_cost,
  parameters_dropadd,
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  filter(
    (dropout == "baseline" & force_inclusion_prob == 0) | # A/B
      (setting == "Ethiopia" & force_inclusion_prob == 0 & dropout != "baseline") | # C
      (setting == "Ethiopia" & dropout == "baseline" & design=="SSR_12") # D
  ) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

parameters |>
  vary_n_analysis(cl=10, iters=iterations) ->
  res
# qsave(res, "notebooks/paper_2025/fig2_res.rqs")

LETTERS[1:4] |>
  lapply(\(x){
    if(x=="A"){
      res |>
        filter(setting == "Ethiopia", dropout == "baseline", force_inclusion_prob == 0) |>
        mutate(Panel = "A: Ethiopian cost")
    }else if(x=="B"){
      res |>
        filter(setting == "Tanzania", dropout == "baseline", force_inclusion_prob == 0) |>
        mutate(Panel = "B: Tanzanian cost")
    }else if(x=="C"){
      res |>
        filter(setting == "Ethiopia", force_inclusion_prob == 0, dropout != "baseline") |>
        mutate(Panel = "C: Drop-outs")
    }else if(x=="D"){
      res |>
        filter(setting == "Ethiopia", dropout == "baseline", design=="SSR_12", force_inclusion_prob%in%c(0,0.1,0.2)) |>
        mutate(Panel = "D: Assessing multiple STH")
    }else{
      stop("ERROR")
    }
  }) |>
  bind_rows() |>
  filter(endemicity==2) |>
  plot_data_cost() |>
  filter(name=="Performance") |> #, value>0.5, value<0.95) |>
  {function(x){
    x |>
      distinct(design, force_inclusion_prob, Panel) |>
      mutate(MeanCost = 70*1e3, value=1) |>
      bind_rows(x)
  }}() |>
  mutate(design = fct(design, levels=c("NS_11","NS_12","SSR_11","SSR_12"))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*1e2, col=design, lty=str_c(force_inclusion_prob*100,"%")), parse=TRUE) +
  geom_line() +
  facet_wrap(~Panel, scales="fixed") +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100), xlim = c(0,30)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Performance (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1))
ggsave("notebooks/paper_2025/fig2.pdf", width=7, height=6)


############################################
## Re-create figures 3 and S3
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm", endemicity!=2),
  parameters_fixed |> filter(min_positive%in%c(1)),
  parameters_cost,
  parameters_dropadd |> filter(force_inclusion_prob == 0),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  filter(dropout=="baseline" | setting=="Ethiopia") |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

parameters |>
  vary_n_analysis(cl=10, iters=iterations, increment=1) ->
  res
# qsave(res, "notebooks/paper_2025/fig3_res.rqs")

res |>
  plot_data_cost() |>
  filter(name=="Performance", dropout=="baseline", setting=="Ethiopia") |> #, value>0.5, value<0.95) |>
  mutate(Panel = str_c(
    factor(endemicity, levels=c(5,15,35,65), labels=LETTERS[1:4]) |> as.character(),
    ": ", round(mean_epg,1), " epg; ", endemicity, "% prev."
  )
  ) |>
  {function(x){
    x |>
      distinct(design, Panel) |>
      mutate(MeanCost = 70*1e3, value=1) |>
      bind_rows(x)
  }}() |>
  filter(!is.na(Panel)) |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, col=design)) +
  geom_line() +
  facet_wrap(~Panel, scales="fixed") +
  ylab("Performance (%)") +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100), xlim = c(0,10)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Performance (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1))
ggsave("notebooks/paper_2025/fig3.pdf", width=7, height=6)


res |>
  plot_data_cost() |>
  filter(name=="Performance") |> #, value>0.5, value<0.95) |>
  {function(x){
    x |>
      distinct(design, setting, dropout, endemicity) |>
      mutate(MeanCost = 70*1e3, value=1) |>
      bind_rows(x)
  }}() |>
  mutate(Row = fct(case_when(
    dropout=="baseline" ~ str_c(setting),
    TRUE ~ str_c(setting, " w/ drop-outs")
  ), levels=c("Ethiopia", "Tanzania", "Ethiopia w/ drop-outs"))) |>
  mutate(Col = fct(str_c(endemicity,"% prev."))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, col=design)) +
  geom_line() +
  facet_grid(Row~Col, scales="fixed") +
  ylab("Performance (%)") +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100), xlim = c(0,10)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Performance (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1)) +
  scale_x_continuous(breaks=seq(0,10,by=2))
ggsave("notebooks/paper_2025/figS3.pdf", width=9, height=7)




############################################
## Re-create figure S4
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="hookworm", endemicity!=2),
  parameters_fixed |> filter(min_positive%in%c(1), design=="SSR_12"),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "baseline"),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

parameters |>
  vary_n_analysis(cl=10, iters=iterations, increment=1) ->
  res


res |>
  plot_data_cost() |>
  filter(name=="Performance") |> #, value>0.5, value<0.95) |>
  {function(x){
    x |>
      distinct(design, force_inclusion_prob, endemicity) |>
      mutate(MeanCost = 70*1e3, value=1) |>
      bind_rows(x)
  }}() |>
  mutate(Col = fct(str_c(endemicity,"% prev."))) |>
  mutate(Padd = fct(str_c(force_inclusion_prob*100,"%"), levels=str_c(c(0,0.05,0.1,0.15,0.2)*100,"%"))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, col=Padd)) +
  geom_line() +
  facet_wrap(~Col, scales="fixed") +
  ylab("Performance (%)") +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100), xlim = c(0,10)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Performance (%)") +
  guides(color = guide_legend(title=bquote(P[add]))) +
  scale_x_continuous(breaks=seq(0,10,by=2))
ggsave("notebooks/paper_2025/figS4.pdf", width=7, height=6)


############################################
## Re-create figure 4
############################################

expand_grid(
  parameters_scenario |> filter(endemicity==15),
  parameters_fixed |> filter(min_positive%in%c(1)),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "baseline", force_inclusion_prob==0),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(framework=="FHT", parasite=="hookworm" | drug=="ALB") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

parameters |>
  vary_n_analysis(cl=10, iters=iterations, increment=1) ->
  res

res |>
  plot_data_cost() |>
  filter(name=="Performance") |> #, value>0.5, value<0.95) |>
  {function(x){
    return(x)
    x |>
      distinct(parasite, drug, design) |>
      mutate(MeanCost = case_when(
        parasite=="ascaris" ~ 1250,
        parasite=="trichuris" ~ 100000,
        drug=="ALB" ~ 7000,
        drug=="MEB" ~ 60000
      ), value=1) |>
      bind_rows(x)
  }}() |>
  mutate(Col = str_c(drug, " vs. ", parasite)) |>
  mutate(Col = factor(str_c(drug," against ", parasite), levels=c(
    "ALB against hookworm", "MEB against hookworm", "ALB against ascaris", "ALB against trichuris"
  ), labels=c(
    "ALB against Hookworm", "MEB against Hookworm", "ALB against Ascaris", "ALB against Trichuris"
  ))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, col=design)) +
  geom_line() +
  facet_wrap(~Col, scales="free_x") +
  ylab("Performance (%)") +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100)) +
  #coord_cartesian(ylim = c(50,100), xlim = c(0, 60)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Performance (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(color = guide_legend(title="Survey design"))
ggsave("notebooks/paper_2025/fig4.pdf", width=7, height=6)



############################################
## Re-create Table S2
############################################

## TODO: Table S2:  add all padd values, and ethiopia & tanzanian costs, all 4 - should be 1300 rows
## Add colum for difference in cost and sample size compared to optimal design for that combo
## Also add variance in costs

expand_grid(
  parameters_scenario |> filter(endemicity%in%c(5,15,35,65)),
  parameters_fixed |> filter(min_positive%in%c(1)),
  parameters_cost, # |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(force_inclusion_prob==0) |> select(starts_with("dropout")),
  parameters_dropadd |> filter(dropout=="baseline") |> select(!starts_with("dropout")),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters


## For Table 3:
parameters |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  vary_n_analysis(cl=10L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resA

parameters |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  vary_n_analysis(cl=10L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resB

parameters |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  vary_n_analysis(cl=10L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resC

bind_rows(
  resA |> mutate(Replicate = "A"),
  resB |> mutate(Replicate = "B"),
  resC |> mutate(Replicate = "C"),
) |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  select(drug, parasite, design, Target, Replicate, n_individ) |>
  mutate(parasite = fct(parasite, levels=c("hookworm","ascaris","trichuris"))) |>
  arrange(drug, parasite, Target, design, Replicate) ->
  res

res |>
  writexl::write_xlsx("notebooks/paper_2025/table_3.xlsx")

# qsave(res, "notebooks/paper_2025/table3_res.rqs")


## For Table S2:

## Takes around 7.5 hours:
parameters |>
  vary_n_analysis(cl=10L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  res

res |>
  group_by(endemicity, dropout, force_inclusion_prob, setting, drug, parasite, Target) |>
  mutate(n_individ_min = min(n_individ), cost_mean_min = min(cost_mean)) |>
  ungroup() |>
  mutate(n_individ_delta = n_individ-n_individ_min, cost_mean_delta = cost_mean-cost_mean_min) |>
  select(drug, parasite, setting, endemicity, dropout, force_inclusion_prob, Target, design, n_individ, n_individ_min, n_individ_delta, cost_mean, cost_mean_min, cost_mean_delta, cost_variance) |>
  arrange(drug, parasite, setting, endemicity, dropout, force_inclusion_prob, Target, design) ->
  res

qsave(res, "~/Desktop/tables2_res.rqs")

stopifnot(nrow(res)==(nrow(parameters)*2L))

qsave(res, "notebooks/paper_2025/tables2_res.rqs")
# res <- qread("notebooks/paper_2025/tables2_res.rqs")


res |>
  filter(cost_mean==cost_mean_min) |>
  count(drug, parasite, endemicity, design) |>
  print(n=Inf)


## Additional plots (not for the paper):

ggplot(res, aes(x=design, y=cost_mean_delta+1)) +
  geom_violin() +
  facet_grid(endemicity ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(trans="log10")

ggplot(res, aes(x=design, y=cost_mean_delta+1)) +
  geom_boxplot() +
  facet_grid(endemicity ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(trans="log10")

ggplot(res, aes(col=design, y=cost_mean_delta+1)) +
  stat_ecdf() +
  facet_grid(endemicity ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(trans="log10")

ggplot(res, aes(col=design, y=cost_mean+1)) +
  stat_ecdf() +
  facet_grid(endemicity ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(trans="log10")

ggplot(res, aes(col=design, y=n_individ_delta+1)) +
  stat_ecdf() +
  facet_grid(endemicity ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(trans="log10")

ggplot(res, aes(col=design, y=n_individ+1)) +
  stat_ecdf() +
  facet_grid(endemicity ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(trans="log10")
