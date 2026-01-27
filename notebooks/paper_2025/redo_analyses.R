############################################
##
## Script to re-generate analysis/results
## Matt Denwood, 2025-07-29
## This file is distributed as part of eggSim
## License:  GPL-3
##
############################################


## source("~/Documents/GitHub/eggSim/notebooks/paper_2025/redo_analyses.R")

## Create a results folder:
reswd <- file.path("~/Desktop", paste0("eggsimres_", strftime(Sys.Date(), "%Y-%m-%d")))
reswd <- file.path("~/Documents/Research/Papers/Deo paper", "eggsimres_2025-07-29")
if(dir.exists(reswd)){
  # stop("Path ", reswd, " already exists")
}else{
  dir.create(reswd)
}
cwd <- getwd()
on.exit(setwd(cwd))
setwd(reswd)


## The tidyverse, qs and remotes packages are available from CRAN
library("tidyverse")
theme_set(theme_light())
library("qs")
stopifnot(requireNamespace(c("ggh4x","mgcv","writexl")))
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


## The eggSim package currently must be installed from github
## (bayescount-link branch, which requires a currently in-development version of bayescount ... sorry)
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
  dropout = c("no dropouts", "with dropouts"),
  dropout_screen = c(0,0.1),
  dropout_pre = c(0,0.2),
) |>
  expand_grid(
    force_inclusion_prob = c(0, 0.05, 0.1, 0.15, 0.2)
  ) |>
  filter(dropout=="with dropouts" | force_inclusion_prob==0.1) ->
  parameters_dropadd

## Parameters for analysis type
parameters_analysis <- tibble(analysis_type = c("mean","delta"))

## Parameters for simulated drug efficacy
parameters_efficacy <- tibble(true_efficacy = seq(25,100,by=0.25)/100)

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

## Parameters for drug efficacy with Moderate targets:
tribble(~parasite, ~drug, ~WHO.efficacy_lower_target, ~WHO.efficacy_expected, ~FHT.efficacy_lower_target, ~FHT.efficacy_expected,
        "ascaris", "ALB", 85.0, 95.0, 99.6, 99.9,
        "ascaris", "MEB", 85.0, 95.0, 94.0, 98.0,
        "trichuris", "ALB", 40.0, 50.0, 30.0, 64.5,
        "trichuris", "MEB", 40.0, 50.0, 26.0, 62.7,
        "hookworm", "ALB", 80.0, 90.0, 91.0, 96.2,
        "hookworm", "MEB", 60.0, 70.0, 56.0, 80.6
) |>
  pivot_longer(cols=c(-parasite, -drug)) |>
  separate_wider_delim(name, delim=".", names=c("framework", "name")) |>
  pivot_wider(names_from=name, values_from=value) |>
  mutate(efficacy_lower_target = efficacy_lower_target / 100) |>
  mutate(efficacy_expected = efficacy_expected / 100) ->
  parameters_thresholds

## Add "cross-over" targets:
parameters_thresholds |>
  group_by(parasite, drug) |>
  summarise(efficacy_lower_target = efficacy_expected[framework=="WHO"],
            efficacy_expected = efficacy_expected[framework=="FHT"],
            .groups="drop") ->
  ptx

bind_rows(
  parameters_thresholds |> full_join(parameters_analysis |> mutate(framework=if_else(analysis_type=="mean", "WHO","FHT")), by=join_by(framework)),
  ptx |> mutate(framework="WHO-X") |> mutate(analysis_type="mean"),
  ptx |> mutate(framework="FHT-X") |> mutate(analysis_type="delta"),
) |>
  arrange(parasite, drug, efficacy_expected, efficacy_lower_target, desc(framework)) ->
  parameters_thresholds_fig1


## Other thresholds:
parameters_all_thresholds <- structure(list(parasite = c("ascaris", "ascaris", "ascaris",
"ascaris", "ascaris", "ascaris", "ascaris", "ascaris", "hookworm",
"hookworm", "hookworm", "hookworm", "hookworm", "hookworm", "hookworm",
"hookworm", "trichuris", "trichuris", "trichuris", "trichuris",
"trichuris", "trichuris", "trichuris", "trichuris"), drug = c("ALB",
"ALB", "ALB", "ALB", "MEB", "MEB", "MEB", "MEB", "ALB", "ALB",
"ALB", "ALB", "MEB", "MEB", "MEB", "MEB", "ALB", "ALB", "ALB",
"ALB", "MEB", "MEB", "MEB", "MEB"), Effort = structure(c(1L,
2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L,
2L, 3L, 4L, 1L, 2L, 3L, 4L), levels = c("Easy", "Moderate", "Hard",
"Extreme"), class = "factor"), efficacy_expected = c(0.999, 0.999,
0.999, 0.999, 0.98, 0.98, 0.98, 0.98, 0.962, 0.962, 0.962, 0.962,
0.806, 0.806, 0.806, 0.806, 0.645, 0.645, 0.645, 0.645, 0.627,
0.627, 0.627, 0.627), End_05 = c(0.98000740818078, 0.990830333208262,
NA, NA, 0.907538011164544, 0.936058685854414, NA, NA, 0.862464471953256,
0.90089659602472, NA, NA, 0.42727999992868, 0.563776557080032,
NA, NA, 0.0957487418508169, 0.28209873850991, NA, NA, 0.0440084450946148,
0.246110163849649, NA, NA), End_15 = c(0.993427204213393, 0.995832551456478,
0.99709279710343, NA, 0.906668158482266, 0.935102071504428, 0.951675872849942,
NA, 0.881197573685397, 0.910349001551201, 0.928682198443676,
NA, 0.404699254311741, 0.550244987526767, 0.639713006438012,
NA, 0.081105144023385, 0.280693100509651, 0.406002195814648,
NA, 0.0331521092157447, 0.244545084143181, 0.376149884509504,
NA), End_35 = c(0.993740885108069, 0.996446127680408, 0.997447414327719,
0.998016565762301, 0.905246839752115, 0.940564308207604, 0.954958654624399,
0.963757972690939, 0.872912975723795, 0.913865489566585, 0.930850285141987,
0.941501220168466, 0.38099323191139, 0.574970944797789, 0.656810059076236,
0.706876671171625, 0.0502330987117909, 0.316874757614194, 0.430342491212717,
0.501998357976262, -0.00212824491934049, 0.28227106091769, 0.401273363505762,
0.47814432931775), End_65 = c(NA, 0.996356377666566, 0.997611892167323,
0.998109874148231, NA, 0.937508701725692, 0.956934537742248,
0.965087140677914, NA, 0.90930774318794, 0.932695956938442, 0.942862838269598,
NA, 0.555575673512074, 0.665789315056734, 0.713020040470797,
NA, 0.287535509047836, 0.443030857185861, 0.51027300243537, NA,
0.252368640757529, 0.415161036565098, 0.48562532324896), Average = c(0.989058499167414,
0.994866347502929, 0.997384034532824, 0.998063219955266, 0.906484336466308,
0.937308441823035, 0.954523021738863, 0.964422556684426, 0.872191673787483,
0.908604707582611, 0.930742813508035, 0.942182029219032, 0.404324162050604,
0.561142040729166, 0.654104126856994, 0.709948355821211, 0.0756956615286643,
0.291800526420398, 0.426458514737742, 0.506135680205816, 0.0250107697970063,
0.256323737417012, 0.397528094860121, 0.481884826283355), Using = c(0.989,
0.996, 0.997, 0.998, 0.91, 0.94, 0.95, 0.96, 0.87, 0.91, 0.93,
0.94, 0.4, 0.56, 0.65, 0.71, 0.08, 0.3, 0.43, 0.51, 0.03, 0.26,
0.4, 0.48)), row.names = c(NA, -24L), class = c("tbl_df", "tbl",
"data.frame"))

with(full_join(
  parameters_all_thresholds |>
    filter(Effort=="Moderate") |>
    select(parasite, drug, full=Using),
  parameters_thresholds |>
    filter(framework=="FHT") |>
    select(parasite, drug, partial=efficacy_lower_target),
  by = join_by(parasite, drug)
), stopifnot(full==partial))


## Save parameters:
qsavem(iterations, parameters_scenario, parameters_cv, parameters_dropadd, parameters_analysis, parameters_efficacy, parameters_cost, parameters_fixed, parameters_thresholds, parameters_thresholds_fig1, parameters_all_thresholds, file="parameters.rqm")


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

fix_n_analysis <- function(parameters, iters=iterations, min = 10L, max = 500L, increment=1L, cl=NULL){

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


vary_nim_analysis <- function(parameters, performance=0.8, min = 100L, max = 1000L, NIMlen=25L, cl=NULL){

  stopifnot(length(performance)==1L)

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

        tibble(
          Iterations = c(100,100,1000,10000),
          Scaling = c(0.2, 0.15, 0.1, 0.05),
          GAM = c(FALSE, TRUE, TRUE, TRUE)
        ) ->
          pgrs

        nimrange <- c(10^-4, 5)
        for(r in seq_len(nrow(pgrs))){

          (10^seq(log10(nimrange[2]),log10(nimrange[1]),length=NIMlen)) |>
            as.list() |>
            lapply(\(x){
              pp |>
                mutate(efficacy_lower_target = efficacy_expected - x) |>
                fix_n_analysis(iters = pgrs$Iterations[r], min=min, max=max, increment=10, cl=NULL) |>
                mutate(NIM = x)
            }) |>
            bind_rows() ->
            pilot

          pilot |>
            group_by(n_individ) |>
            group_split() |>
            lapply(\(x){
              obs <- x |> arrange(NIM) |> filter(Performance > performance) |> slice(1)
              if(pgrs$GAM[r]){
                ss <- try(suppressWarnings(mod <- mgcv::gam(cbind(Positive, Negative) ~ s(log(NIM)), family="binomial", data=x)))
                if(inherits(ss, "try-error")) browser()

                afun <- with(x |> mutate(logNIM = log(NIM)), approxfun(logNIM ~ predict(mod)))
                tibble(
                  SampleSize = x$n_individ[1],
                  Performance = performance,
                  MinNIM = afun(qlogis(performance*(1-pgrs$Scaling[r]))) |> exp(),
                  BestNIM = afun(qlogis(performance)) |> exp(),
                  MaxNIM = afun(qlogis(performance*(1+pgrs$Scaling[r]))) |> exp(),
                  ObsPerf = obs$Performance,
                  ObsNIM = obs$NIM
                )
              }else{
                tibble(
                  SampleSize = x$n_individ[1],
                  Performance = performance,
                  MinNIM = x |> arrange(NIM) |> filter(Performance > 0.5) |> slice(1) |> pull(NIM),
                  BestNIM = NA_real_,
                  MaxNIM = x |> arrange(desc(NIM)) |> filter(Performance < 1) |> slice(1) |> pull(NIM),
                  ObsPerf = obs$Performance,
                  ObsNIM = obs$NIM
                )
              }
            }) |>
            bind_rows() ->
            nims

          nimrange <- c(min(nims$MinNIM, na.rm=TRUE), max(nims$MaxNIM, na.rm=TRUE))
        }

        bind_cols(
          pp |> select(ParSet, parasite, drug, endemicity, design, efficacy_expected),
          nims
        ) ->
          out

      }) -> ss

      if(inherits(ss, "try-error")) return(pp |> mutate(Status = "Failed", MSG = as.character(ss)))

      return(out)

    }) ->
    out

  if(any(sapply(out, \(x) "Status" %in% names(x)))){
    warning("One or more error encountered")
    return(out)
  }

  ss <- try({
    out |>
      bind_rows() |>
      ungroup() ->
      out
  })

  return(out)
}


vary_n_analysis <- function(parameters, iters=iterations, performance=NULL, performance_max=0.999, min = 10L, max = 500L, increment=1L, cl=NULL){

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
        #if(is.null(cl)) cat("Parameter cluster ", i, " of ", length(pars), "...\n", sep="")

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
## Calibrate non-inferiority margins
############################################

## NB: non-inferiority margins are without dropouts!
expand_grid(
  parameters_scenario |> distinct(parasite) |> expand_grid(endemicity = seq(5,65,by=2.5)),
  parameters_fixed |> filter(design == "NS_11", min_positive == 1),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "no dropouts", force_inclusion_prob == 0.1) |> mutate(force_inclusion_prob = 0),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |>
      filter(framework=="FHT") |>
      mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

set.seed(2025-07-28)
fn <- "figS1_res.rqs"
if(file.exists(fn)){
  perfout <- qread(fn)
}else{
  ## Takes 40 mins:
  parameters |>
    vary_nim_analysis(performance=0.8, min=100, max=1000, cl=8) ->
    perfout
  qsave(perfout, fn)
}

if(FALSE){
perfout |>
  bind_rows() |>
  select(parasite, drug, endemicity, SampleSize, efficacy_expected, BestNIM) |>
  mutate(Lower = efficacy_expected - BestNIM) |>
  right_join(effort, by = join_by(endemicity, SampleSize)) |>
  group_by(parasite, drug, Effort, efficacy_expected) |>
  summarise(Lower = mean(Lower), .groups='drop') |>
  bind_rows(
    perfout |>
      bind_rows() |>
      mutate(Lower = efficacy_expected - BestNIM) |>
      right_join(effort, by = join_by(endemicity, SampleSize)) |>
      select(parasite, drug, endemicity, Effort, efficacy_expected, Lower)
  ) |>
  mutate(endemicity = if_else(is.na(endemicity), "Average", str_c("End_", endemicity |> format() |> str_replace(" ", "0")))) |>
  bind_rows(
    parameters_thresholds |>
      filter(framework=="FHT") |>
      mutate(Effort=fct("Moderate", levels=levels(effort$Effort)), endemicity = "Using") |>
      select(Effort, parasite, drug, efficacy_expected, Lower=efficacy_lower_target, endemicity)
  ) |>
  mutate(endemicity = fct(endemicity, levels=c(str_c("End_",format(c(5,15,35,65))|>str_replace(" ", "0")), "Average", "Using"))) |>
  pivot_wider(names_from=endemicity, values_from=Lower, names_sort=TRUE) |>
  arrange(parasite, drug, Effort) |>
  mutate(Using = case_when(
    Effort=="Moderate" ~ Using,
    TRUE ~ round(Average, if_else(parasite=="ascaris"&drug=="ALB", 3, 2))
    #TRUE ~ round(End_35, if_else(parasite=="ascaris", 3, 2))
  )) ->
  parameters_all_thresholds_recalc
## Pre-computed/fixed version at the top of this script
writexl::write_xlsx(parameters_all_thresholds, "tableS1_thresholds.xlsx")
}

parameters_all_thresholds |>
  select(parasite, drug, Effort, efficacy_expected, Using) |>
  mutate(Effort = factor(Effort |> as.character(), levels=c("Expected","Extreme","Hard","Moderate","Easy"), labels=c("Expected","Very small","Small","Moderate","Large")), Using=format(Using*100)) |>
  pivot_wider(names_from=Effort, values_from=Using) |>
  arrange(parasite, drug) |>
  writexl::write_xlsx("table_1.xlsx")


### Figures
perfout |>
  mutate(Lower = efficacy_expected - BestNIM) |>
#  filter(Lower >= 1-(1-efficacy_expected)*10, endemicity %in% parameters_scenario[["endemicity"]]) |>
  filter(endemicity %in% parameters_scenario[["endemicity"]]) |>
  arrange(drug, parasite) |>
  mutate(panel = case_when(
    parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
    parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
    parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
  ) |> fct()) ->
  pltdt
pltdt |>
  mutate(Endemicity = fct(str_c(endemicity,"%"))) |>
  ggplot(aes(x=SampleSize, y=Lower*100, col=Endemicity)) +
  #  geom_hline(aes(yintercept = efficacy_expected*100), lty="solid", col="black") +
  #  geom_hline(yintercept = 100, lty="dotted") +
  #  geom_hline(aes(yintercept = 100 * (1-(1-efficacy_expected)*2)), lty="dotted") +
  geom_hline(data=bind_rows(
    parameters_all_thresholds |> select(parasite, drug, Effort, Using),
    parameters_all_thresholds |> mutate(Effort="Expected") |> select(parasite, drug, Effort, Using=efficacy_expected)
  ) |>
    mutate(panel = case_when(
      parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
      parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
      parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
    )) |>
    mutate(Effort = fct(Effort, levels=c("Expected","Extreme","Hard","Moderate","Easy"))),
  aes(yintercept = 100*Using, lty=Effort)) +
  geom_line() +
  scale_x_continuous(breaks=c(100,250,500,1000), minor_breaks=NULL, limits=c(0,1000)) +
  facet_wrap(~ panel, scales="free_y", labeller=label_parsed) +
  scale_linetype_manual(values=c(Expected="solid",Extreme="dotted",Hard="dotdash",Moderate="dashed",Easy="longdash")) ->
  plt
eval(parse(text=str_c(lapply(0:6, \(i){
  if(i==0) return("plt")
  pltdt |>
    distinct(panel, drug, parasite, efficacy_expected) |>
    slice(i) ->
    ii
  parameters_all_thresholds |>
    filter(Effort=="Easy") |>
    distinct(parasite, drug, Effort, Using) |>
    mutate(panel = case_when(
      parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
      parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
      parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
    )) |>
    slice(i) ->
    ii
  ll <- max(0, 100*(1-(1-ii$Using)*1.2))
  if(ii[["panel"]]=="ALB vs. ascaris") ll <- 98.5
  if(ii[["panel"]]=="MEB vs. ascaris") ll <- 90
  if(ii[["panel"]]=="ALB vs. hookworm") ll <- 85
  if(ii[["panel"]]=="MEB vs. hookworm") ll <- 20
  str_c("ggh4x::scale_y_facet(panel=='", ii[["panel"]], "', limits=c(", ll, ", 100))")
}), collapse=" + "))) +
  theme(legend.position="right") +
  guides(linetype = guide_legend(element_blank(), order=1), col = guide_legend(order=2, reverse=TRUE)) +
  ylab("Efficacy (%)") + xlab("Number of Children")


############################################
## New figure S1
############################################

perfout |>
  mutate(Lower = efficacy_expected - BestNIM) |>
  filter(!is.na(Lower)) |>
  #  filter(Lower >= 1-(1-efficacy_expected)*10, endemicity %in% parameters_scenario[["endemicity"]]) |>
  filter(SampleSize %in% c(100,250,500,1000)) |>
  arrange(drug, parasite) |>
  mutate(panel = case_when(
    parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
    parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
    parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
  ) |> fct()) ->
  pltdt
pltdt |>
  mutate(SampleSize = fct(str_c(SampleSize))) |>
  ggplot(aes(x=endemicity, y=Lower*100, col=SampleSize)) +
  #  geom_hline(aes(yintercept = efficacy_expected*100), lty="solid", col="black") +
  #  geom_hline(yintercept = 100, lty="dotted") +
  #  geom_hline(aes(yintercept = 100 * (1-(1-efficacy_expected)*2)), lty="dotted") +
  geom_hline(data=bind_rows(
    parameters_all_thresholds |> select(parasite, drug, Effort, Using),
    #parameters_all_thresholds |> mutate(Effort="Expected") |> select(parasite, drug, Effort, Using=efficacy_expected)
  ) |>
    mutate(panel = case_when(
      parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
      parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
      parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
    ) |> fct()) |>
    mutate(Effort = factor(Effort |> as.character(), levels=c("Expected","Extreme","Hard","Moderate","Easy"), labels=c("Expected","Very small","Small","Moderate","Large"))),
  aes(yintercept = 100*Using, lty=Effort)) +
  geom_hline(data=parameters_all_thresholds |> mutate(Effort="Expected") |> select(parasite, drug, Effort, Using=efficacy_expected) |> mutate(panel = case_when(
    parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
    parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
    parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
  ) |> fct()), aes(yintercept = 100*Using), lty="solid", col="grey50", lwd=2) +
  geom_line() +
  scale_x_continuous(breaks=unique(parameters_scenario[["endemicity"]]), minor_breaks=NULL, limits=c(0,70)) +
  facet_wrap(~ panel, scales="free_y", labeller=label_parsed) +
  scale_linetype_manual(values=c(Expected="solid",`Very small`="dotted",Small="dotdash",Moderate="dashed",Large="longdash")) ->
  plt
eval(parse(text=str_c(lapply(0:6, \(i){
  if(i==0) return("plt")
  pltdt |>
    distinct(panel, drug, parasite, efficacy_expected) |>
    slice(i) ->
    ii
  parameters_all_thresholds |>
    filter(Effort=="Easy") |>
    distinct(parasite, drug, Effort, Using) |>
    arrange(drug, parasite) |>
    mutate(panel = case_when(
      parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
      parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
      parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
    ) |> fct()) |>
    slice(i) ->
    ii
  ll <- max(0, 100*(1-(1-ii$Using)*1.2))
  if(as.numeric(ii[["panel"]])==1) ll <- 98.5
  if(as.numeric(ii[["panel"]])==2) ll <- 85
  if(as.numeric(ii[["panel"]])==4) ll <- 90
  if(as.numeric(ii[["panel"]])==5) ll <- 20
  str_c("ggh4x::scale_y_facet(panel=='", ii[["panel"]], "', limits=c(", ll, ", 100))")
}), collapse=" + "))) +
  theme(legend.position="right") +
  guides(linetype = guide_legend("NIM", order=1), col = guide_legend("SAC", order=2, reverse=TRUE)) +
  ylab("Efficacy (%)") + xlab("Endemicity (%)")
ggsave("figS1.pdf", height=6, width=10)
ggsave("figS1.eps", height=6, width=10)




############################################
## Recreate figure 1
############################################

cols <- c(gg_colour_hue(3),"grey50")
names(cols) <- c("Reduced","Adequate","Inconclusive","Failed")
cols <- c("Reduced" = "#F8766D", "Inconclusive" = "#00BFC4", "Adequate" = "#7CAE00")

parameters_thresholds |>
  distinct(parasite, drug) |>
  left_join(
    parameters_scenario,
    by = join_by(parasite),
    relationship="many-to-many"
  ) |>
  #filter(parasite=="ascaris", drug=="ALB") |>
  expand_grid(
    #n_individ = c(20, 50, 100, 250, 300, 500),
    n_individ = c(100, 500),
    min_positive = c(1)
  ) |>
  filter(parasite=="hookworm" | drug=="ALB", endemicity==15) |>
  # filter(n_individ==250) |>
  identity() ->
  all

## Takes 30 mins:
set.seed(2025-07-28)
fn <- "fig1_res.rqs"
if(file.exists(fn)){
  plots <- qread(fn)
}else{

all |>
  arrange(n_individ) |>
  mutate(Row = row_number()) |>
  rowwise() |>
  group_split() |>
  pblapply(function(x){

    #cat(x$Row, "of", nrow(all), "-", as.numeric(Sys.time()-st, "mins"), "\n")
    #print(x)

    ## NOTE: no dropouts for this!
    expand_grid(
      parameters_fixed |> filter(design == "NS_11", min_positive == x$min_positive),
      parameters_cost |> filter(setting == "Ethiopia"),
      parameters_dropadd |> filter(dropout == "no dropouts") |> mutate(force_inclusion_prob = 0.0),
      parameters_efficacy
    ) |>
      bind_cols(x |> select(-min_positive)) |>
      add_mean_and_cv() |>
      left_join(
        parameters_thresholds_fig1,
        by = join_by(parasite, drug),
        relationship="many-to-many"
      ) |>
      filter({
        (parasite=="hookworm" & drug=="MEB" & true_efficacy >= (0.55-0.05)) |
          (parasite=="hookworm" & drug=="ALB" & true_efficacy >= (0.8-0.05)) |
          (parasite=="ascaris" & drug=="ALB" & true_efficacy >= (0.85-0.05)) |
          (parasite=="trichuris" & drug=="ALB" & true_efficacy >= (0.3-0.05))
      }) |>
      fix_n_analysis(iters=iterations, min = x$n_individ, max=x$n_individ, cl=NULL) ->
      fig_1_data

    fig_1_data |>
      mutate(
        Adequate = if_else(analysis_type=="delta", n_Susceptible, n_above_cutoffs),
        Reduced = if_else(analysis_type=="delta", n_Resistant + n_LowResistant, n_below_cutoffs),
        Inconclusive = if_else(analysis_type=="delta", n_Inconclusive, n_between_cutoffs),
        Failed = n_FailZeroPre + n_FailPositiveScreen + n_FailPositivePre + if_else(analysis_type=="delta", n_ClassifyFail, 0)
      ) |>
      mutate(Total = Failed + Adequate + Reduced + Inconclusive) |>
      select(true_efficacy, efficacy_expected, efficacy_lower_target, analysis=analysis_type, framework, Adequate, Reduced, Inconclusive, Failed) |>
      pivot_longer(Adequate:Failed, names_to="classification", values_to="tally") |>
      mutate(analysis = if_else(analysis=="delta", "hypothesis", analysis)) |>
      bind_cols(x)

  }, cl=8) ->
  plots
qsave(plots, fn)
}


plots |> bind_rows() |> distinct(Row, parasite, drug, n_individ, endemicity)

make_type_factor <- function(x, reorder=TRUE){

  #lvs <- c("mean - WHO", "hypothesis - FHT", "mean - WHO-X", "hypothesis - FHT-X")
  #if(reorder) lvs <- lvs[c(1,3,2,4)]

  lvs <- c("mean - WHO", "mean - WHO-X", "hypothesis - FHT-X", "hypothesis - FHT")
  reorder <- FALSE

  x |>
    mutate(type = str_c(analysis, " - ", framework) |> fct(levels=lvs)) ->
    x

  x |>
    distinct(type, efficacy_lower_target, efficacy_expected) |>
    select(-type) |>
    as.matrix() ->
    vals
  stopifnot(nrow(vals)==length(lvs))
  vals[] <- format(vals*100)

  makelab <- function(char, type, index){
    type <- case_match(type, "WHO" ~ "Point estimate", "FHT" ~ "Hypothesis testing")
    parse(text=str_c('expression(atop("', char, ': ', type, ' with", "T"[l] *"="* " ', vals[index,1], '% T"[u] *"="* " ', vals[index,2], '%"))')) |> eval()
  }

  lin <- 1:4
  if(reorder) lin <- c(1,3,2,4)

  x |>
    mutate(type = factor(type, levels=lvs,
                         labels=c(makelab("A", "WHO", 1),
                                  makelab("B", "WHO", 2),
                                  #makelab("B", "FHT", 2),
                                  #makelab("C", "WHO", 3),
                                  makelab("C", "FHT", 3),
                                  makelab("D", "FHT", 4)
                                  )[lin]
                         ))
                                  #                         expression("A: Point estimate with T"[l] *"="* " 95%"),
                                  #                         expression("C: Point estimate with T"[l] *"="* " 99.6%"),
                                  #                         expression("B: Hypothesis testing with T"[l] *"="* " 95%"),
                                  #                         expression("D: Hypothesis testing with T"[l] *"="* " 99.6%")
}


## Figure 1:
plots |>
  bind_rows() |>
  filter(parasite=="ascaris", drug=="ALB", n_individ==500, endemicity==15) |>
  #filter(n_individ==250) |>
  filter(classification != "Failed") |>
  mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
  group_by(efficacy_expected, analysis, framework, efficacy_lower_target, true_efficacy, n_individ) |>
  arrange(classification) |>
  mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
  ungroup() |>
  make_type_factor() |>
  ggplot(aes(x=true_efficacy*100, ymin=ymin*100, ymax=ymax*100, fill=classification)) +
  geom_ribbon() +
  facet_wrap( ~ type, labeller = label_parsed) +
  # geom_hline(yintercept=c(0.05,0.95)) +
#  geom_segment(aes(y=5, x=50, xend=efficacy_expected*100-10), lty="dotted") +
#  geom_segment(aes(y=95, x=efficacy_expected*100, xend=100), lty="dotted") +
  geom_vline(aes(xintercept=efficacy_expected*100)) +
  geom_vline(aes(xintercept=efficacy_lower_target*100), lty="dashed") +
  theme_minimal() +
  scale_fill_manual(values=cols, guide = guide_legend(reverse = TRUE)) +
  labs(x = "True efficacy (%)",
       y = "Proportion of iterations (%)",
       fill = "Efficacy classification") +
  theme_minimal() +
  # theme(strip.text = element_text(size = 12), legend.title = element_text(size = 12), legend.position = "bottom")
  theme(strip.text = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom") +
  xlim(85,100)
ggsave("fig1.pdf", height=8, width=10)
ggsave("fig1.eps", height=8, width=10)

expand_grid(
  parameters_scenario,
  parameters_fixed,
  parameters_cost |> filter(setting == "Ethiopia")
) |>
  add_mean_and_cv() |>
  distinct(parasite, endemicity, mean_epg)



## New figure S2:
plots |>
  bind_rows() |>
  filter(parasite=="hookworm" | drug=="ALB", endemicity==15) |>
  filter({
    (parasite=="hookworm" & drug=="MEB" & true_efficacy >= 0.55) |
    (parasite=="hookworm" & drug=="ALB" & true_efficacy >= 0.8) |
    (parasite=="ascaris" & drug=="ALB" & true_efficacy >= 0.85) |
    (parasite=="trichuris" & drug=="ALB" & true_efficacy >= 0.3)
  }
  #, true_efficacy >= efficacy_lower_target
  ) |>
  filter(n_individ %in% c(100, 500)) |>
  filter(classification != "Failed") |>
  mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
  group_by(parasite, drug, efficacy_expected, analysis, framework, efficacy_lower_target, true_efficacy, n_individ) |>
  arrange(classification) |>
  mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
  ungroup() |>
  mutate(type = case_when(
    str_c(analysis, " - ", framework) == "mean - WHO" ~ "A",
    str_c(analysis, " - ", framework) == "mean - WHO-X" ~ "B",
    str_c(analysis, " - ", framework) == "hypothesis - FHT-X" ~ "C",
    str_c(analysis, " - ", framework) == "hypothesis - FHT" ~ "D",
  )) |>
  mutate(type_ss = str_c(type, ": ", n_individ, " SAC")) |>
  arrange(parasite, drug) |>
  mutate(drug_parasite = case_when(
    parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
    parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
    parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
  ) |> fct()) |>
  ungroup() |>
  mutate(parasite = if_else(parasite=="hookworm", parasite, str_to_title(parasite))) |>
  ggplot(aes(x=true_efficacy*100, ymin=ymin*100, ymax=ymax*100, fill=classification)) +
  geom_ribbon() +
  facet_grid(type_ss ~ drug_parasite, scales="free_x", labeller = labeller(type_ss = label_value, drug_parasite = label_parsed)) +
#  geom_segment(aes(y=5, x=50, xend=efficacy_lower_target*100), lty="dotted") +
#  geom_segment(aes(y=95, x=efficacy_expected*100, xend=100), lty="dotted") +
  geom_vline(aes(xintercept=efficacy_expected*100)) +
  geom_vline(aes(xintercept=efficacy_lower_target*100), lty="dashed") +
  theme_minimal() +
  scale_fill_manual(values=cols, guide = guide_legend(reverse = TRUE)) +
  labs(x = "True efficacy (%)",
       y = "Proportion of iterations (%)",
       fill = "Efficacy classification") +
  theme_minimal() +
  #theme(strip.text = element_text(size = 12), legend.title = element_text(size = 12))
  theme(strip.text = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom")
ggsave("figS2.pdf", height=10, width=9)
ggsave("figS2.eps", height=10, width=9)


## Old new figure S2 (not using):
plots |>
  bind_rows() |>
  filter(parasite=="hookworm", drug=="MEB", endemicity==15) |>
  #filter(n_individ!=300, n_individ>=100) |>
  filter(classification != "Failed") |>
  mutate(classification = factor(classification, levels=c("Adequate","Inconclusive","Reduced"))) |>
  group_by(efficacy_expected, analysis, framework, efficacy_lower_target, true_efficacy, n_individ) |>
  arrange(classification) |>
  mutate(total = sum(tally), ymax = cumsum(tally/total), ymin = lag(ymax, default=0)) |>
  ungroup() |>
  make_type_factor(reorder=FALSE) |>
  mutate(n_individ = factor(n_individ, levels=c(20,50,100,250,500), labels=c(
    expression("N "*"="*" 20"), expression("N "*"="*" 50"), expression("N "*"="*" 100"), expression("N "*"="*" 250"), expression("N "*"="*" 500")))
  ) |>
  ungroup() |>
  ggplot(aes(x=true_efficacy*100, ymin=ymin*100, ymax=ymax*100, fill=classification)) +
  geom_ribbon() +
  facet_grid(n_individ ~ type, labeller = label_parsed) +
  #  geom_segment(aes(y=5, x=50, xend=efficacy_lower_target*100), lty="dotted") +
  #  geom_segment(aes(y=95, x=efficacy_expected*100, xend=100), lty="dotted") +
  geom_vline(aes(xintercept=efficacy_expected*100)) +
  geom_vline(aes(xintercept=efficacy_lower_target*100), lty="dashed") +
  theme_minimal() +
  scale_fill_manual(values=cols, guide = guide_legend(reverse = TRUE)) +
  labs(x = "True efficacy (%)",
    y = "Proportion of iterations (%)",
    fill = "Efficacy classification") +
  theme_minimal() +
  #theme(strip.text = element_text(size = 12), legend.title = element_text(size = 12))
  theme(strip.text = element_text(size = 12), legend.title = element_blank(), legend.position = "bottom")


############################################
## New figure S3
############################################

set.seed(2025-07-28)
expand_grid(
  parameters_scenario |> filter(parasite=="ascaris", endemicity==15),
  parameters_fixed |> filter(min_positive == 1),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "with dropouts", force_inclusion_prob == 0.1),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) |>
  mutate(min_positive_pre = 1) |>
  fix_n_analysis(iters=iterations, min = 10, max=1000, cl=8) ->
  fig_data
qsave(fig_data, "figS3_res.rqs")
# fig_data <- qread("figS3_res.rqs")

fig_data |>
  mutate(design = fct(design, levels=c("NS_11","NS_12","SSR_11","SSR_12"))) |>
  mutate(cost_sd = sqrt(cost_variance) / 1e3, cost_mean = cost_mean/1e3, Performance=Performance*1e2) |>
  pivot_longer(c(Performance, cost_mean, cost_sd)) |>
  mutate(Panel = factor(
    name,
    levels = c("Performance","cost_mean","cost_sd"),
    labels = c(
      expression("A: Power (" * "%" * ")"),
      expression("B: Mean cost"[total]~ "(x1000 US$)"),
      expression("C: Standard deviation of cost"[total]~ "(x1000 US$)")
    )
  )) ->
  psd

psd |>
  ggplot(aes(x=n_individ, y=value, col=design)) +
  geom_hline(data=tibble(Panel=factor(levels(psd$Panel)[1], levels=levels(psd$Panel)), yi=80), aes(yintercept=yi), lty="dashed") +
#  geom_hline(data=tibble(Panel=factor(levels(psd$Panel)[1], levels=levels(psd$Panel)), yi=90), aes(yintercept=yi), lty="dotted") +
  geom_line() +
  facet_wrap(~Panel, scales="free_y", ncol=1, labeller = label_parsed) +
  coord_cartesian(xlim=c(10,1000)) +
  scale_x_continuous(breaks=c(10,250,500,750,1000)) +
  labs(x="Number of Children", y=NULL) +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1)) +
  ggh4x::scale_y_facet(PANEL==1, limits=c(0,100), breaks=seq(0,100,by=20)) +
  ggh4x::scale_y_facet(PANEL==2, limits=c(0,7.5), breaks=seq(0,8,by=1)) +
  ggh4x::scale_y_facet(PANEL==3, limits=c(0,0.16), breaks=seq(0,0.2,by=0.05))
ggsave("figS3.pdf", width=7, height=6)
ggsave("figS3.eps", width=7, height=6)


############################################
## Re-create figure 2
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="ascaris", endemicity==15),
  parameters_fixed |> filter(min_positive==1),
  parameters_cost,
  parameters_dropadd,
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  filter(
    (dropout == "with dropouts" & force_inclusion_prob == 0.1) | # A/B
      (setting == "Ethiopia" & force_inclusion_prob == 0.1 & dropout == "no dropouts") | # C
      (setting == "Ethiopia" & dropout == "with dropouts" & design=="SSR_12") # D
  ) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

set.seed(2025-07-28)
parameters |>
  vary_n_analysis(cl=8, iters=iterations) ->
  res
qsave(res, "fig2_res.rqs")
# res <- qread("fig2_res.rqs")

LETTERS[1:4] |>
  lapply(\(x){
    if(x=="A"){
      res |>
        filter(setting == "Ethiopia", dropout == "with dropouts", force_inclusion_prob == 0.1) |>
        mutate(Panel = "A: Ethiopian cost")
    }else if(x=="B"){
      res |>
        filter(setting == "Tanzania", dropout == "with dropouts", force_inclusion_prob == 0.1) |>
        mutate(Panel = "B: Tanzanian cost")
    }else if(x=="C"){
      res |>
        filter(setting == "Ethiopia", force_inclusion_prob == 0.1, dropout == "no dropouts") |>
        mutate(Panel = "C: No drop-outs")
    }else if(x=="D"){
      res |>
        filter(setting == "Ethiopia", dropout == "with dropouts", design=="SSR_12", force_inclusion_prob%in%c(0,0.1,0.2)) |>
        mutate(Panel = "D: Assessing multiple STH")
    }else{
      stop("ERROR")
    }
  }) |>
  bind_rows() |>
  #mutate(force_inclusion_prob = case_when(
  #  str_detect(design, "NS") ~ 0.0,
  #  .default = force_inclusion_prob,
  #)) |>
  plot_data_cost() |>
  filter(name=="Performance") |> #, value>0.5, value<0.95) |>
  {function(x){
    return(x)
    x |>
      filter(Panel=="A: Ethiopian cost", design == "NS_12") |>
      mutate(Panel = "D: Assessing multiple STH") |>
      bind_rows(
        x
      )
    #x |>
    #  distinct(design, force_inclusion_prob, Panel) |>
    #  mutate(MeanCost = 70*1e3, value=1) |>
    #  bind_rows(x)
  }}() |>
  mutate(design = fct(design, levels=c("NS_11","NS_12","SSR_11","SSR_12"))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*1e2, col=design, lty=str_c(force_inclusion_prob*100,"%")), parse=TRUE) +
  geom_line() +
  facet_wrap(~Panel, scales="fixed") +
  scale_linetype_manual(values=c(`0%`="dotdash",`10%`="solid",`20%`="dotted")) +
  geom_hline(yintercept = c(80), lty="dashed") +
#  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100), xlim = c(0,7)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Power (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(
    title=bquote(P[add]), order=2, position="inside", theme = theme(legend.position.inside = c(0.92, 0.12)), override.aes = list(col = gg_colour_hue(4)[4])
    ),
    color = guide_legend(title="Survey design", order=1)
  )
ggsave("fig2.pdf", width=7, height=6)
ggsave("fig2.eps", width=7, height=6)


############################################
## Re-create figures 3 and S4
############################################

expand_grid(
  parameters_scenario |> filter(parasite=="ascaris", endemicity!=2),
  parameters_fixed |> filter(min_positive%in%c(1)),
  parameters_cost,
  parameters_dropadd |> filter(force_inclusion_prob == 0.1),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  filter(dropout=="with dropouts" | setting=="Ethiopia") |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

set.seed(2025-07-28)
parameters |>
  vary_n_analysis(cl=8, iters=iterations, increment=1) ->
  res
qsave(res, "fig3_res.rqs")
#res <- qread("fig3_res.rqs")

res |>
  plot_data_cost() |>
  filter(name=="Performance", dropout=="with dropouts", setting=="Ethiopia") |> #, value>0.5, value<0.95) |>
  mutate(Panel = str_c(
    factor(endemicity, levels=c(5,15,35,65), labels=LETTERS[1:4]) |> as.character(),
    ": ", endemicity, "% prevalence; ", round(mean_epg,1), " epg"
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
  geom_hline(yintercept = c(80), lty="dashed") +
#  geom_hline(yintercept = c(90), lty="dotted") +
  scale_x_log10() +
  coord_cartesian(ylim = c(50,100), xlim = c(1.0,30)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Power (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1))
ggsave("fig3.pdf", width=7, height=6)
ggsave("fig3.eps", width=7, height=6)


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
    dropout=="with dropouts" ~ str_c(setting),
    TRUE ~ str_c(setting, " without drop-outs")
  ), levels=c("Ethiopia", "Tanzania", "Ethiopia without drop-outs"))) |>
  mutate(Col = fct(str_c(endemicity,"% prevalence"))) |>
  #filter(dropout=="with dropouts") |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, col=design)) +
  geom_line() +
  facet_grid(Row~Col, scales="fixed") +
  geom_hline(yintercept = c(80), lty="dashed") +
#  geom_hline(yintercept = c(90), lty="dotted") +
  scale_x_log10() +
  coord_cartesian(ylim = c(50,100), xlim = c(1.0,30)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Power (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(lty = guide_legend(title=bquote(P[add]), order=2), color = guide_legend(title="Survey design", order=1))
ggsave("figS4.pdf", width=9, height=7)
ggsave("figS4.eps", width=9, height=7)




############################################
## NOW REMOVED!!!! Re-create figure S4
############################################


expand_grid(
  parameters_scenario |> filter(parasite=="ascaris", endemicity!=2),
  parameters_fixed |> filter(min_positive%in%c(1), design=="SSR_12"),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "with dropouts"),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(drug=="ALB", framework=="FHT") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

set.seed(2025-07-28)
parameters |>
  vary_n_analysis(cl=8, iters=iterations, increment=1) ->
  res
qsave(res, "figOS4_res.rqs")
#res <- qread("figOS4_res.rqs")

res |>
  plot_data_cost() |>
  filter(name=="Performance") |> #, value>0.5, value<0.95) |>
  {function(x){
    x |>
      distinct(design, force_inclusion_prob, endemicity) |>
      mutate(MeanCost = 70*1e3, value=1) |>
      bind_rows(x)
  }}() |>
  mutate(Col = fct(str_c(endemicity,"% prevalence"))) |>
  mutate(Padd = fct(str_c(force_inclusion_prob*100,"%"), levels=str_c(c(0,0.05,0.1,0.15,0.2)*100,"%"))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, lty=Padd)) +
  geom_line(col=gg_colour_hue(4)[4]) +
  facet_wrap(~Col, scales="fixed") +
  ylab("Performance (%)") +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  coord_cartesian(ylim = c(50,100), xlim = c(0,6)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Performance (%)") +
  guides(color = guide_legend(title=bquote(P[add]))) +
  scale_x_continuous(breaks=seq(0,10,by=2))
ggsave("figOS4.pdf", width=7, height=6)



############################################
## Re-create figure 4
############################################

expand_grid(
  parameters_scenario |> filter(endemicity==15),
  parameters_fixed |> filter(min_positive%in%c(1)),
  parameters_cost |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(dropout == "with dropouts", force_inclusion_prob==0.1),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  left_join(
    parameters_thresholds |> filter(framework=="FHT", parasite=="hookworm" | drug=="ALB") |> mutate(true_efficacy = efficacy_expected),
    by = "parasite", relationship="many-to-many"
  ) ->
  parameters

set.seed(2025-07-28)
parameters |>
  vary_n_analysis(cl=8, iters=iterations, increment=1) ->
  res
qsave(res, "fig4_res.rqs")
# res <- qread("fig4_res.rqs")

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
    "ALB against ascaris", "ALB against hookworm", "MEB against hookworm", "ALB against trichuris"
  ), labels=c(
    "ALB~vs~italic(Ascaris)", "ALB~vs.~hookworms", "MEB~vs.~hookworms", "ALB~vs.~italic(Trichuris)"
  ))) |>
  ggplot(aes(x=MeanCost/1e3, y=value*100, col=design)) +
  geom_line() +
  facet_wrap(~Col, scales="free_x", labeller=label_parsed) +
  geom_hline(yintercept = c(80), lty="dashed") +
  geom_hline(yintercept = c(90), lty="dotted") +
  scale_x_log10() +
  coord_cartesian(ylim = c(50,100), xlim = c(1.0,30)) +
  #coord_cartesian(ylim = c(50,100), xlim = c(0, 60)) +
  labs(x=bquote("Mean cost"[total]~ "(x1000 US$)"), y="Power (%)") +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  guides(color = guide_legend(title="Survey design"))
ggsave("fig4.pdf", width=7, height=6)
ggsave("fig4.eps", width=7, height=6)



############################################
## Re-create Table 3 and S2
############################################

expand_grid(
  parameters_scenario,
  parameters_fixed |> filter(min_positive%in%c(1)),
  parameters_cost, # |> filter(setting == "Ethiopia"),
  parameters_dropadd |> filter(force_inclusion_prob==0.1) |> select(starts_with("dropout")),
  parameters_dropadd |> filter(dropout=="with dropouts") |> select(!starts_with("dropout")),
  parameters_analysis |> filter(analysis_type=="delta")
) |>
  add_mean_and_cv() |>
  #  left_join(
  #    parameters_thresholds |>
  #      filter(framework=="FHT") |>
  #      mutate(true_efficacy = efficacy_expected) |>
  #      select(parasite, drug, framework, efficacy_expected, true_efficacy),
  #    by = "parasite", relationship="many-to-many"
  #  ) |>
  left_join(
    parameters_all_thresholds |> select(parasite, drug, Effort, efficacy_expected, efficacy_lower_target=Using),
    by = join_by(parasite),
    relationship="many-to-many"
  ) |>
  mutate(true_efficacy = efficacy_expected) ->
  parameters


## For Table 3:

# Note force_inclusion_prob = 0 is correct
parameters |>
  filter(endemicity==15, Effort=="Moderate", dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") ->
  subp

set.seed(2025-07-28)
subp |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resA
subp |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resB
subp |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
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
qsave(res, "table3_res.rqs")
# res <- qread("table3_res.rqs")

res |>
  filter(Target==0.8) |>
  group_by(drug, parasite, design) |>
  summarise(n_individ = ceiling(mean(n_individ)/5)*5, cost_mean = ceiling(mean(cost_mean)/500)*500, .groups="drop") |>
  mutate(across(c(n_individ,cost_mean), \(x) format(x, big.mark=","))) |>
  # identity()
  #  pivot_wider(names_from="design", values_from="n_individ") |>
  writexl::write_xlsx("table_3.xlsx")


#### Table S2

## Take out p_add for NS then duplicate them later
parameters |>
  filter(design %in% c("SSR_11","SSR_12") | force_inclusion_prob==0) ->
  parameters_subselected

## Takes 15 hours:
set.seed(2025-07-28)
parameters_subselected |>
  slice_sample(prop=1) |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  res
qsave(res, "table3_res.rqs")
# res <- qread("table3_res.rqs")

## Add missing padd and then do group comparisons:
res |>
  filter(design %in% c("NS_11","NS_12")) |>
  select(-force_inclusion_prob) |>
  expand_grid(force_inclusion_prob = unique(res$force_inclusion_prob)) |>
  bind_rows(res |> filter(design %in% c("SSR_11","SSR_12"))) |>
  group_by(endemicity, dropout, force_inclusion_prob, setting, drug, parasite, efficacy_expected, efficacy_lower_target, Effort, Target) |>
  mutate(ndes = n(), n_individ_min = min(n_individ), cost_mean_min = min(cost_mean)) |>
  mutate(n_individ_delta = n_individ-n_individ_min, cost_mean_delta = cost_mean-cost_mean_min) |>
  mutate(n_individ_rel = n_individ-n_individ[design=="NS_11"], cost_mean_rel = cost_mean-cost_mean[design=="NS_11"]) |>
  ungroup() |>
  filter(Target==0.8) |>
  select(drug, parasite, setting, endemicity, dropout, force_inclusion_prob, efficacy_expected, efficacy_lower_target, Effort, ndes, design, n_individ, n_individ_min, n_individ_rel, n_individ_delta, cost_mean, cost_mean_rel, cost_mean_min, cost_mean_delta, cost_variance) |>
  arrange(drug, parasite, setting, endemicity, dropout, force_inclusion_prob, Effort, design) ->
  res

stopifnot(nrow(res)==nrow(parameters), res$ndes==4)

res |>
  mutate(Effort = factor(Effort |> as.character(), levels=c("Expected","Extreme","Hard","Moderate","Easy"), labels=c("Expected","Very small","Small","Moderate","Large"))) |>
  select(-ndes) |>
  rename(NIM = Effort) |>
  select(setting, drug, parasite, Tu = efficacy_expected, endemicity, dropouts=dropout, padd=force_inclusion_prob, NIM, Tl = efficacy_lower_target, design, SAC = n_individ, SACdelta = n_individ_delta, meancost = cost_mean, meancostdelta = cost_mean_delta, varcost = cost_variance) |>
  arrange(setting, drug, parasite, Tu, endemicity, dropouts, padd, NIM, Tl, design) |>
  group_by(setting) |>
  group_split() |>
  set_names(c("ethiopian_costs", "tanzanian_costs")) |>
  writexl::write_xlsx("table_S2.xlsx")


res |>
  #filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  mutate(dropout = factor(dropout, levels=c("baseline","with dropouts"), labels=c("Baseline","With Dropouts"))) |>
  mutate(xloc = fct(case_when(
    design%in%c("NS_11","NS_12") ~ design,
    TRUE ~ str_c(design, ": ", format(force_inclusion_prob))
  ), levels=c("NS_11","sp1","NS_12","sp2",str_c("SSR_11: ", format(seq(0,0.2,by=0.05))),"sp3",str_c("SSR_12: ", format(seq(0,0.2,by=0.05))))) |> as.numeric()) ->
  plotres


### Figure 5

aa <- 0.75

plotres |>
  filter(Effort=="Moderate", setting=="Ethiopia", dropout=="With Dropouts") |>
  mutate(cost_mean = cost_mean/1e3, force_inclusion_prob=str_c(force_inclusion_prob*100,"%") |> fct(levels=str_c(c(0,5,10,15,20),"%"))) |>
  arrange(parasite, drug) |>
  mutate(drug_parasite = case_when(
    parasite=="ascaris" ~ str_c(drug, "~vs.~italic(Ascaris)"),
    parasite=="trichuris" ~ str_c(drug, "~vs.~italic(Trichuris)"),
    parasite=="hookworm" ~ str_c(drug, "~vs.~hookworms"),
  ) |> fct()) |>
  mutate(targets = str_c("(", format(efficacy_lower_target*100), " - ", format(efficacy_expected*100), "%)")) |>
  mutate(endemicity = fct(str_c(endemicity, "% prev."))) |>
  identity() ->
  fig5dat

fig5dat |>
  ggplot(aes(x=force_inclusion_prob, y=cost_mean, col=design, group=design)) +
  #  geom_rect(xmin=0,xmax=2,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[1],alpha=aa) +
  #  geom_rect(xmin=2,xmax=4,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[2],alpha=aa) +
  #  geom_rect(xmin=4,xmax=10,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[3],alpha=aa) +
  #  geom_rect(xmin=10,xmax=16,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[4],alpha=aa) +
#  geom_hline(data=fig5dat |>
#               filter(design=="NS_11") |>
#               distinct(endemicity,drug,parasite,efficacy_lower_target,efficacy_expected,cost_mean,design),
#             aes(yintercept=cost_mean, col=design), lty="solid", alpha=aa) +
#  geom_hline(data=fig5dat |>
#               filter(design=="NS_12") |>
#               distinct(endemicity,drug,parasite,efficacy_lower_target,efficacy_expected,cost_mean,design),
#             aes(yintercept=cost_mean, col=design), lty="dashed", alpha=aa) +
  geom_line(lty="dotted") +
  geom_point() +
  facet_grid(endemicity ~ drug_parasite + targets, scales="free_y", labeller = labeller(endemicity = label_value, drug_parasite = label_parsed, targets = label_value)) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  #  geom_hline(yintercept=0, lty="dashed")+
  #  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
  #  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  # scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  # scale_x_continuous(breaks=seq(1,15,by=2), labels=c(rep(0,2),rep(seq(0,0.2,by=0.1),2))) +
  # scale_color_manual(values=c(Ethiopia="black", Tanzania="grey50")) +
  xlab(bquote("P"[add])) +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  ylab(expression("Mean cost"[total]~ "(x1000 US$) for 80% power")) +
  theme(strip.text.x = element_text(margin = margin(t=1.5, b=1.5)), strip.background = element_rect(fill = "grey70", colour = "grey70"))
ggsave("fig5.pdf", width=12, height=9)
ggsave("fig5.eps", width=12, height=9)


fig5dat |>
  ggplot(aes(x=force_inclusion_prob, y=n_individ, col=design, group=design)) +
  #  geom_rect(xmin=0,xmax=2,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[1],alpha=aa) +
  #  geom_rect(xmin=2,xmax=4,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[2],alpha=aa) +
  #  geom_rect(xmin=4,xmax=10,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[3],alpha=aa) +
  #  geom_rect(xmin=10,xmax=16,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[4],alpha=aa) +
  geom_line(lty="dotted") +
  geom_point() +
  facet_grid(endemicity ~ drug_parasite + targets, scales="free_y", labeller = labeller(endemicity = label_value, drug_parasite = label_parsed, targets = label_value)) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  #  geom_hline(yintercept=0, lty="dashed")+
  #  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
  #  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  # scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  # scale_x_continuous(breaks=seq(1,15,by=2), labels=c(rep(0,2),rep(seq(0,0.2,by=0.1),2))) +
  # scale_color_manual(values=c(Ethiopia="black", Tanzania="grey50")) +
  xlab(bquote("P"[add])) +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  ylab("Required sample size for 80% power") +
  theme(strip.text.x = element_text(margin = margin(t=1.5, b=1.5)), strip.background = element_rect(fill = "grey70", colour = "grey70"))
ggsave("figS5.pdf", width=12, height=9)
ggsave("figS5.eps", width=12, height=9)





stop("FINISHED!")

### OLD figures S5-S7

res |>
  filter(Target == 0.8, force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  mutate(dropout = factor(dropout, levels=c("baseline","with dropouts"), labels=c("Baseline","With Dropouts"))) |>
  mutate(xloc = fct(case_when(
    design%in%c("NS_11","NS_12") ~ design,
    TRUE ~ str_c(design, ": ", format(force_inclusion_prob))
  ), levels=c("NS_11","sp1","NS_12","sp2",str_c("SSR_11: ", format(seq(0,0.2,by=0.05))),"sp3",str_c("SSR_12: ", format(seq(0,0.2,by=0.05))))) |> as.numeric()) ->
  plotres

aa <- 0.01

get_plot <- function(effort,type){
  dt <- plotres |> filter(Effort==effort)
  if(type=="n_individ"){
    pt <- ggplot(dt, aes(x=xloc, y=n_individ, col=setting, pch=dropout, group=str_c(setting,dropout,design)))
  }else if(type=="cost"){
    pt <- ggplot(dt, aes(x=xloc, y=cost_mean, col=setting, pch=dropout, group=str_c(setting,dropout,design)))
  }else if(type=="cost_rel"){
    pt <- ggplot(dt, aes(x=xloc, y=cost_mean_rel, col=setting, pch=dropout, group=str_c(setting,dropout,design)))
  }else{
    stop("Unrecognised type")
  }
  pt +
    geom_rect(xmin=0,xmax=2,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[1],alpha=aa) +
    geom_rect(xmin=2,xmax=4,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[2],alpha=aa) +
    geom_rect(xmin=4,xmax=10,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[3],alpha=aa) +
    geom_rect(xmin=10,xmax=16,ymin=-Inf,ymax=Inf,col="transparent",fill=gg_colour_hue(4)[4],alpha=aa) +
    geom_line() +
    geom_point() +
    facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug, " vs. ", parasite, "\n(", format(efficacy_lower_target*100), " - ", format(efficacy_expected*100), "%)"), scales="free_y") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    #  geom_hline(yintercept=0, lty="dashed")+
    #  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
    #  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
    # scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
    scale_x_continuous(breaks=seq(1,15,by=2), labels=c(rep(0,2),rep(seq(0,0.2,by=0.1),2))) +
    scale_color_manual(values=c(Ethiopia="black", Tanzania="grey50")) +
    xlab(bquote("Survey Design & P"[add]))
}

pdf("notebooks/paper_2025/figS5.pdf", width=12, height=9)
lapply(unique(res$Effort), function(effort){
   get_plot(effort, "n_individ") +
    ylab(str_c("Required sample size for 80% power (", effort, " effort)"))
})
dev.off()

pdf("notebooks/paper_2025/figS6.pdf", width=12, height=9)
lapply(unique(res$Effort), function(effort){
  get_plot(effort, "cost") +
    ylab(str_c("Mean cost for 80% power (", effort, " effort)"))
})
dev.off()

pdf("notebooks/paper_2025/figS7.pdf", width=12, height=9)
lapply(unique(res$Effort), function(effort){
  get_plot(effort, "cost_rel") +
    geom_hline(yintercept=0, lty="solid", col="white") +
    ylab(str_c("Absolute difference in mean cost for 80% power (", effort, " effort)"))
})
dev.off()





##### GRAVEYARD

res |>
  group_by(endemicity, dropout, force_inclusion_prob, setting, drug, parasite, Target) |>
  mutate(cost_mean_min = cost_mean[design=="NS_11"]) |>
  ungroup() |>
  mutate(cost_mean_delta = cost_mean-cost_mean_min) |>
  filter(Target==0.8, design!='NS_11') |>
  filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=cost_mean_delta, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
  #  geom_boxplot() +
  # geom_violin() +
  #  geom_point() +
  geom_point(position = position_dodge(width=0.5)) +
  #  facet_grid(str_c(drug," vs. ", parasite) ~ fct(str_c(endemicity, "% prev.")), scales="fixed") +
  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug, " vs. ", parasite, "\n(", efficacy_lower_target*100, " - ", efficacy_expected*100, ")"), scales="free_y") +
  #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
  labs(y="Absolute difference in mean cost relative to NS_11 (for 80% power)", x="Padd") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=0, lty="dashed")+
  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
#  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
ggsave("notebooks/paper_2025/figS7.pdf", width=12, height=9)





## For Table 3:
parameters |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resA

parameters |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  resB

parameters |>
  filter(endemicity==15, dropout=="with dropouts", force_inclusion_prob==0, setting=="Ethiopia") |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
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
  group_by(drug, parasite, design, Effort, Target) |>
  summarise(n_individ = ceiling(mean(n_individ)/5)*5, .groups="drop") |>
  pivot_wider(names_from="design", values_from="n_individ") |>
  writexl::write_xlsx("notebooks/paper_2025/table_3.xlsx")


## For Table S2:

## Takes around 7.5 hours:
parameters |>
  slice_sample(prop=1) |>
  vary_n_analysis(cl=8L, iters=iterations, performance=c(0.8,0.9), increment=1) ->
  res

res |>
  group_by(endemicity, dropout, force_inclusion_prob, setting, drug, parasite, efficacy_expected, efficacy_lower_target, NIM, Target) |>
  mutate(n_individ_min = min(n_individ), cost_mean_min = min(cost_mean)) |>
  ungroup() |>
  mutate(n_individ_delta = n_individ-n_individ_min, cost_mean_delta = cost_mean-cost_mean_min) |>
  select(drug, parasite, setting, endemicity, dropout, force_inclusion_prob, Target, design, n_individ, n_individ_min, n_individ_delta, cost_mean, cost_mean_min, cost_mean_delta, cost_variance) |>
  arrange(drug, parasite, setting, endemicity, dropout, force_inclusion_prob, Target, design) ->
  res

stopifnot(nrow(res)==(nrow(parameters)*2L))
# qsave(res, "notebooks/paper_2025/tables2_res_t3.rqs")
# res <- qread("notebooks/paper_2025/tables2_res.rqs")


res |>
  filter(cost_mean==cost_mean_min) |>
  count(drug, parasite, endemicity, design) |>
  print(n=Inf)


res |>
  rename(power=Target) |>
  writexl::write_xlsx("notebooks/paper_2025/table_S2.xlsx")


## TODO: just 80% for Figures S5A/B, and remove duplicates padd NS

## Additional plots:

pdf("notebooks/paper_2025/figX_mean_cost.pdf", width=6, height=6)
for(dp in unique(with(res, str_c(drug," vs. ", parasite)))){
  print({
    res |>
      filter(Target==0.8, str_c(drug," vs. ", parasite)==dp) |>
      filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
      ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=cost_mean, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
      #  geom_boxplot() +
      # geom_violin() +
      #  geom_point() +
      geom_point(position = position_dodge(width=0.5)) +
      facet_wrap( ~ fct(str_c(endemicity, "% prev.")), scales="free_y", nrow=2) +
      #  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug," vs. ", parasite) , scales="fixed") +
      #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
      labs(y="Mean cost of survey design for 80% power", x="Padd") +
      theme(legend.position="bottom", legend.title=element_blank()) +
      #  geom_hline(yintercept=0, lty="dashed") +
      ggtitle(dp) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      guides(color = "none", fill = "none")
    #  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
    #  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
  })
}
dev.off()

pdf("notebooks/paper_2025/figX_total_children.pdf", width=6, height=6)
for(dp in unique(with(res, str_c(drug," vs. ", parasite)))){
  print({
    res |>
      filter(Target==0.8, str_c(drug," vs. ", parasite)==dp) |>
      filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
      ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=n_individ, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
      #  geom_boxplot() +
      # geom_violin() +
      #  geom_point() +
      geom_point(position = position_dodge(width=0.5)) +
      facet_wrap( ~ fct(str_c(endemicity, "% prev.")), scales="free_y", nrow=2) +
      #  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug," vs. ", parasite) , scales="fixed") +
      #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
      labs(y="Sample size for 80% power", x="Padd") +
      theme(legend.position="bottom", legend.title=element_blank()) +
      #  geom_hline(yintercept=0, lty="dashed") +
      ggtitle(dp) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      guides(color = "none", fill = "none")
    #  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
    #  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
  })
}
dev.off()

pdf("notebooks/paper_2025/figX_relative_cost.pdf", width=6, height=6)
for(dp in unique(with(res, str_c(drug," vs. ", parasite)))){
  print({
    res |>
      group_by(endemicity, dropout, force_inclusion_prob, setting, drug, parasite, Target) |>
      mutate(cost_mean_min = cost_mean[design=="NS_11"]) |>
      ungroup() |>
      mutate(cost_mean_delta = cost_mean-cost_mean_min) |>
      filter(Target==0.8, design!='NS_11') |>
      filter(Target==0.8, str_c(drug," vs. ", parasite)==dp) |>
      filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
      ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=cost_mean_delta, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
      #  geom_boxplot() +
      # geom_violin() +
      #  geom_point() +
      geom_point(position = position_dodge(width=0.5)) +
      facet_wrap( ~ fct(str_c(endemicity, "% prev.")), scales="free_y", nrow=2) +
      #  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug," vs. ", parasite) , scales="fixed") +
      #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
      labs(y="Absolute difference in cost relative to NS_11 (for 80% power)", x="Padd") +
      theme(legend.position="bottom", legend.title=element_blank()) +
      geom_hline(yintercept=0, lty="dashed") +
      ggtitle(dp) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      guides(col = "none", fill = "none") +
      scale_colour_manual(values=gg_colour_hue(4)[-1])
    #  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
    #  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
  })
}
dev.off()

res |>
  filter(Target==0.8) |>
  filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), design), y=n_individ, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
  #  geom_boxplot() +
  # geom_violin() +
  #  geom_point() +
  geom_point(position = position_dodge(width=0.5)) +
  facet_grid(str_c(drug," vs. ", parasite) ~ fct(str_c(endemicity, "% prev.")), scales="free_y") +
  #  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug," vs. ", parasite) , scales="fixed") +
  #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
  labs(y="Absolute difference in cost relative to NS_11 (for 80% power)", x="Padd") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=0, lty="dashed")



res |>
  filter(Target==0.8) |>
  filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=n_individ, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
  #  geom_boxplot() +
  # geom_violin() +
#  geom_point() +
  geom_point(position = position_dodge(width=0.5)) +
#  facet_grid(str_c(drug," vs. ", parasite) ~ fct(str_c(endemicity, "% prev.")), scales="fixed") +
  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug, " vs. ", parasite, "\n(", efficacy_lower_target*100, " - ", efficacy_expected*100, ")"), scales="free_y") +
#  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
  labs(y="Required sample size for 80% power", x="Padd") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=0, lty="dashed")+
#  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
#  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
ggsave("notebooks/paper_2025/figS5.pdf", width=12, height=9)

res |>
  filter(Target==0.8) |>
  filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=cost_mean, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
  #  geom_boxplot() +
  # geom_violin() +
  #  geom_point() +
  geom_point(position = position_dodge(width=0.5)) +
  #  facet_grid(str_c(drug," vs. ", parasite) ~ fct(str_c(endemicity, "% prev.")), scales="fixed") +
  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug, " vs. ", parasite, "\n(", efficacy_lower_target*100, " - ", efficacy_expected*100, ")"), scales="free_y") +
  #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
  labs(y="Mean cost for 80% power", x="Padd") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=0, lty="dashed")+
  #  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
#  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
ggsave("notebooks/paper_2025/figS6.pdf", width=12, height=9)


res |>
  group_by(endemicity, dropout, force_inclusion_prob, setting, drug, parasite, Target) |>
  mutate(cost_mean_min = cost_mean[design=="NS_11"]) |>
  ungroup() |>
  mutate(cost_mean_delta = cost_mean-cost_mean_min) |>
  filter(Target==0.8, design!='NS_11') |>
  filter(force_inclusion_prob==0 | design%in%c("SSR_11","SSR_12")) |>
  ggplot(aes(x=str_c(format(force_inclusion_prob) |> str_replace(" ", "0"), " (", design, ")"), y=cost_mean_delta, col=design, fill=design, pch=str_c(setting, " ", dropout))) +
  #  geom_boxplot() +
  # geom_violin() +
  #  geom_point() +
  geom_point(position = position_dodge(width=0.5)) +
  #  facet_grid(str_c(drug," vs. ", parasite) ~ fct(str_c(endemicity, "% prev.")), scales="fixed") +
  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug, " vs. ", parasite, "\n(", efficacy_lower_target*100, " - ", efficacy_expected*100, ")"), scales="free_y") +
  #  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
  labs(y="Absolute difference in mean cost relative to NS_11 (for 80% power)", x="Padd") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=0, lty="dashed")+
  scale_colour_manual(values=gg_colour_hue(4)[-1]) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
#  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
ggsave("notebooks/paper_2025/figS7.pdf", width=12, height=9)


res |>
  filter(Target==0.8, force_inclusion_prob==0.2) |>
  ggplot(aes(x=design, y=log10(cost_mean_delta+1), col=design, fill=design)) +
  #  geom_boxplot() +
  geom_violin() +
  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(breaks=c(0,1,2,3,4), labels=c("$0","$10","$100","$1k","$10k")) +
  scale_x_discrete(breaks=NULL) +
  labs(y="Absolute difference in cost", x=NULL) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
ggsave("notebooks/paper_2025/figS5a.pdf", width=9, height=7)

ggplot(res, aes(x=design, y=((cost_mean/cost_mean_min)-1)*100, col=design, fill=design)) +
  geom_boxplot() +
  facet_grid(fct(str_c(endemicity, "% prev.")) ~ str_c(drug," vs. ", parasite), scales="fixed") +
  scale_y_continuous(breaks=c(0,25,50,75,100), labels=c("0%","25%","50%","75%","100%")) +
  scale_x_discrete(breaks=NULL) +
  labs(y="Relative increase in cost", x=NULL) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  scale_colour_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"]))) +
  scale_fill_discrete(labels=c(bquote(NS["1x1/1x1"]),bquote(NS["1x1/1x2"]),bquote(SSR["1x1/1x1"]),bquote(SSR["1x1/1x2"])))
ggsave("notebooks/paper_2025/figS5b.pdf", width=9, height=7)


ggplot(res, aes(x=design, y=cost_mean_delta+1, col=design, fill=design)) +
  geom_violin() +
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
