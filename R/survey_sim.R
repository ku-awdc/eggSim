#' Simulate one or more survey design
#'
#' @param design survey design(s) to use.  Ignored if parameters and scenario are specified.
#' @param parasite parasite(s) to use.  Ignored if parameters and scenario are specified.
#' @param method method(s) to use.  Ignored if parameters and scenario are specified.
#' @param n_individ number of individuals (can be a vector)
#' @param scenario a data frame of scenarios (see \code{\link{survey_scenario}})
#' @param parameters a list of parameter sets (see \code{\link{survey_parameters}})
#' @param iterations number of iterations per simulation
#' @param cl option to specify either a number of cores or an existing cluster for parallel computation (passed to \code{\link[pbapply]{pblapply}})
#' @param output type of output:  one of summarised, full or extended
#' @param analysis type of ERR/FECRT analysis:  one of mean or delta
#'
#' @importFrom pbapply pblapply
#' @importFrom parallel makeForkCluster makePSOCKcluster stopCluster clusterSetRNGStream clusterExport
#' @importFrom tidyr expand_grid everything
#' @importFrom dplyr group_by group_split select bind_rows bind_cols
#' @importFrom rlang .data
#' @import tidyverse
#'
#' @examples
#' results <- survey_sim(design="NS_11", n_individ = c(100,200,300,400,500))
#'
#' @export
survey_sim <- function(design = c("NS_11","SS_11","SSR_11"),
                       parasite = "hookworm", method = "kk",
                       n_individ = seq(100,1000,by=10),
                       scenario = survey_scenario(parasite),
                       parameters = survey_parameters(design, parasite, method),
                       iterations = 1e3, cl=NULL,
                       output="summarised", analysis="mean"){

  # TODO: pmatching for string arguments
  stopifnot(length(output)==1L, output %in% c("rsummarised","summarised","full","extended"))
  summarise <- output == "summarised"

  stopifnot(length(analysis)==1L, analysis %in% c("mean","delta"))
  if(analysis=="delta" && output=="summarised") stop("Summarised output is not (yet) available with the delta analysis")

  design <- check_design(design)
  parasite <- check_parasite(parasite)
  method <- check_method(method)

  stopifnot(is.numeric(n_individ), all(n_individ > 0L), all(n_individ%%1 == 0L))
  stopifnot(is.numeric(iterations), length(iterations)==1L, iterations > 0L, iterations%%1==0L)

  # Could be specified as a list:
  if(is.data.frame(parameters)) parameters <- list(parameters)
  # TODO: check all needed parameter values are present and non-missing
  check_parameters(parameters, iterations)

  scenario <- check_scenario(scenario)

  ## If cl is an int then set up a cluster:
  if(!is.null(cl) && is.numeric(cl)){
    stopifnot(length(cl)==1L, !is.na(cl), cl > 0L, cl%%1 == 0.0)
    if(.Platform$OS.type=="unix"){
      cl <- makeForkCluster(cl)
    }else{
      cl <- makePSOCKcluster(cl)
    }
    # Required to honour set.seed() on PSOCK clusters:
    clusterSetRNGStream(cl, NULL)
    on.exit({
      stopCluster(cl)
    })
  }

  ## Run the parameter/scenario/n_individ combos:
  cat("Running simulations for ", length(parameters), " parameter sets...\n", sep="")
  st <- Sys.time()

  parameters |>
    # Potentially split different count_parameter sets (should be rare):
    lapply(function(x){
      x |>
        group_split(min_positive_screen, min_positive_pre, count_add, count_mult,
          count_intercept, count_coefficient,
          tail, target_efficacy, target_lower)
    }) |>
    do.call("c", args=_) |>
    # Use of cl argument means we always should use pblapply:
    pblapply(function(x){
      # Remove the design and ns from the parameters:
      x |>
        count(design, n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post, min_positive_screen, min_positive_pre) |>
        select(-n) ->
        all_ns

      # Ensure that N_* are consistent:
      stopifnot(nrow(all_ns)==1L)
      des <- all_ns$design
      all_ns |> select(-design) |> unlist() -> all_ns
      stopifnot(length(all_ns)==8L, all(all_ns>=0L), all(all_ns%%1L == 0L))
      # TODO: nicer error messages

      # Remove count_parameters from the parameters:
      x |>
        count(min_positive_screen, min_positive_pre, count_add, count_mult,
          count_intercept, count_coefficient,
          tail, target_efficacy, target_lower) |>
        select(-n) ->
        count_parameters
      stopifnot(nrow(count_parameters)==1L)

      # Ensure that parameter_set is consistent:
      stopifnot(all(x$parameter_set == x$parameter_set[1L]))

      # Do whatever cost calculations can be done before expanding:
      x |>
        mutate(
          # Note: these costs/times are per aliquot!
          time_consumables_screen = case_when(
            n_aliquot_screen == 0L ~ 0.0,
            TRUE ~ (time_prep_screen + time_record*n_aliquot_screen) / n_aliquot_screen
          ),
          time_consumables_pre = case_when(
            n_aliquot_pre == 0L ~ 0.0,
            TRUE ~ (time_prep_pre + time_record*n_aliquot_pre) / n_aliquot_pre
          ),
          time_consumables_post = case_when(
            n_aliquot_post == 0L ~ 0.0,
            TRUE ~ (time_prep_post + time_record*n_aliquot_post) / n_aliquot_post
          ),

          # Add demography time as a single cost per individual for SSR:
          time_consumables_screen = case_when(
            n_aliquot_screen > 0L ~ time_consumables_screen + time_demography / (n_day_screen*n_aliquot_screen),
            TRUE ~ time_consumables_screen
          ),
          # Add demography time as a single cost per individual for SS and NS:
          time_consumables_pre = case_when(
            n_aliquot_screen == 0L ~ time_consumables_pre + time_demography / (n_day_pre*n_aliquot_pre),
            TRUE ~ time_consumables_pre
          ),

          cost_consumables_screen = case_when(
            n_aliquot_screen == 0L ~ 0.0,
            TRUE ~ (cost_sample + cost_aliquot_screen) / n_aliquot_screen
          ),
          cost_consumables_pre = case_when(
            n_aliquot_pre == 0L ~ 0.0,
            TRUE ~ (cost_sample + cost_aliquot_pre) / n_aliquot_pre
          ),
          cost_consumables_post = case_when(
            n_aliquot_post == 0L ~ 0.0,
            TRUE ~ (cost_sample + cost_aliquot_post) / n_aliquot_post
          )
        ) |>
        # This supercedes the following variables:
        # time_demography, time_prep_*, time_record, cost_sample, cost_aliquot_*
        identity() ->
        x

      # Add a scenario_int that is guaranteed to be 0...S:
      scenario |>
        filter(parasite %in% x$parasite) |>
        mutate(scenario_int = (1:n())-1L) |>
        # Changed parameter name:
        mutate(reduction = 1-true_efficacy) ->
        tscenario

      # Replicate the parameters over iterations (if necessary),
      # and then inner_join with scenario and add replicate ID:
      if(nrow(x)>1L){
        if(nrow(x)!=iterations) stop(paste0("The number of rows of the parameter data frame (", nrow(x), ") does not match the number of iterations (", iterations, ")"))
        if(!"iteration" %in% names(x)) x$iteration <- 1:iterations
        x <- x |> arrange(iteration)
        stopifnot(all(x$iteration == 1:iterations))
      }else{
        x <- bind_cols(x, iteration = 1:iterations)
      }
      x <- inner_join(x, tscenario, by="parasite") |>
        mutate(replicateID = 1:n(), mu_pre = mean_epg * weight * recovery)
      stopifnot(min(x$scenario_int)==0L)

      n_individ <- sort(n_individ)
      stopifnot(all(n_individ > 0L), all(n_individ%%1 == 0))

      if(analysis=="mean"){
        dist_string <- "cs_ga_ga_nb_be_mean"
        if(all(x$aliquot_cv <= 0)) dist_string <- "cs_ga_ga_po_be_mean"
      }else if(analysis=="delta"){
        dist_string <- "cs_ga_ga_nb_be_delta"
        if(all(x$aliquot_cv <= 0)) dist_string <- "cs_ga_ga_po_be_delta"
      }else{
        stop("Unhandled analysis type")
      }

      y <- Rcpp_survey_sim(des, dist_string, all_ns, as.data.frame(x), as.data.frame(count_parameters), n_individ, summarise)

      if(output=="extended"){

        rv <- full_join(y, x, by="replicateID") |>
          select(-replicateID, -scenario_int) |>
          mutate(result = factor(result, levels=c(0,1,2,3), labels=c("Success","FailPositivePre","FailPositiveScreen","ZeroMeanPre"))) |>
          select(-reduction) |>
          select(design, parasite, scenario, mean_epg, true_efficacy, cutoff, method, n_individ, parameter_set, iteration, result, efficacy, total_days, total_cost, everything())
        stopifnot(nrow(rv)==nrow(y))

      }else if(output=="full"){

        rv <- full_join(y,
                        x |> select(design, parasite, cutoff, method, parameter_set, iteration, scenario, mean_epg, true_efficacy, replicateID),
                        by="replicateID"
          ) |>
          mutate(result = factor(result, levels=c(0,1,2,3), labels=c("Success","FailPositivePre","FailPositiveScreen","ZeroMeanPre"))) |>
          select(design, parasite, scenario, mean_epg, true_efficacy, cutoff, method, n_individ, parameter_set, iteration, result, efficacy, total_days, total_cost)
        stopifnot(nrow(rv)==nrow(y))

      }else if(output=="rsummarised"){

        rv <- full_join(y,
                        x |> select(design, parasite, cutoff, method, parameter_set, iteration, scenario, mean_epg, true_efficacy, replicateID),
                        by="replicateID"
          ) |>
          mutate(result = factor(result, levels=c(0,1,2,3), labels=c("Success","FailPositivePre","FailPositiveScreen","ZeroMeanPre"))) |>
          group_by(design, parasite, scenario, mean_epg, true_efficacy, cutoff, method, n_individ, parameter_set) |>

          summarise(below_cutoff = sum(efficacy < cutoff, na.rm=TRUE), above_cutoff = sum(efficacy >= cutoff, na.rm=TRUE), failure_n = sum(result!="Success"), total_n = n(),
            pre_mean = mean(pre_mean, na.rm=TRUE), post_mean = mean(post_mean, na.rm=TRUE),
            pre_imean = mean(pre_imean, na.rm=TRUE), post_imean = mean(post_imean, na.rm=TRUE),
            efficacy_mean = mean(efficacy, na.rm=TRUE), efficacy_variance = var(efficacy, na.rm=TRUE),
            proportion_below = below_cutoff/total_n, cost_mean = mean(total_cost), cost_variance = var(total_cost), .groups="drop") |>
          replace_na(list(below_cutoff = 0L, above_cutoff = 0L))
        with(rv, stopifnot(all((below_cutoff+above_cutoff+failure_n) == total_n)))

      }else if(output=="summarised"){

        rv <- x |>
          count(design, parasite, method, n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post, min_positive_screen, min_positive_pre, scenario, cutoff, mean_epg, true_efficacy, parameter_set, scenario_int) |>
          full_join(y, by="scenario_int") |>
          mutate(efficacy_variance = var_efficacy, failure_n = n_result_1+n_result_2+n_result_3,
                 total_n = n_total, proportion_below = n_below_cutoff/total_n) |>
          select(-scenario_int) |>
          select(design, parasite, scenario, mean_epg, true_efficacy, cutoff, method, n_individ, parameter_set,
                 below_cutoff = n_below_cutoff, above_cutoff = n_above_cutoff, failure_n, total_n,
                 pre_mean = mean_pre, post_mean = mean_post, pre_imean = imean_pre, post_imean = imean_post,
                 efficacy_mean = mean_efficacy, efficacy_variance, proportion_below, days_mean = mean_days, cost_mean = mean_cost, cost_variance = var_cost)
        stopifnot(nrow(y)==nrow(rv), nrow(rv)==(nrow(tscenario)*length(n_individ)))
        with(rv, stopifnot(all((below_cutoff+above_cutoff+failure_n) == total_n)))
      }

      return(rv)
    }, cl=cl) |>
    bind_rows() ->
    results

  cat("Done in ", round(as.numeric(Sys.time()-st, units='mins'), 1), " minutes\n", sep="")

  return(as_tibble(results))
}
