#' Simulate one or more survey design
#'
#' @param design survey design(s) to use.  Ignored if parameters and scenario are specified.
#' @param parasite parasite(s) to use.  Ignored if parameters and scenario are specified.
#' @param method method(s) to use.  Ignored if parameters and scenario are specified.
#' @param n_individ number of individuals (can be a vector)
#' @param scenario a data frame of scenarios (see \code{\link{survey_scenario}})
#' @param parameters a list of parameter sets (see \code{\link{survey_parameters}})
#' @param iterations number of iterations per simulation
#' @param cl option to specify a cluster for parallel computation
#' @param output type of output:  one of summarised, full or extended
#'
#' @importFrom pbapply pblapply
#' @importFrom tidyr expand_grid everything
#' @importFrom dplyr group_by group_split select bind_rows bind_cols
#' @importFrom rlang .data
#'
#' @examples
#' results <- survey_sim()
#'
#' @export
survey_sim <- function(design = c("NS_11","SS_11","SSR_11"),
                       parasite = "hookworm", method = "kk",
                       n_individ = c(100, 200, 1000),
                       scenario = survey_scenario(parasite),
                       parameters = survey_parameters(design, parasite, method),
                       iterations = 1e3,
                       cl=NULL, output="summarised"){

  # TODO: pmatching for string arguments
  stopifnot(length(output)==1L, output %in% c("rsummarised","summarised","full","extended"))
  summarise <- output == "summarised"

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

  ## Run the parameter/scenario/n_individ combos:
  parameters |>
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

      # Do whatever cost calculations can be done before expanding:
      x |>
        mutate(
          # Note: these costs/times are per aliquot!
          time_consumables_screen = case_when(
            n_aliquot_screen == 0L ~ 0.0,
            TRUE ~ (time_demography + time_prep_screen + time_record*n_aliquot_screen) / n_aliquot_screen
          ),
          time_consumables_pre = case_when(
            n_aliquot_pre == 0L ~ 0.0,
            TRUE ~ (time_demography + time_prep_pre + time_record*n_aliquot_pre) / n_aliquot_pre
          ),
          time_consumables_post = case_when(
            n_aliquot_post == 0L ~ 0.0,
            TRUE ~ (time_demography + time_prep_post + time_record*n_aliquot_post) / n_aliquot_post
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
        mutate(scenario_int = (1:n())-1L) ->
        tscenario

      # Replicate the parameters over iterations (if necessary),
      # and then inner_join with scenario and add replicate ID:
      if(nrow(x)>1L){
        stopifnot("iteration" %in% names(x))
        x <- x |> arrange(iteration)
        stopifnot(all(x$iteration == 1:iterations))
      }else{
        x <- bind_cols(x, iteration = 1:iterations)
      }
      x <- inner_join(x, tscenario, by="parasite") |>
        mutate(replicateID = 1:n(), mu_pre = mean_epg * weight * recovery)

      n_individ <- sort(n_individ)
      stopifnot(all(n_individ > 0L), all(n_individ%%1 == 0))

      # TODO: don't expand parameters by iteration unless needed??

      dist_string <- "cs_ga_ga_nb_be"
      if(all(x$aliquot_cv <= 0)) dist_string <- "cs_ga_ga_po_be"

      y <- Rcpp_survey_sim(des, dist_string, all_ns, as.data.frame(x), n_individ, summarise)

      if(output=="extended"){
        rv <- full_join(y, x, by="replicateID") |>
          select(-replicateID, -scenario_int) |>
          mutate(result = factor(result, levels=c(0,1,2,3), labels=c("Success","FailPositivePre","FailPositiveScreen","ZeroMeanPre"))) |>
          select(design, parasite, scenario, mean_epg, reduction, method, n_individ, parameter_set, iteration, result, efficacy, total_cost, everything())
        stopifnot(nrow(rv)==nrow(y))
      }else if(output=="full"){
        rv <- full_join(y,
                        x |> select(design, parasite, cutoff, method, parameter_set, iteration, scenario, mean_epg, reduction, replicateID),
                        by="replicateID"
          ) |>
          mutate(result = factor(result, levels=c(0,1,2,3), labels=c("Success","FailPositivePre","FailPositiveScreen","ZeroMeanPre"))) |>
          select(design, parasite, scenario, mean_epg, reduction, cutoff, method, n_individ, parameter_set, iteration, result, efficacy, total_cost)
        stopifnot(nrow(rv)==nrow(y))
      }else if(output=="rsummarised"){
        rv <- full_join(y,
                        x |> select(design, parasite, cutoff, method, parameter_set, iteration, scenario, mean_epg, reduction, replicateID),
                        by="replicateID"
          ) |>
          mutate(result = factor(result, levels=c(0,1,2,3), labels=c("Success","FailPositivePre","FailPositiveScreen","ZeroMeanPre"))) |>
          group_by(design, parasite, scenario, mean_epg, reduction, method, n_individ, parameter_set) |>
          summarise(below_cutoff = sum(efficacy < cutoff, na.rm=TRUE), above_cutoff = sum(efficacy >= cutoff, na.rm=TRUE), failure = sum(result!="Success"), total_n = n(),
            efficacy_mean = mean(efficacy, na.rm=TRUE), efficacy_precision = 1/var(efficacy, na.rm=TRUE), bias_mean = mean(efficacy/(1-reduction), na.rm=TRUE),
            proportion_below = below_cutoff/total_n, cost_mean = mean(total_cost), .groups="drop") |>
          replace_na(list(below_cutoff = 0L, above_cutoff = 0L)) |>
          mutate(efficacy_mean = case_when((below_cutoff+above_cutoff)<10L ~ NA_real_, TRUE ~ efficacy_mean)) |>
          mutate(efficacy_precision = case_when((below_cutoff+above_cutoff)<10L ~ NA_real_, TRUE ~ efficacy_precision))
        with(rv, stopifnot(all(abs((below_cutoff+above_cutoff+failure) == total_n))))

      }else if(output=="summarised"){
        rv <- x |>
          count(design, parasite, method, n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post, min_positive_screen, min_positive_pre, scenario, cutoff, mean_epg, reduction, scenario_int) |>
          full_join(y, by="scenario_int") |>
          select(-scenario_int)
        stopifnot(nrow(y)==nrow(rv), nrow(rv)==(nrow(tscenario)*length(n_individ)))
      }

      return(rv)
    }, cl=cl) |>
    bind_rows() ->
    results

  warning("TODO: Tidy up names of summarised and rsummarised and write test to check summarised outputs are the same with same PRNG")
  warning("TODO: Move cutoff to scenario rather than parameters")
  warning("TODO: allow cl to specify number of cores")

  return(as_tibble(results))
}
