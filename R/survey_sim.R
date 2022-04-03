#' Title
#'
#' @param design
#' @param iterations
#' @param n_individ
#' @param mu_pre
#' @param reduction
#' @param n_samples
#' @param weight
#' @param performance
#' @param cost_sample
#' @param cost_aliquot
#' @param individ_cv
#' @param day_cv
#' @param aliquot_cv
#' @param reduction_cv
#' @param family currently ignored
#' @param pb
#' @param parameter_output
#'
#' @importFrom pbapply pbsapply
#' @importFrom tidyr expand_grid everything
#' @importFrom dplyr group_by group_split select bind_rows bind_cols
#' @importFrom rlang .data
#'
#' @examples
#' results <- survey_sim(design = c("NS_1x1", "NS_2x2"))
#'
#' @export
survey_sim <- function(design = c("NS_11","SS_11","SSR_11"),
                       parasite = "hookworm", method = "kk",
                       n_individ = c(100, 200, 1000),
                       scenario = survey_scenario(parasite),
                       parameters = survey_parameters(design, parasite, method),
                       iterations = 1e3,
                       pb=NA, output="full"){

  # TODO: pmatching for string arguments
  stopifnot(length(output)==1L, output %in% c("summarised","full","extended"))
  # Disable summarised output option for now:
  stopifnot(output!="summarised")
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

  # TODO: possibility for parallel apply?
  stopifnot(is.logical(pb), length(pb)==1L)
  if(is.na(pb)) pb <- length(parameters) > 1L
  if(pb) appfun <- pbapply::pblapply else appfun <- base::lapply


  ## Run the parameter/scenario/n_individ combos:
  parameters |>
    appfun(function(x){

      # Remove the design and ns from the parameters:
      x |>
        count(design, n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post) |>
        select(-n) ->
        all_ns

      # Ensure that N_* are consistent:
      stopifnot(nrow(all_ns)==1L)
      des <- all_ns$design
      all_ns |> select(-design) |> unlist() -> all_ns
      stopifnot(length(all_ns)==6L, all(all_ns>=0L), all(all_ns%%1L == 0L))
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

      # Replicate the parameters over iterations (if necessary),
      # and then inner_join with scenario and add replicate ID:
      if(nrow(x)>1L){
        stopifnot("iteration" %in% names(x))
        x <- x |> arrange(iteration)
        stopifnot(all(x$iteration == 1:iterations))
      }else{
        x <- bind_cols(x, iteration = 1:iterations)
      }
      x <- inner_join(x, scenario, by="parasite") |>
        mutate(replicateID = 1:n(), mu_pre = mean_epg * weight * recovery)

      n_individ <- sort(n_individ)
      stopifnot(all(n_individ > 0L), all(n_individ%%1 == 0))

      # TODO: don't expand parameters by iteration unless needed??

      stopifnot(!summarise)
      all_dists <- c(individ="rgamma", day="rgamma", aliquot="rnbinom", reduction="rbeta")
      if(all(x$aliquot_cv <= 0)) all_dists["aliquot"] <- "rpois"

      y <- Rcpp_survey_sim(des, all_dists, all_ns, as.data.frame(x), n_individ, summarise)

      if(output=="extended"){
        rv <- full_join(y, x, by="replicateID")
        stopifnot(nrow(rv)==nrow(y))
      }else if(output=="full"){
        rv <- full_join(y,
                        x |> select(design, parasite, method, parameter_set, iteration, scenario, mean_epg, reduction, replicateID),
                        by="replicateID")
        stopifnot(nrow(rv)==nrow(y))
      }else{
        stop("Unimplemented output argument", call.=FALSE)
      }

      return(rv)
    }) |>
    bind_rows() ->
    results

  return(as_tibble(results))
}
