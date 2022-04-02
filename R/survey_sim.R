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
                       pb=NA, output=NA){

  # output can be without parameters, with parameters or summarised
  # for now only without parameters
  output <- "extended"
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

      stopifnot(nrow(all_ns)==1L)
      des <- all_ns$design
      all_ns |> select(-design) |> unlist() -> all_ns

      stopifnot(length(all_ns)==6L, all(all_ns>=0L), all(all_ns%%1L == 0L))

      # Do whatever cost calculations can be done before expanding:
      x |>
        mutate(
          time_consumables_screen = n_day_screen * (time_demography + time_prep_screen + time_record*n_aliquot_screen),
          time_consumables_pre = n_day_pre * (time_demography + time_prep_pre + time_record*n_aliquot_pre),
          time_consumables_post = n_day_post * (time_demography + time_prep_post + time_record*n_aliquot_post),

          cost_consumables_screen = n_day_screen * (cost_sample + cost_aliquot_screen),
          cost_consumables_pre = n_day_pre * (cost_sample + cost_aliquot_pre),
          cost_consumables_post = n_day_post * (cost_sample + cost_aliquot_post)
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
      y <- Rcpp_survey_sim(des, all_ns, as.data.frame(x), n_individ, summarise)

      if(output=="extended"){
        rv <- full_join(y, x, by="replicateID")
        stopifnot(nrow(rv)==nrow(y))
      }else if(output=="full"){
        rv <- full_join(y,
                        x |> select(design, parasite, method, parameter_set, iteration, scenario, mean_epg, reduction, replicateID),
                        by="replicateID")
        stopifnot(nrow(rv)==nrow(y))
        # Do nothing: just return y
      }else{
        stop("Unimplemented output argument", call.=FALSE)
      }

      return(rv)
    }) |>
    bind_rows() ->
    results

  return(results)

  parameters |>
    mutate(# Note: this supercedes the following variables:
           # time_demography, time_prep_*, time_record, cost_sample, cost_aliquot_*
           # But I will leave them in for now

           time_consumables_screen = n_day_screen * (time_demography + time_prep_screen + time_record*n_aliquot_screen),
           time_consumables_pre = n_day_pre * (time_demography + time_prep_pre + time_record*n_aliquot_pre),
           time_consumables_screen = n_day_post * (time_demography + time_prep_post + time_record*n_aliquot_post),

           cost_consumables_screen = n_day_screen * (cost_sample + cost_aliquot_screen),
           cost_consumables_pre = n_day_pre * (cost_sample + cost_aliquot_pre),
           cost_consumables_post = n_day_post * (cost_sample + cost_aliquot_post)
        )

  # Could be specified as a data frame or a list of data frames:
  if(is.data.frame(params_design)){
    params_design <- list(params_design)
  }else{
    stopifnot(is.list(params_design))
  }
  if(is.data.frame(params_parasite)){
    params_parasite <- list(params_parasite)
  }else{
    stopifnot(is.list(params_parasite))
  }
  if(is.data.frame(params_method)){
    params_method <- list(params_method)
  }else{
    stopifnot(is.list(params_method))
  }

  # We will do every combo of design/parasite/method list specified, so if any
  # design, parasite and/or method have >1 rows then they must have the same rows:
  rows <- sapply(c(params_design,params_parasite,params_method), nrow)

  if(any(!rows %in% c(1L, max(rows)))) stop("The parameter sets for params_design, params_parasite and params_method must have either the same number of rows or 1 row", call.=FALSE)
  if(max(rows)>1L && iterations!=max(rows)) stop("The number of iterations must match the vectorisation over parameter sets for params_design, params_parasite and params_method", call.=FALSE)

  # Expand grid and also by iterations:
  expand_grid(pd = 1:length(params_design), pp = 1:length(params_parasite), pm=1:length(params_method)) |>
    group_by(pd, pp, pm) |>
    group_split() |>
    lapply(function(x){
      bind_cols(params_design[[x$pd]], params_parasite[[x$pp]], params_method[[x$pm]], iteration=1:iterations)
    }) |>
    bind_rows() ->
    parset

  # The remaining parameters can be expand grid'ed:
  expand_grid(n_individ=n_individ, mu_pre=mu_pre, reduction=reduction) |>
    expand_grid(parset) ->
    parset

  stddsgn <- expand_grid(c("NS","SS","SSR"), c("1x1","1x2","2x1","2x2")) |>
    apply(1,paste,collapse="_")

  parset |>
    # group_by(design_preset, parasite_preset, method_preset, n_individ, mu_pre, reduction) |>
    group_by(design_preset, parasite_preset, method_preset) |>
    group_split() ->
    parset

  stopifnot(length(pb)==1)
  if(is.na(pb)) pb <- length(parset) > 1L
  if(pb) appfun <- pbapply::pblapply else appfun <- base::lapply

  parset |>
    appfun(function(x){
      des <- x$design_preset[1]
      if(des %in% stddsgn){
        # TODO: implement standard specialisations
        # y <- Rcpp_survey_sim_std(des, as.data.frame(x))
        des <- str_replace(des, "_.*$", "")
        stopifnot(des%in%c("NS","SS","SSR"))
        y <- Rcpp_survey_sim_nstd(des, as.data.frame(x))
      }else{
        stopifnot(des%in%c("NS","SS","SSR"))
        y <- Rcpp_survey_sim_nstd(des, as.data.frame(x))
      }

      if(output=="full"){
        y <- bind_cols(y, x)
      }else{
        stop("Unimplemented output argument", call.=FALSE)
      }

      return(y)
    }) |>
    bind_rows() ->
    results

  stopifnot(length(parameter_output)==1)
  if(is.na(parameter_output)) parameter_output <- length(params_design)>1

  # Calculate costs:
  results$cost <- results$time_count

  if(parameter_output){
    results <- results |> select(design_preset, iteration, efficacy, cost, everything())
  }else{
    if(!all(results$design_preset %in% stddsgn)){
      results <- results |> select(design_preset, iteration, efficacy, cost)
    }else{
      results <- results |> select(design_preset, iteration, efficacy, cost, n_day_screen, n_aliquot_screen, n_day_pre, n_aliquot_pre, n_day_post, n_aliquot_post)
    }
  }

  return(as_tibble(results))
}
