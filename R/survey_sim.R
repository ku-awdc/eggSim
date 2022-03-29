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
survey_sim <- function(design = "NS_2x2", parasite = "HW", method = "KK",
                       iterations = 1e3,
                       n_individ=100, mu_pre=100, reduction=0.9,
                       params_design = parameters_design(design),
                       params_parasite = parameters_parasite(parasite),
                       params_method = parameters_method(method),
                       cost = survey_cost(method),
                       pb=NA, parameter_output=NA){

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

  browser()

  stopifnot(length(pb)==1)
  if(is.na(pb)) pb <- length(unique(parset$design_preset)) > 1L
  if(pb) appfun <- pbapply::pblapply else appfun <- base::lapply

  stddsgn <- expand_grid(c("NS","SS","SSR"), c("1x1","1x2","2x1","2x2")) |>
    apply(1,paste,collapse="_")

  parset |>
    group_by(design_preset) |>
    group_split() |>
    appfun(function(x){
      des <- x$design_preset[1]
      if(des %in% stddsgn){
        des <- str_replace(des, "_.*$", "")
        if(des!="NS") stop("Currently only NS is implemented", call.=FALSE)
        # TODO: implement standard specialisations
        # y <- Rcpp_survey_sim_std(x$design, as.data.frame(parset))
        y <- Rcpp_survey_sim_nstd(des, as.data.frame(x))
      }else{
        if(des!="NS") stop("Currently only NS is implemented", call.=FALSE)
        y <- Rcpp_survey_sim_nstd(des, as.data.frame(x))
      }
      bind_cols(x |> select(-iteration), y, parset)
    }) |>
    bind_rows() ->
    results

  browser()

  stopifnot(length(parameter_output)==1)
  if(is.na(parameter_output)) parameter_output <- nparsets>1

  # Calculate costs:
  results$cost <- results$time_count

  if(parameter_output){
    results <- results |> select(design, iteration, efficacy, cost, everything())
  }else{
    if(!all(results$design %in% stddsgn)){
      results <- results |> select(design, n_label, iteration, efficacy, cost)
    }else{
      results <- results |> select(design, iteration, efficacy, cost)
    }
  }
  results
}
