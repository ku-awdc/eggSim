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
survey_sim <- function(design = "NS_2x2", iterations = 1e3, n_individ=100,
                       mu_pre=100, reduction=0.9,
                       n_samples = data.frame(n_label = "example",
                                              n_day_screen=1, n_aliquot_screen=1,
                                              n_day_pre=1, n_aliquot_pre=1,
                                              n_day_post=1, n_aliquot_post=1),
                       weight=1, performance=1,
                       cost_sample=1, cost_aliquot=1,
                       extra_eggs_mult=0, # A multiplier for the cost of counting for other spp
                       extra_eggs_add=0, # An additive (expected) cost of counting for other spp
                       individ_cv=1, day_cv=1,
                       aliquot_cv=1, reduction_cv=0.1,
                       family = "gamma", pb=NA,
                       parameter_output=NA){

  # Allow design (including n_samples) to be replicated for each parameter set:
  parset <- list(n_individ=n_individ, mu_pre=mu_pre, reduction=reduction, weight=weight, performance=performance, cost_sample=cost_sample, cost_aliquot=cost_aliquot, individ_cv=individ_cv, day_cv=day_cv, aliquot_cv=aliquot_cv, reduction_cv=reduction_cv)

  # TODO: nicer error message:
  nparsets <- max(sapply(parset, length))
  stopifnot(all(sapply(parset, length) %in% c(1, nparsets)))

  if(nparsets==iterations){
    parset <- as.data.frame(c(parset, list(iteration=1:iterations)))
  }else{
    parset <- as.data.frame(parset) |> expand_grid(iteration=1:iterations)
  }

  stopifnot(all(c("n_label","n_day_screen","n_aliquot_screen","n_day_pre","n_aliquot_pre","n_day_post","n_aliquot_post") %in% names(n_samples)))
  stopifnot(all(table(n_samples$n_label)==1))

  # Ignore n_samples for standard designs:
  stddsgn <- expand_grid(c("NS","SS","SSR"), c("1x1","1x2","2x1","2x2")) |>
    apply(1,paste,collapse="_")
  designs <- data.frame(design = design[design %in% stddsgn]) |>
    mutate(nums = apply(str_extract_all(.data$design, "[[:digit:]]", simplify=TRUE),1,paste,collapse="_")) |>
    separate(nums, c("n_day_pre", "n_aliquot_pre")) |>
    mutate(n_day_pre=as.numeric(n_day_pre), n_aliquot_pre=as.numeric(n_aliquot_pre)) |>
    mutate(n_day_post=n_day_pre, n_aliquot_post=n_aliquot_pre) |>
    bind_rows(
      expand_grid(design = design[!design %in% stddsgn], n_samples)
    ) |>
    mutate(set = 1:n()) |>
    group_by(set) |>
    group_split()

  stopifnot(length(pb)==1)
  if(is.na(pb)) pb <- length(designs)>1
  if(pb) appfun <- pbapply::pblapply else appfun <- base::lapply

  appfun(designs, function(x){
    if(x$design %in% stddsgn){
      des <- str_replace(x$design, "_.*$", "")
      if(des!="NS") stop("Currently only NS is implemented", call.=FALSE)
      # TODO: implement standard specialisations
      # y <- Rcpp_survey_sim_std(x$design, as.data.frame(parset))
      y <- Rcpp_survey_sim_nstd(des, as.integer(x$n_day_pre), as.integer(x$n_aliquot_pre), as.integer(x$n_day_post), as.integer(x$n_aliquot_post), as.data.frame(parset))
    }else{
      if(x$design!="NS") stop("Currently only NS is implemented", call.=FALSE)
      y <- Rcpp_survey_sim_nstd(x$design, as.integer(x$n_day_pre), as.integer(x$n_aliquot_pre), as.integer(x$n_day_post), as.integer(x$n_aliquot_post), as.data.frame(parset))
    }
    bind_cols(x |> select(-set), y, parset)
  }) |>
    bind_rows() ->
    results

  stopifnot(length(parameter_output)==1)
  if(is.na(parameter_output)) parameter_output <- nparsets>1

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
