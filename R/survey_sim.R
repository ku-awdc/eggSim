#' Title
#'
#' @param design
#' @param iterations
#' @param n_individ
#' @param n_day_pre
#' @param n_aliquot_pre
#' @param n_day_post
#' @param n_aliquot_post
#' @param mu_pre
#' @param weight
#' @param performance
#' @param cost_sample
#' @param cost_aliquot
#' @param individ_k
#' @param day_k
#' @param aliquot_k
#' @param efficacy_a
#' @param efficacy_b
#' @param pb
#' @param parameter_output
#'
#' @importFrom pbapply pbsapply
#' @importFrom tidyr expand_grid everything
#' @importFrom dplyr group_by group_split select bind_rows bind_cols
#' @importFrom rlang .data
#'
#' @export
survey_sim <- function(design = "NS", iterations = 10, n_individ=10,
                       mu_pre=10, reduction=0.9,
                       n_samples = data.frame(n_label = "example",
                                              n_day_pre=1, n_aliquot_pre=1,
                                              n_day_post=1, n_aliquot_post=1),
                       weight=1, performance=1,
                       cost_sample=1, cost_aliquot=1,
                       individ_cv=1, day_cv=1,
                       aliquot_cv=1, reduction_cv=1,
                       pb=NA,
                       parameter_output=NA){

  # Allow design (including n_samples) to be replicated for each parameter set:
  parset <- list(n_individ=n_individ, mu_pre=mu_pre, reduction=reduction, weight=weight, performance=performance, cost_sample=cost_sample, cost_aliquot=cost_aliquot, individ_cv=individ_cv, day_cv=day_cv, aliquot_cv=aliquot_cv, reduction_cv=reduction_cv)

  # TODO: nicer error message:
  nparsets <- max(sapply(parset, length))
  stopifnot(all(sapply(parset, length) %in% c(1, nparsets)))

  stop("ALLOW nparsets==iterations AS ALT")

  parset <- as.data.frame(c(parset, list(iteration=1:iterations)))
  stopifnot(all(c("n_label","n_day_pre","n_aliquot_pre","n_day_post","n_aliquot_post") %in% names(n_samples)))
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
      y <- Rcpp_survey_sim_std(x$design, as.data.frame(parset))
    }else{
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


# Old version:
survey_ns <- function(design = "NS", iterations = 10, n_individ=10, n_day_pre=1,
                       n_aliquot_pre=1, n_day_post=1,
                       n_aliquot_post=1, mu_pre=10, weight=1, performance=1,
                       cost_sample=1, cost_aliquot=1, individ_k=1, day_k=1,
                       aliquot_k=1, efficacy_a=1, efficacy_b=1,
                       pb=(length(design)*iterations)>1e4,
                       parameter_output=FALSE){

  if(pb) appfun <- pbapply::pblapply else appfun <- base::lapply

  # Allow design to be replicated for each parameter set:
  parset <- list(Iteration = 1:iterations, n_individ=n_individ, n_day_pre=n_day_pre, n_aliquot_pre=n_aliquot_pre, n_day_post=n_day_post, n_aliquot_post=n_aliquot_post, mu_pre=mu_pre, weight=weight, performance=performance, cost_sample=cost_sample, cost_aliquot=cost_aliquot, individ_k=individ_k, day_k=day_k, aliquot_k=aliquot_k, efficacy_a=efficacy_a, efficacy_b=efficacy_b)
  # TODO: check recycling i.e. all length 1 or iterations

  # TODO: do this in C++ (split and combine has overhead):
  expand_grid(Design = design, as.data.frame(parset)) |>
    group_by(Design, Iteration) |>
    group_split() |>
    appfun(function(x){
      if(x[["Design"]]=="NS"){
        rv <- with(x, Rcpp_survey_ns(as.integer(n_individ), as.integer(n_day_pre), as.integer(n_aliquot_pre), as.integer(n_day_post), as.integer(n_aliquot_post), as.double(mu_pre), as.double(weight), as.double(performance), as.double(cost_sample), as.double(cost_aliquot), as.double(individ_k), as.double(day_k), as.double(aliquot_k), as.double(efficacy_a), as.double(efficacy_b)), simplify=TRUE)
      }else{
        stop("Unrecognised survey design", call.=FALSE)
      }
      bind_cols(data.frame(Efficacy = rv[1], Cost = rv[2]), x)
    }) |>
    bind_rows() |>
    identity() ->
    results

  if(parameter_output){
    results <- results |> select(Design, Iteration, Efficacy, Cost, everything())
  }else{
    results <- results |> select(Design, Iteration, Efficacy, Cost)
  }
  results
}
