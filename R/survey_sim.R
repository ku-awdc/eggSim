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
#'
#' @export
survey_sim <- function(design = "NS", iterations = 10, n_individ=10, n_day_pre=1,
                       n_aliquot_pre=1, n_day_post=1,
                       n_aliquot_post=1, mu_pre=10, weight=1, performance=1,
                       cost_sample=1, cost_aliquot=1, individ_k=1, day_k=1,
                       aliquot_k=1, efficacy_a=1, efficacy_b=1,
                       pb=length(design)>1,
                       parameter_output=FALSE){

  if(pb) appfun <- pbapply::pblapply else appfun <- base::lapply

  # Allow design to be replicated for each parameter set:
  parset <- list(Iteration = 1:iterations, n_individ=n_individ, n_day_pre=n_day_pre, n_aliquot_pre=n_aliquot_pre, n_day_post=n_day_post, n_aliquot_post=n_aliquot_post, mu_pre=mu_pre, weight=weight, performance=performance, cost_sample=cost_sample, cost_aliquot=cost_aliquot, individ_k=individ_k, day_k=day_k, aliquot_k=aliquot_k, efficacy_a=efficacy_a, efficacy_b=efficacy_b)
  # TODO: check recycling i.e. all length 1 or iterations

  # rv <- appfun(design, function(x) Rcpp_survey_sim("NS", as.data.frame(parset)))
  # return(rv)

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
