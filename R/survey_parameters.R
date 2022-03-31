#' @name survey_parameters
#' @title Parameter sets for survey simulations
#'
#' @param n_label
#' @param n_day_screen
#' @param n_aliquot_screen
#' @param n_day_pre
#' @param n_aliquot_pre
#' @param n_day_post
#' @param n_aliquot_post
#'

#' @rdname survey_parameters
#' @export
parameters_design <- function(design_preset = "custom",
                           n_day_screen=1, n_aliquot_screen=1,
                           n_day_pre=2, n_aliquot_pre=2,
                           n_day_post=2, n_aliquot_post=2
){

  if(FALSE){
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
  }

  rv <- lapply(design_preset, function(x)
  data.frame(design_preset=x, n_day_screen=n_day_screen, n_aliquot_screen=n_aliquot_screen, n_day_pre=n_day_pre, n_aliquot_pre=n_aliquot_pre, n_day_post=n_day_post, n_aliquot_post=n_aliquot_post)
  )

  if(length(rv)==1L) rv <- rv[[1]]

  return(rv)

}


#' @rdname survey_parameters
#' @export
parameters_parasite <- function(parasite_preset = "custom",
                      individ_cv=1, day_cv=1,
                      aliquot_cv=1, reduction_cv=0.1,
                      count_add = 0, count_mult = 1
){

  rv <- lapply(parasite_preset, function(x)
  # TODO: add distributions gamma/lognormal/beta
  data.frame(parasite_preset=x, individ_cv=individ_cv, day_cv=day_cv, aliquot_cv=aliquot_cv, reduction_cv=reduction_cv, count_add=count_add, count_mult=count_mult)
  )

  if(length(rv)==1L) rv <- rv[[1]]

  return(rv)

}


#' @rdname survey_parameters
#' @export
parameters_method <- function(method_preset = "custom",
                       weight=1, performance=1,
                       count_intercept = 2.38961691,
                       count_coefficient = 0.06612497,
                       cost_sample=1, cost_aliquot=1
){

  rv <- lapply(method_preset, function(x)
    data.frame(method_preset=x, weight=weight, performance=performance, count_intercept=count_intercept, count_coefficient=count_coefficient, cost_sample=cost_sample, cost_aliquot=cost_aliquot)
  )

  if(length(rv)==1L) rv <- rv[[1]]

  return(rv)

}
