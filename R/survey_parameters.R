#' @name survey_parameters
#' @title Parameter sets for survey simulations
#'
#' @param design survey design(s) to use - must be one of the standard setups
#' @param parasite parasite(s) to use - must be one of the standard options
#' @param method method(s) to use - must be one of the standard options
#'
#' @rdname survey_parameters
#'
#' @examples
#' pars <- survey_parameters(design=c("SSR_11","SSR_12"), parasite="ascaris", method="kk")
#' pars
#'
#' @export
survey_parameters <- function(design = c("SS_11","SS_12","NS_11","NS_12","SSR_11","SSR_12"),
                              parasite = c("ascaris","trichuris","hookworm"),
                              method = c("kk","miniflotac","fecpak")){

  if(length(design)==1L && design=="CHECK"){
    ## Hack to allow this to be used to get required parameter names from check_parameters:
    is_check <- TRUE
    design <- "SS_11"
    parasite <- "ascaris"
    method <- "kk"
  }else{
    is_check <- FALSE
  }
  stopifnot(all(design %in% stddsgn))

  design <- check_design(design)
  parasite <- check_parasite(parasite)
  method <- check_method(method)

  # Add in default parameter values:
  expand_grid(template=design, parasite, method) |>
    separate(template, c("design","variant"), remove=FALSE) |>
    mutate(
      parameter_set = str_c(template, parasite, method, sep="_"),
      n_day_screen = case_when(
        design %in% c("NS","SS") ~ 0L,
        variant %in% c("11","12") ~ 1L
      ),
      n_aliquot_screen = case_when(
        design %in% c("NS","SS") ~ 0L,
        variant %in% c("11","12") ~ 1L
      ),
      n_day_pre = case_when(
        variant %in% c("11","12") ~ 1L
      ),
      # Note: SS_12 is different to the others!
      n_aliquot_pre = case_when(
        design == "SS" & variant == "12" ~ 2L,
        variant %in% c("11","12") ~ 1L
      ),
      n_day_post = case_when(
        variant %in% c("11","12") ~ 1L
      ),
      n_aliquot_post = case_when(
        variant %in% c("11") ~ 1L,
        variant %in% c("12") ~ 2L
      ),
      min_positive_screen = case_when(
        design %in% c("NS","SS") ~ 0L,
        design %in% c("SSR") ~ 50L
      ),
      min_positive_pre = case_when(
        design %in% c("NS","SS") ~ 50L,
        design %in% c("SSR") ~ 1L
      ),
      individ_cv = case_when(
        parasite == "ascaris" ~ 3.0,
        parasite == "trichuris" ~ 1.1,
        parasite == "hookworm" ~ 3.0
      ),
      day_cv = case_when(
        parasite == "ascaris" ~ 1.25,
        parasite == "trichuris" ~ 0.75,
        parasite == "hookworm" ~ 1.25
      ),
      aliquot_cv = case_when(
        parasite == "ascaris" & method == "kk" ~ 0.0,
        parasite == "ascaris" & method == "miniflotac" ~ 1/0.579^0.5,
        parasite == "ascaris" & method == "fecpak" ~ 1/0.520^0.5,
        parasite == "trichuris" & method == "kk" ~ 0.0,
        parasite == "trichuris" & method == "miniflotac" ~ 1/3.022^0.5,
        parasite == "trichuris" & method == "fecpak" ~ 1/0.706^0.5,
        parasite == "hookworm" & method == "kk" ~ 0.0,
        parasite == "hookworm" & method == "miniflotac" ~ 1/1.465^0.5,
        parasite == "hookworm" & method == "fecpak" ~ 1/0.574^0.5
      ),
      reduction_cv = case_when(
        # Note: alpha and beta deliberately swapped for compatability with gamma
        parasite == "ascaris" ~ beta_to_cv(12.4, 49.5),
        parasite == "trichuris" ~ beta_to_cv(56.1, 30.2),
        parasite == "hookworm" ~ beta_to_cv(17.9, 53.8)
      ),
      weight = case_when(
        method == "kk" ~ 1/24.0,
        method == "miniflotac" ~ 1/10.0,
        method == "fecpak" ~ 1/34.0
      ),
      recovery = case_when(
        parasite == "ascaris" & method == "kk" ~ 1.0,
        parasite == "ascaris" & method == "miniflotac" ~ 0.645,
        parasite == "ascaris" & method == "fecpak" ~ 0.248,
        parasite == "trichuris" & method == "kk" ~ 1.0,
        parasite == "trichuris" & method == "miniflotac" ~ 1.005,
        parasite == "trichuris" & method == "fecpak" ~ 0.152,
        parasite == "hookworm" & method == "kk" ~ 1.0,
        parasite == "hookworm" & method == "miniflotac" ~ 0.801,
        parasite == "hookworm" & method == "fecpak" ~ 0.569
      ),
      count_add = case_when(
        parasite == "ascaris" ~ 0.0,
        parasite == "trichuris" ~ 0.0,
        parasite == "hookworm" ~ 0.0
      ),
      count_mult = case_when(
        parasite == "ascaris" ~ 1.0,
        parasite == "trichuris" ~ 1.0,
        parasite == "hookworm" ~ 1.0
      ),
#      count_intercept = case_when(
#        method == "kk" ~ 2.38961691, #2.38,
#        method == "miniflotac" ~ 2.51542336, #2.52,
#        method == "fecpak" ~ 1.8348640, #1.84
#      ),
#      count_coefficient = case_when(
#        method == "kk" ~ 0.06612497, #0.066,
#        method == "miniflotac" ~ 0.06609474, #0.066,
#        method == "fecpak" ~ 0.1730919, #0.172
#      ),
      count_intercept = case_when(
        method == "kk" ~ 2.38,
        method == "miniflotac" ~ 2.52,
        method == "fecpak" ~ 1.84
      ),
      count_coefficient = case_when(
        method == "kk" ~ 0.066,
        method == "miniflotac" ~ 0.066,
        method == "fecpak" ~ 0.172
      ),
      tail = 0.025,
      target_efficacy = case_when(
        parasite == "ascaris" ~ 0.95,
        parasite == "trichuris" ~ 0.5,
        parasite == "hookworm" ~ 0.7
      ),
      target_lower = case_when(
        parasite == "ascaris" ~ 0.9,
        parasite == "trichuris" ~ 0.45,
        parasite == "hookworm" ~ 0.65
      ),
      time_demography = case_when(
        method == "kk" ~ 15,
        method == "miniflotac" ~ 15,
        method == "fecpak" ~ 34
      ),
      time_prep_variable = case_when(
        method == "kk" ~ 134-67,
        method == "miniflotac" ~ 197-131,
        method == "fecpak" ~ 1050-596
      ),
      time_prep_fixed = case_when(
        method == "kk" ~ 67 - time_prep_variable,
        method == "miniflotac" ~ 131 - time_prep_variable,
        method == "fecpak" ~ 596 - time_prep_variable
      ),
      time_prep_screen = case_when(
        n_aliquot_screen == 0L ~ 0,
        TRUE ~ time_prep_fixed + time_prep_variable * n_aliquot_screen
      ),
      time_prep_pre = case_when(
        n_aliquot_pre == 0L ~ 0,
        TRUE ~ time_prep_fixed + time_prep_variable * n_aliquot_pre
      ),
      time_prep_post = case_when(
        n_aliquot_post == 0L ~ 0,
        TRUE ~ time_prep_fixed + time_prep_variable * n_aliquot_post
      ),
      time_record = case_when(
        method == "kk" ~ 9,
        method == "miniflotac" ~ 9,
        method == "fecpak" ~ 0
      ),
      cost_sample = case_when(
        method == "kk" ~ 0.57,
        method == "miniflotac" ~ 0.57,
        method == "fecpak" ~ 0.57
      ),
      cost_aliquot_variable = case_when(
        method == "kk" ~ 1.51 - 1.37,
        method == "miniflotac" ~ 1.87 - 1.51,
        method == "fecpak" ~ 2.73 - 1.69
      ),
      cost_aliquot_fixed = case_when(
        method == "kk" ~ 1.37 - cost_aliquot_variable,
        method == "miniflotac" ~ 1.51 - cost_aliquot_variable,
        method == "fecpak" ~ 1.69 - cost_aliquot_variable
      ),
      cost_aliquot_screen = case_when(
        n_aliquot_screen == 0L ~ 0,
        TRUE ~ cost_aliquot_fixed + cost_aliquot_variable * n_aliquot_screen
      ),
      cost_aliquot_pre = case_when(
        n_aliquot_pre == 0L ~ 0,
        TRUE ~ cost_aliquot_fixed + cost_aliquot_variable * n_aliquot_pre
      ),
      cost_aliquot_post = case_when(
        n_aliquot_post == 0L ~ 0,
        TRUE ~ cost_aliquot_fixed + cost_aliquot_variable * n_aliquot_post
      ),
      cost_salary = case_when(
        TRUE ~ 22.50
      ),
      cost_travel = case_when(
        TRUE ~ 90.0
      ),
      # The team is made up of 3 technicians and 1 nurse
      n_technicians = case_when(
        TRUE ~ 3.0
      ),
      n_team = case_when(
        TRUE ~ 4.0
      )
    ) |>
    mutate(design = template) |>
    select(-template, -variant, -cost_aliquot_fixed, -cost_aliquot_variable, -time_prep_fixed, -time_prep_variable) |>
    group_by(design, parasite, method, parameter_set) |>
    group_split() ->
    parameters

  if(!is_check) check_parameters(parameters, 1L)
  if(length(parameters)==1L) parameters <- parameters[[1]]

  return(parameters)
}

#' @rdname survey_parameters
#' @export
survey_scenario <- function(parasite = c("ascaris","trichuris","hookworm")){

  parin <- check_parasite(parasite)

  mean_epg <- tribble(~parasite, ~scenario, ~mean_epg,
                   "ascaris", 1L, 9.6,
                   "ascaris", 2L, 85.2,
                   "ascaris", 3L, 360.0,
                   "ascaris", 4L, 2195.5,
                   "trichuris", 1L, 2.8,
                   "trichuris", 2L, 12.9,
                   "trichuris", 3L, 49.7,
                   "trichuris", 4L, 124.7,
                   "hookworm", 1L, 3.7,
                   "hookworm", 2L, 23.7,
                   "hookworm", 3L, 61.7,
                   "hookworm", 4L, 210.3
  )

  reduction <- tribble(~parasite, ~true_efficacy,
                    "ascaris", 1-beta_to_mu(12.4, 49.5),
                    "trichuris", 1-beta_to_mu(56.1, 30.2),
                    "hookworm", 1-beta_to_mu(17.9, 53.8)
  )

  full_join(mean_epg, reduction, by="parasite") |>
    filter(.data$parasite %in% parin) |>
    mutate(cutoff = case_when(
      parasite == "ascaris" ~ 0.85,
      parasite == "trichuris" ~ 0.40,
      parasite == "hookworm" ~ 0.80
    ))

}

beta_to_mu <- function(a, b) (a/(a+b))
beta_to_cv <- function(a, b) sqrt((a*b)/((a+b)^2 * (a+b+1))) / beta_to_mu(a,b)

check_parameters <- function(pp, iters=1L){
  parnames <- survey_parameters(design="CHECK") |> names()
  lapply(pp, function(x){
    stopifnot(is.data.frame(x))
    stopifnot(all(!is.na(x)))
    stopifnot(nrow(x) %in% c(1L, iters))
    if(any(!parnames %in% names(x))){
      stop("One or more missing parameters: ", str_c(parnames[!parnames %in% names(x)], collapse=", "))
    }
  })
  ps <- lapply(pp, function(x){
    if(!length(unique(x$parameter_set))==1L){
      stop("Invalid parameters:  parameter_set must be consistent", call.=FALSE)
    }
  })
  invisible(pp)
}

check_scenario <- function(scenario){
  stopifnot(is.data.frame(scenario))
  stopifnot(all(c("parasite","scenario","mean_epg","true_efficacy","cutoff") %in% names(scenario)))
  scenario$parasite <- check_parasite(scenario$parasite)
  scenario
}

check_design <- function(x) x

check_parasite <- function(x){
  stopifnot(all(x %in% c("ascaris","hookworm","trichuris")))
  x
}

check_method <- function(x){
  stopifnot(all(x %in% c("kk","miniflotac","fecpak")))
  x
}

stddsgn <- expand_grid(c("NS","SS","SSR"), c("11","12")) |>
  apply(1,paste,collapse="_")
