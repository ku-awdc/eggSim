#' @name survey_parameters
#' @title Parameter sets for survey simulations
#'
#' @param design survey design(s) to use - must be one of the standard setups
#' @param parasite parasite(s) to use - must be one of the standard options
#' @param method method(s) to use - must be one of the standard options
#'
#' @rdname survey_parameters
#' @export
survey_parameters <- function(design = c("SS_11","SS_12","NS_11","NS_12","SSR_11","SSR_12"),
                              parasite = c("ascaris","trichuris","hookworm"),
                              method = c("kk","miniflotac","fecpak")){

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
      n_aliquot_pre = case_when(
        variant %in% c("11","12") ~ 1L
      ),
      n_day_post = case_when(
        variant %in% c("11","12") ~ 1L
      ),
      n_aliquot_post = case_when(
        variant %in% c("11") ~ 1L,
        variant %in% c("12") ~ 2L
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
      cutoff = case_when(
        parasite == "ascaris" ~ 0.85,
        parasite == "trichuris" ~ 0.40,
        parasite == "hookworm" ~ 0.80
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
      time_demography = case_when(
        method == "kk" ~ 12,
        method == "miniflotac" ~ 12,
        method == "fecpak" ~ 31
      ),
      time_prep_1 = case_when(
        method == "kk" ~ 67,
        method == "miniflotac" ~ 128,
        method == "fecpak" ~ 590
      ),
      time_prep_2 = case_when(
        method == "kk" ~ 133,
        method == "miniflotac" ~ 192,
        method == "fecpak" ~ 1004
      ),
      time_prep_screen = case_when(
        n_aliquot_screen == 0L ~ 0,
        n_aliquot_screen == 1L ~ time_prep_1,
        n_aliquot_screen == 2L ~ time_prep_2
      ),
      time_prep_pre = case_when(
        n_aliquot_pre == 0L ~ 0,
        n_aliquot_pre == 1L ~ time_prep_1,
        n_aliquot_pre == 2L ~ time_prep_2
      ),
      time_prep_post = case_when(
        n_aliquot_post == 0L ~ 0,
        n_aliquot_post == 1L ~ time_prep_1,
        n_aliquot_post == 2L ~ time_prep_2
      ),
      time_record = case_when(
        method == "kk" ~ 8,
        method == "miniflotac" ~ 8,
        method == "fecpak" ~ 0
      ),
      cost_sample = case_when(
        method == "kk" ~ 0.57,
        method == "miniflotac" ~ 0.57,
        method == "fecpak" ~ 0.57
      ),
      cost_aliquot_1 = case_when(
        method == "kk" ~ 1.37,
        method == "miniflotac" ~ 1.51,
        method == "fecpak" ~ 1.69
      ),
      cost_aliquot_2 = case_when(
        method == "kk" ~ 1.51,
        method == "miniflotac" ~ 1.87,
        method == "fecpak" ~ 2.73
      ),
      cost_aliquot_screen = case_when(
        n_aliquot_screen == 0L ~ 0,
        n_aliquot_screen == 1L ~ cost_aliquot_1,
        n_aliquot_screen == 2L ~ cost_aliquot_2
      ),
      cost_aliquot_pre = case_when(
        n_aliquot_pre == 0L ~ 0,
        n_aliquot_pre == 1L ~ cost_aliquot_1,
        n_aliquot_pre == 2L ~ cost_aliquot_2
      ),
      cost_aliquot_post = case_when(
        n_aliquot_post == 0L ~ 0,
        n_aliquot_post == 1L ~ cost_aliquot_1,
        n_aliquot_post == 2L ~ cost_aliquot_2
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
    select(-template, -variant, -cost_aliquot_1, -cost_aliquot_2, -time_prep_1, -time_prep_2) |>
    group_by(design, parasite, method, parameter_set) |>
    group_split() ->
    parameters

  check_parameters(parameters, 1L)
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

  reduction <- tribble(~parasite, ~reduction,
                    "ascaris", beta_to_mu(12.4, 49.5),
                    "trichuris", beta_to_mu(56.1, 30.2),
                    "hookworm", beta_to_mu(17.9, 53.8)
  )

  full_join(mean_epg, reduction, by="parasite") |> filter(.data$parasite %in% parin)

}

beta_to_mu <- function(a, b) (a/(a+b))
beta_to_cv <- function(a, b) sqrt((a*b)/((a+b)^2 * (a+b+1))) / beta_to_mu(a,b)

check_parameters <- function(pp, iters=1L){
  lapply(pp, function(x){
    stopifnot(is.data.frame(x))
    stopifnot(all(!is.na(x)))
    stopifnot(nrow(x) %in% c(1L, iters))
  })
  ps <- sapply(pp, function(x) x$parameter_set)
  stopifnot(all(table(ps)==1L))
  invisible(pp)
}

check_scenario <- function(scenario){
  stopifnot(is.data.frame(scenario))
  stopifnot(all(c("parasite","scenario","mean_epg","reduction") %in% names(scenario)))
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
