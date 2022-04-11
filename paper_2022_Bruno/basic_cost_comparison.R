## Only use NS11 (for now)

## Scenario 1:
N <- 100
pre <- rep(20, N)
post <- rep(1, N)

## Scenario 2:
N <- 100
pre <- 51:150
post <- rep(0:5, each=20)

## Use all 3 methods

## Here are the needed parameters:
parameters <- structure(list(method = c("fecpak", "kk", "miniflotac"), count_intercept = c(1.84,  2.38, 2.52), count_coefficient = c(0.172, 0.066, 0.066), time_demography = c(31,  12, 12), time_prep_screen = c(0, 0, 0), time_prep_pre = c(590,  67, 128), time_prep_post = c(590, 67, 128), time_record = c(0,  8, 8), cost_sample = c(0.57, 0.57, 0.57), cost_aliquot_screen = c(0,  0, 0), cost_aliquot_pre = c(1.69, 1.37, 1.51), cost_aliquot_post = c(1.69,  1.37, 1.51), cost_salary = c(22.5, 22.5, 22.5), cost_travel = c(90,  90, 90), n_technicians = c(3, 3, 3), n_team = c(4, 4, 4)), row.names = c(NA,  -3L), class = c("tbl_df", "tbl", "data.frame"))




stop()

# To get parameters:

library("eggSim")
survey_parameters(design="NS_11", parasite="ascaris") %>%
  bind_rows() |>
  select(method, count_intercept, count_coefficient, starts_with("time"), starts_with("cost"), n_technicians, n_team) |>
  dput()

