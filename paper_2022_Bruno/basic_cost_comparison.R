## Only use NS11 (for now)

## Scenario 1:
N <- 100
pre <- rep(20, N)
post <- rep(1, N)
stopifnot(length(pre)==length(post))

## Scenario 2:
N <- 100
pre <- 51:150
post <- rep(0:4, each=20)
stopifnot(length(pre)==length(post))

## Use all 3 methods

## Here are the needed parameters:
parameters <- structure(list(method = c("fecpak", "kk", "miniflotac"), count_intercept = c(1.84,  2.38, 2.52), count_coefficient = c(0.172, 0.066, 0.066), time_demography = c(31,  12, 12), time_prep_screen = c(0, 0, 0), time_prep_pre = c(590,  67, 128), time_prep_post = c(590, 67, 128), time_record = c(0,  8, 8), cost_sample = c(0.57, 0.57, 0.57), cost_aliquot_screen = c(0,  0, 0), cost_aliquot_pre = c(1.69, 1.37, 1.51), cost_aliquot_post = c(1.69,  1.37, 1.51), cost_salary = c(22.5, 22.5, 22.5), cost_travel = c(90,  90, 90), n_technicians = c(3, 3, 3), n_team = c(4, 4, 4)), row.names = c(NA,  -3L), class = c("tbl_df", "tbl", "data.frame"))


library("tidyverse")


parameters |>
  expand_grid(scenario = 1:2) |>
  group_by(method, scenario) |>
  group_split() |>
  lapply(
    function(x){
      N <- 100
      if(x$scenario==1L){
        pre <- rep(20, N)
        post <- rep(1, N)
      }else if(x$scenario==2L){
        pre <- 51:150
        post <- rep(0:4, each=20)
      }else{
        stop("Unrecognised scenario")
      }

      ctime <- c(
        10^(x$count_intercept + x$count_coefficient*log10(pre+1)^2) |> sum(),
        10^(x$count_intercept + x$count_coefficient*log10(post+1)^2) |> sum()
      )

      stopifnot(x$time_prep_pre==x$time_prep_post)
      ttime <- N * (x$time_demography + x$time_prep_pre + x$time_record) + ctime
      ndays <- ceiling(ttime / (x$n_technicians * 4*60*60)) |> sum()

      cost_s <- ndays*x$n_team*x$cost_salary
      stopifnot(x$cost_aliquot_pre==x$cost_aliquot_post)
      cost_c <- 2L * N * (x$cost_sample + x$cost_aliquot_pre)
      cost_t <- ndays * x$cost_travel

      tibble(method = x$method, scenario = x$scenario, r_cost=cost_s+cost_c+cost_t)
    }
  ) |>
  bind_rows() ->
  r_costs


library("eggSim")

pars <- survey_parameters(design="NS_11", parasite="ascaris")
scen <- survey_scenario(parasite="ascaris") |> slice(1L)

# Note: recompilation needed for each!
c_costs1 <- survey_sim(scenario=scen, parameters=pars, n_individ=100L, iterations = 1) |> mutate(scenario=1L) |> select(method, scenario, c_cost=total_cost)
c_costs1 |> dput()
c_costs2 <- survey_sim(scenario=scen, parameters=pars, n_individ=100L, iterations = 1) |> mutate(scenario=2L) |> select(method, scenario, c_cost=total_cost)
c_costs2 |> dput()

c_costs1 <- structure(list(method = c("fecpak", "kk", "miniflotac"), scenario = c(1L,  1L, 1L), c_cost = c(1172, 748, 1136)), row.names = c(NA, -3L), class = c("tbl_df",  "tbl", "data.frame"))
c_costs2 <- structure(list(method = c("fecpak", "kk", "miniflotac"), scenario = c(2L,  2L, 2L), c_cost = c(1352, 928, 1136)), row.names = c(NA, -3L), class = c("tbl_df",  "tbl", "data.frame"))

all_costs <- cbind(tibble(design="NS_11"), bind_rows(c_costs1, c_costs2)) |> full_join(r_costs)

all_costs

stop()

# To get parameters:

library("eggSim")
survey_parameters(design="NS_11", parasite="ascaris") %>%
  bind_rows() |>
  select(method, count_intercept, count_coefficient, starts_with("time"), starts_with("cost"), n_technicians, n_team) |>
  dput()

