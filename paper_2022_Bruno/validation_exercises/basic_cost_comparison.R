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

pars <- survey_parameters(parasite="ascaris")
scen <- survey_scenario(parasite="ascaris") |> slice(1L)

# Note: recompilation needed for each!
scen1 <- survey_sim(scenario=scen, parameters=pars, n_individ=100L, iterations = 1, output="extended") |> mutate(scenario=1L)
scen2 <- survey_sim(scenario=scen, parameters=pars, n_individ=100L, iterations = 1, output="extended") |> mutate(scenario=2L)
scenarios <- bind_rows(scen1, scen2)
# save(scenarios, file="paper_2022_Bruno/scenarios.rda")

# Costs:
scenarios |>
  select(scenario, method, design, consumables_cost) |>
  mutate(method = factor(method, c("kk","miniflotac","fecpak"))) |>
  mutate(design = factor(design, c("SS_11","SS_12","SSR_11","SSR_12","NS_11","NS_12"))) |>
  arrange(scenario, method, design) ->
  costs
write_csv(costs, file="paper_2022_Bruno/comp_costs.csv")

# Time for counting:
scenarios |>
  mutate(pre = time_screen_count + time_pre_count, post = time_post_count) |>
  mutate(method = factor(method, c("kk","miniflotac","fecpak"))) |>
  filter(design=="SS_11") |>
  arrange(scenario, method, design) |>
  select(scenario, method, pre, post) ->
  counting
write_csv(counting, file="paper_2022_Bruno/comp_counting.csv")


# Time total:
scenarios |>
  mutate(method = factor(method, c("kk","miniflotac","fecpak"))) |>
  mutate(design = factor(design, c("SS_11","SS_12","SSR_11","SSR_12","NS_11","NS_12"))) |>
  arrange(method, design, scenario) |>
  select(method, design, scenario, time_screen, time_screen_count, time_pre, time_pre_count, time_post, time_post_count) ->
  time
write_csv(time, file="paper_2022_Bruno/comp_time.csv")


stop()


# Simple evaluation of bias:

mu <- c(20, 1)
k <- c(0.5, 0.2)
N <- 50
iters <- 1e5

# Note that the following are equivalent:
rnbinom(1, k[1]*N, mu=mu[1]*N)/N
mean(rnbinom(N, k[1], mu=mu[1]))

# So we can use the former to calculate the mean of iters datasets quickly:
pre <- rnbinom(iters, k[1]*N, mu=mu[1]*N)/N
post <- rnbinom(iters, k[2]*N, mu=mu[2]*N)/N

# These are unbiased:
mean(pre) - mu[1]
mean(post) - mu[2]
mean(post)/mean(pre) - mu[2]/mu[1]
# But these are not the same:
mean(post)/mean(pre); mean(post/pre)
# And the latter is biased:
mean(post/pre) - mu[2]/mu[1]

library("tidyverse")
theme_set(theme_light())

# Bias for different values of N:
expand_grid(N = seq(10,1000,by=10), ii=1:it) |>
  mutate(pre = rnbinom(n(), k[1]*N, mu=mu[1]*N)/N) |>
  mutate(post = rnbinom(n(), k[2]*N, mu=mu[2]*N)/N) |>
  group_by(N) |>
  summarise(red1 = mean(post/pre), red2 = mean(post)/mean(pre)) |>
  pivot_longer(c("red1", "red2"), names_to="type", values_to="reduction") |>
  mutate(type = factor(type, levels=c("red1","red2"), labels=c("mean(post/pre)", "mean(post)/mean(pre)"))) ->
  bias

ggplot(bias, aes(x=N, y=100 * (1-reduction), col=type)) +
  geom_line() +
  geom_hline(yintercept=100 * (1-mu[2]/mu[1]), lty="dashed") +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ylab("Efficacy") +
  ggtitle("Plot of bias (black dashed line is true value)")
ggsave("simple_bias.pdf")

mean(mean2)/mean(mean1)
mean(mean2/mean1)



# To get parameters:

library("eggSim")
survey_parameters(design="NS_11", parasite="ascaris") %>%
  bind_rows() |>
  select(method, count_intercept, count_coefficient, starts_with("time"), starts_with("cost"), n_technicians, n_team) |>
  dput()

