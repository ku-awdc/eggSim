### Relationship between power and expected efficacy

library("eggSim")
library("tidyverse")
library("pbapply")
library("TeachingDemos")

theme_set(theme_light())

basepars <- survey_parameters(design="NS_11", parasite="trichuris", method="kk")

expand_grid(
  Target = seq(0.75,0.95,by=0.01),
  DeltaType = c("absolute","relative"),
) |>
  mutate(
    Delta = case_when(
      DeltaType == "absolute" ~ 0.05,
      DeltaType == "relative" ~ (1-Target)*0.25,
    ),
    Lower = Target-Delta
  ) |>
  rowwise() |>
  group_split() |>
  pblapply(function(x){
    basepars |>
      bind_cols(x) |>
      mutate(efficacy_expected = Target, efficacy_lower_target=Lower) ->
      pars

    scen <- survey_scenario("trichuris") |> mutate(true_efficacy = pars$Target)

    capture.output(res <- survey_sim(n_individ = 1e3, parameters = pars, scenario = scen, analysis="delta"))
    res |> bind_cols(x)
  }, cl=10L) |>
  bind_rows() ->
  results

add_ci <- function(x){
  x |>
    rowwise() |>
    group_split() |>
    lapply(function(y){
      ci <- with(y, hpd(qbeta, shape1=n_LowResistant+n_Susceptible+1, shape2=n_total-(n_LowResistant+n_Susceptible)+1))
      y |> mutate(Lower = ci[1], Upper=ci[2])
    }) |>
    bind_rows()
}

results |>
  mutate(Power = (n_LowResistant+n_Susceptible)/n_total, MeanEPG=str_c("MeanEPG: ", format(mean_epg))) |>
  add_ci() |>
  select(Target, MeanEPG, DeltaType, Delta, Lower, Upper, Power) |>
  ggplot(aes(x=Target, y=Power, ymin=Lower, ymax=Upper, col=DeltaType, fill=DeltaType)) +
  geom_ribbon(alpha=0.5, col="transparent") +
  geom_line() +
  facet_wrap(~MeanEPG) +
  xlab("Target Efficacy") +
  ylab("Power (non-inferiority)")
ggsave("powerplot.pdf")
