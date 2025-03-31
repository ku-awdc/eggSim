scenario <- survey_scenario("hookworm") |> mutate(true_efficacy=1)
parameters <- survey_parameters(parasite="hookworm", method="kk", design=c("NS_11","NS_12","SSR_11","SSR_12")) |> bind_rows() |> mutate(min_positive_screen=pmin(1,min_positive_screen), min_positive_pre=1) |> rowwise() |> group_split()

res <- survey_sim(scenario=scenario, parameters=parameters, output="full", analysis="delta", n_individ=seq(10,100,by=10))

ggplot(res, aes(x=n_individ, y=lower_stat, group=n_individ)) +
  geom_boxplot() +
  facet_wrap(~design)# + scale_y_continuous(trans="log10")
