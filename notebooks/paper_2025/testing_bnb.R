survey_sim("NS_11", scenario=survey_scenario("hookworm") |> mutate(true_efficacy=1), output="full", analysis="delta") |> View()
survey_sim("NS_12", scenario=survey_scenario("hookworm") |> mutate(true_efficacy=1), output="full", analysis="delta") |> View()
