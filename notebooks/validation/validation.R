## Code to check output between eggSim versions

library("eggSim")

st <- system.time({
sp <- survey_parameters(parasite="hookworm", method="kk")
ss <- survey_scenario(parasite="hookworm")
results <- survey_sim(scenario=ss, parameters=sp)
})

si <- sessionInfo()

results_long <- survey_sim(scenario=ss[1,], parameters=sp[[5]], output="extended")


#save(results, results_long, sp, ss, st, si, file="results.rda")

stop()

rl <- survey_sim(scenario=ss[1,], parameters=sp[[5]], output="extended")
rl |> View()
orl |> View()

rl |>
  group_by(n_individ) |>
#  summarise(cost_mean = mean(consumables_cost+travel_cost+salary_cost), .groups="drop") |>
  summarise(cost_mean = mean(total_cost), .groups="drop") |>
  ggplot(aes(x=n_individ, y=cost_mean)) +
  geom_point()


env <- new.env()
(load("notebooks/validation/validation_results_0.9.5-1_arm.rda", envir=env))
results_old <- env$results

cc <- "parameter_set"
identical(results[[cc]], results_old[[cc]])

cc <- "below_cutoff"
plot(results[[cc]], results_old[[cc]])

names(results)[names(results) %in% names(results_old)] |>
  lapply(function(x){
    if(is.numeric(results_old[[x]]))
      tibble(Parameter=x, Design=results[["design"]], Old=results_old[[x]], New=results[[x]])
    else
      NULL
  }) |>
  bind_rows() |>
  ggplot(aes(x=Old, y=New, col=Design)) +
  geom_abline(slope=1, intercept=0) +
  geom_point() +
  facet_wrap(~Parameter, scales="free")


names(results)[str_detect(names(results), "cost") | str_detect(names(results), "days")] |>
  lapply(function(x){
    if(is.numeric(results_old[[x]]))
      tibble(Parameter=x, Individuals=results[["n_individ"]], Design=results[["design"]], Scenario=results[["scenario"]], Old=results_old[[x]], New=results[[x]])
    else
      NULL
  }) |>
  bind_rows() |>
  filter(Design=="SSR_11") |>
  pivot_longer(c("Old","New")) |>
  ggplot(aes(x=Individuals, y=value, col=name)) +
  #  geom_abline(slope=1, intercept=0) +
  geom_point() +
  facet_wrap(~Parameter+Scenario, scales="free")


names(results)[names(results) %in% names(results_old)] |>
  lapply(function(x){
    if(is.numeric(results_old[[x]]))
      tibble(Parameter=x, Individuals=results[["n_individ"]], Design=results[["design"]], Scenario=results[["scenario"]], Old=results_old[[x]], New=results[[x]])
    else
      NULL
  }) |>
  bind_rows() |>
  pivot_longer(c("Old","New"), names_to="ResultSet", values_to="Estimate") ->
  allres

pdf("comparison.pdf")
allres |>
  group_split(Parameter) |>
  lapply(function(x){
    ggplot(x, aes(x=Individuals, y=Estimate, col=ResultSet)) +
      #  geom_abline(slope=1, intercept=0) +
      geom_point() +
      facet_grid(Design~Scenario) +
      ggtitle(x$Parameter[1])
  })
dev.off()



cst <- eggSim:::CountSummariseTest$new()
cst$add_counts(1:100, 0:99)
cst$result

## TODO: if covariance is 1, method fails
