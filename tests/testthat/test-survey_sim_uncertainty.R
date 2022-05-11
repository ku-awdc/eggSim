## Test ability to integrate over uncertainty in input parameters

library("eggSim")

n_individ_us <- c(100,200,500)
params <- survey_parameters()[1:2]
scen <- survey_scenario()
iters <- 1e2
cl <- 1L

# Change one of cvs to have uncertainty:
set.seed(2022-04-14)
params[[2]] <- params[[2]] |> slice(rep(1, iters)) |> mutate(reduction_cv = rgamma(iters, 1, 1))

# Should succeed:
test_that("iterating over uncertainty succeeds", {
  r1 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=cl, output="summarised")
  expect_equal(nrow(r1), length(n_individ_us)*4*2)
})

# Should fail:
test_that("iteration argument is checked", {
  expect_error(survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters*2, cl=cl, output="summarised"))
})

# Check we have some correlation with costs over iterations:
params[[2]] <- params[[2]] |> mutate(cost_salary = 1:n())
test_that("there is correlation in input and output cost", {
  set.seed(2022-04-14)
  r2 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params[[2]], iterations=iters, cl=cl, output="extended")
  expect_true(with(r2, cor.test(total_cost, cost_salary))$estimate > 0.25)
})
