## Test outputs of survey_sim

n_individ_us <- c(100,200,500)
params <- survey_parameters()
scen <- survey_scenario()
iters <- 1e2
cl <- 2L

set.seed(2022-04-14)
r1 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=cl, output="extended")

set.seed(2022-04-14)
r2 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=cl, output="full")

set.seed(2022-04-14)
r3 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=cl, output="summarised")

test_that("full_is_complete", {
  expect_true(all(names(r2) %in% names(r1)))
})

gni <- 1:9
gpnms <- names(r1)[gni]
test_that("extended_gpnames", {
  expect_true(all(names(r1)[gni] == gpnms))
})
test_that("full_gpnames", {
  expect_true(all(names(r2)[gni] == gpnms))
})
test_that("summarised_gpnames", {
  expect_true(all(names(r3)[gni] == gpnms))
})


