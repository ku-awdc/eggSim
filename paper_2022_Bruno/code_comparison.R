## Comparison of code bases for Bruno's paper

library("tidyverse")

list(
  read_csv("paper_2022_Bruno/data_al_17DEC2021_POIS.csv") |> mutate(parasite = "ascaris"),
  read_csv("paper_2022_Bruno/data_hw_17DEC2021_POIS.csv") |> mutate(parasite = "hookworm"),
  read_csv("paper_2022_Bruno/data_tri_17DEC2021_POIS.csv") |> mutate(parasite = "trichuris")
) |>
  bind_rows() |>
  select(parasite, method, mean_epg=mu_pop, n_individ=n, starts_with(c("power","cost","miss"))) |>
  pivot_longer(cols = starts_with(c("power","cost","miss")), names_to=c("type"), values_to="value") |>
  filter(!grepl("_", type)) |>
  mutate(type = gsub("power","power_",type)) |>
  mutate(type = gsub("cost","cost_",type)) |>
  mutate(type = gsub("miss","miss_",type)) |>
  separate(type, c("type", "design"), sep="_") |>
  pivot_wider(names_from="type", values_from="value") |>
  mutate(code="bruno") |>
  select(code, design, everything()) ->
  bruno

library("eggSim")
n_individ_us <- unique(bruno$n_individ)
#n_individ_us <- c(100,200,500,1000)
#n_individ_us <- 100

params <- survey_parameters()
scen <- survey_scenario()

system.time({
  results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=1e3, output="extended")
})

results %>% count(design,result)


if(FALSE){
params <- survey_parameters("NS_11", parasite="ascaris", method="kk")
params$design <- "NS"
scen <- survey_scenario()

system.time({
results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=1e4, output="extended")
})

params1 <- params[[1]]
params2 <- params[[1]]
params2$design <- "NS"
system.time({
  results2 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params2, iterations=1e4, output="extended")
})
system.time({
  results1 <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params1, iterations=1e4, output="extended")
})
results1 |>
  group_by(design, parasite, method, n_individ, scenario, mean_epg, reduction) |>
  summarise(nonmiss = sum(!is.na(efficacy)), proportion = sum(efficacy < cutoff, na.rm=TRUE) / nonmiss, precision = 1/var(efficacy, na.rm=TRUE), miss = 1-(nonmiss/n()), cost = mean(total_cost), .groups="drop") |>
  mutate(design = str_replace(design, "_", ""), code = "matt") |>
  mutate(method = factor(method, levels=c("kk", "fecpak", "miniflotac"), labels=c("KK","FP","MF")) |> as.character()) ->
  matt1
results2 |>
  group_by(design, parasite, method, n_individ, scenario, mean_epg, reduction) |>
  summarise(nonmiss = sum(!is.na(efficacy)), proportion = sum(efficacy < cutoff, na.rm=TRUE) / nonmiss, precision = 1/var(efficacy, na.rm=TRUE), miss = 1-(nonmiss/n()), cost = mean(total_cost), .groups="drop") |>
  mutate(design = str_replace(design, "_", ""), code = "matt") |>
  mutate(method = factor(method, levels=c("kk", "fecpak", "miniflotac"), labels=c("KK","FP","MF")) |> as.character()) ->
  matt2
plot(matt1$proportion, matt2$proportion); abline(0,1)
}


results |>
  group_by(design, parasite, method, n_individ, scenario, mean_epg, reduction) |>
  summarise(nonmiss = sum(!is.na(efficacy)), proportion = sum(efficacy < cutoff, na.rm=TRUE) / n(), precision = 1/var(efficacy, na.rm=TRUE), miss = 1-(nonmiss/n()), cost = mean(total_cost), .groups="drop") |>
  mutate(design = str_replace(design, "_", ""), code = "matt") |>
  mutate(method = factor(method, levels=c("kk", "fecpak", "miniflotac"), labels=c("KK","FP","MF")) |> as.character()) ->
  matt


both <- inner_join(matt |> rename(matt = "proportion") |> select(-code, -cost, -miss),
                  bruno |> rename(bruno = "power") |> select(-code, -cost, -miss)
)
theme_set(theme_light())
ggplot(both, aes(x=bruno, y=matt, group = mean_epg, col=interaction(n_individ, mean_epg, drop=TRUE, sep=" / "))) +
  geom_line() +
  geom_abline(slope=1, intercept=0) +
  geom_point() +
  facet_wrap(~ design + parasite + method, ncol=9) +
  theme(legend.position = "none")
ggsave("comparison_proportion.pdf", width=15, height=15)


both <- inner_join(matt |> rename(matt = "cost") |> select(-code, -miss),
                   bruno |> rename(bruno = "cost") |> select(-code, -miss)
)
theme_set(theme_light())
ggplot(both, aes(x=bruno/1e3, y=matt/1e3, group = mean_epg, col=interaction(n_individ, mean_epg, drop=TRUE, sep=" / "))) +
  geom_line() +
  geom_abline(slope=1, intercept=0) +
  geom_point() +
  facet_wrap(~ design + parasite + method, ncol=9) +
  theme(legend.position = "none")
ggsave("comparison_cost.pdf", width=15, height=15)


## TODO: some of the discrepancy may be due to how pre-mean==0 is handled?

bruno$proportion <- bruno$power

pdf("all_precision.pdf")
for(pp in unique(matt$parasite)){
  for(pn in c("matt")){
    dt <- get(pn)
    pt <- ggplot(dt |> filter(parasite==pp), aes(x=n_individ, y=precision, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn))
    print(pt)
  }
}
dev.off()

pdf("all_tradeoff_po.pdf")
for(pp in unique(matt$parasite)){
  for(pn in c("matt")){
    dt <- get(pn)
    pt <- ggplot(dt |> filter(parasite==pp), aes(x=-cost, y=precision, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn))
    print(pt)
  }
}
dev.off()

pdf("all_proportion.pdf")
for(pp in unique(matt$parasite)){
  for(pn in c("matt","bruno")){
    dt <- get(pn)
    pt <- ggplot(dt |> filter(parasite==pp), aes(x=n_individ, y=proportion, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn))
    print(pt)
  }
}
dev.off()

pdf("all_tradeoff.pdf")
for(pp in unique(matt$parasite)){
  for(pn in c("matt","bruno")){
    dt <- get(pn)
    pt <- ggplot(dt |> filter(parasite==pp), aes(x=-cost, y=proportion, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn))
    print(pt)
  }
}
dev.off()

pdf("all_n_v_cost.pdf")
for(pp in unique(matt$parasite)){
  mc <- 0
  for(pn in c("matt","bruno")){
    dt <- get(pn)
    mc <- max(mc, max(dt$cost))
  }

  for(pn in c("matt","bruno")){
    dt <- get(pn) |> filter(design!="SS12")

    pt <- ggplot(dt |> filter(parasite==pp), aes(x=n_individ, y=cost, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn)) +
      ylim(0, mc)
    print(pt)
  }
}
dev.off()

pdf("all_tradeoff.pdf")
for(pp in unique(matt$parasite)){
  #  for(pn in c("matt","bruno")){
  {
    pn <- "matt"
    dt <- get(pn)
    pt <- ggplot(dt |> filter(parasite==pp), aes(x=-cost, y=proportion, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn))
    print(pt)
  }
}
dev.off()

pdf("all_missingness.pdf")
for(pp in unique(matt$parasite)){
#  for(pn in c("matt","bruno")){
  {
  pn <- "matt"
    dt <- get(pn)
    pt <- ggplot(dt |> filter(parasite==pp), aes(x=n_individ, y=miss, col=design, group=design)) +
      geom_line() +
      #  geom_point() +
      facet_grid(mean_epg ~ method) +
      ggtitle(str_c(pp, " - ", pn))
    print(pt)
  }
}
dev.off()

stop()

# Back of the envelope calculation of cost for biggest discrepancy for NS11 and NS12:
both |> filter(parasite=="trichuris", method=="FP", design%in%c("NS11","NS12"), n_individ==1000, mean_epg==49.7) |> select(design, matt, bruno)

# Conclusion:  I over-estimate costs for NS12, Bruno under-estimates costs for NS11
pars <- survey_parameters(c("NS_12"), "trichuris","fecpak")
scen <- survey_scenario("trichuris")[3,]
nind <- 1000
matt <- survey_sim(n_individ=nind, scenario=scen, parameters=pars, output="extended")

# Cost of consumables:
conscost <- nind*(pars$cost_sample + pars$cost_aliquot_pre) +
  nind*(pars$cost_sample + pars$cost_aliquot_post)

tcount <- matt |> summarise(time = mean(time_count)) |> pull(time)
tcount
ttotal <- nind*(pars$time_demography + pars$time_prep_pre + pars$time_record) +
  nind*(pars$time_demography + pars$time_prep_post + pars$time_record*pars$n_aliquot_post) + tcount
ndays <- ttotal / (3*4*60*60)

# Cost of salary = transport
salcost <- ndays*4*pars$cost_salary

conscost + 2*salcost
matt |> summarise(mean(total_cost))

matt |> summarise(mean(consumables_cost))
conscost
matt |> summarise(mean(salary_cost))
salcost
matt |> summarise(mean(travel_cost))
salcost

# Compare timings for specialisation:
library("eggSim")

n_individ_us <- c(100,200,500,1000)
scen <- survey_scenario()

params1 <- survey_parameters(c("NS_12"))
params2 <- lapply(params1, function(x){
  x$design <- "NS"
  x
})

# Difference is maybe around 5% ??
system.time({
  results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params2, iterations=1e3)
})
system.time({
  results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params1, iterations=1e3)
})

# Test utility functions:
eggSim:::Rcpp_rbeta_cv(5, 0.5, 0)
eggSim:::Rcpp_rbeta_cv(5, 0.5, 0.1)
eggSim:::Rcpp_rgamma_cv(5, 5, 0)
eggSim:::Rcpp_rgamma_cv(5, 5, 1)
eggSim:::Rcpp_rnbinom_cv(50, 5, 0)
eggSim:::Rcpp_rnbinom_cv(50, 5, 5)
eggSim:::Rcpp_count_time(eggSim:::Rcpp_rnbinom_cv(50, 5, 5), 1, 1, 0, 1)
# Maybe I should export these?

results |>
  group_by(design, parasite, method, n_individ, mean_epg, reduction, individ_cv, day_cv, aliquot_cv, reduction_cv, cutoff, weight, recovery) |>
  summarise(.groups="drop") ->
  parsets

rgamma_cv <- function(n, mu, cv)
{
  pow <- function(x,y) x^y
  k = pow(cv, -2.0)
  rv = rgamma(n, k, scale=mu/k);
  rv
}
mean(rgamma_cv(100, 10, 0.1))

rnbinom_cv <- function(n, mu, cv)
{
  pow <- function(x,y) x^y
  if(cv <= 0.0)
  {
    return(rpois(n, mu));
  }
  k = pow(cv, -2.0);
  rv = rnbinom(n, k, mu=mu);

  rv
}
mean(rnbinom_cv(100, 10, 0))

rbeta_cv <- function(n, mu, cv)
{
  pow <- function(x,y) x^y
  sd = cv * mu;
  a = mu * ( (mu*(1.0-mu) / pow(sd,2.0)) - 1.0 );
  b = (1.0 - mu) * ( (mu*(1.0-mu) / pow(sd,2.0)) - 1.0);
  stopifnot(a>0, b>0)
  rv = rbeta(n, a, b);
  rv
}
mean(rbeta_cv(100, 0.5, 0.2))

both <- full_join(matt |> rename(matt = "proportion") |> select(-code, -cost, -miss),
                  bruno |> rename(bruno = "power") |> select(-code, -cost, -miss)
) |>
  filter(n_individ %in% n_individ_us) |>
  select(-scenario)

forluc <- parsets |>
  filter(design=="SS_12") |>
  mutate(method = factor(method, levels=c("kk", "fecpak", "miniflotac"), labels=c("KK","FP","MF")) |> as.character()) |>
  full_join(
    both |>
      filter(design=="SS12") |>
      select(parasite, method, n_individ, mean_epg, reduction, matt, bruno)
    ) |>
  filter(!is.na(reduction)) |>
  arrange(desc(abs(matt-bruno)))

#save(forluc, rgamma_cv, rnbinom_cv, rbeta_cv, file="forluc.RData")
#write_excel_csv(forluc, "parameters.csv")

# Problems:  KK, SS11, SS12,

both$ratio <- with(both, matt/bruno)
lm(ratio ~ 0 + mean_epg + design + parasite + method, data=both)

both <- bind_rows(matt, bruno)
ggplot(both, aes(x=n_individ, y=proportion, col=code)) +
  geom_point() +
  facet_wrap(~ design + parasite + method + mean_epg + reduction)
