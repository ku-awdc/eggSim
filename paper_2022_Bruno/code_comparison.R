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
n_individ_us <- c(100,200,500,1000)
#n_individ_us <- 100

params <- survey_parameters()
scen <- survey_scenario()

results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params)

results |>
  group_by(design, parasite, method, n_individ, scenario, mean_epg, reduction) |>
  summarise(nonmiss = sum(!is.na(efficacy)), proportion = sum(efficacy < cutoff, na.rm=TRUE) / nonmiss, miss = 1-(nonmiss/n()), cost = mean(total_cost), .groups="drop") |>
  mutate(design = str_replace(design, "_", ""), code = "matt") |>
  mutate(method = factor(method, levels=c("kk", "fecpak", "miniflotac"), labels=c("KK","FP","MF")) |> as.character()) ->
  matt


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

save(forluc, rgamma_cv, rnbinom_cv, rbeta_cv, file="forluc.RData")
write_excel_csv(forluc, "parameters.csv")

theme_set(theme_light())
ggplot(both, aes(x=bruno, y=matt, group = mean_epg, col=interaction(n_individ, mean_epg, drop=TRUE, sep=" / "))) +
  geom_line() +
  geom_abline(slope=1, intercept=0) +
  geom_point() +
  facet_wrap(~ design + parasite + method, ncol=9)
ggsave("comparison_proportion.pdf", width=15, height=15)

both <- inner_join(matt |> rename(matt = "cost") |> select(-code, -miss),
                   bruno |> rename(bruno = "cost") |> select(-code, -miss)
)

theme_set(theme_light())
ggplot(both, aes(x=bruno/1e3, y=matt/1e3, group = mean_epg, col=interaction(n_individ, mean_epg, drop=TRUE, sep=" / "))) +
  geom_line() +
  geom_abline(slope=1, intercept=0) +
  geom_point() +
  facet_wrap(~ design + parasite + method, ncol=9)
ggsave("comparison_cost.pdf", width=15, height=15)


# Problems:  KK, SS11, SS12,

both$ratio <- with(both, matt/bruno)
lm(ratio ~ 0 + mean_epg + design + parasite + method, data=both)

both <- bind_rows(matt, bruno)
ggplot(both, aes(x=n_individ, y=proportion, col=code)) +
  geom_point() +
  facet_wrap(~ design + parasite + method + mean_epg + reduction)
