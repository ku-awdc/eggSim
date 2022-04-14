## Comparison of code bases for Bruno's paper

library("tidyverse")
library("eggSim")

n_individ_us <- seq(100,1000,by=5)
n_individ_us <- c(100,200,500)
params <- survey_parameters()
scen <- survey_scenario()

iters <- 1e4
iters <- 1e3

cl <- NULL
cl <- parallel::makePSOCKcluster(5)
system.time({
  results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=cl, output="extended")
})
parallel::stopCluster(cl)

results

if(!grepl("paper_2022_Bruno", getwd())) setwd("paper_2022_Bruno")

theme_set(theme_light())
pdf("graph_failure_rate.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=(failure/total_n), col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_bias_1.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=(efficacy_mean/(1-reduction)), col=design, group=design)) +
    geom_hline(yintercept=1, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_bias_2.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=bias_mean, col=design, group=design)) +
    geom_hline(yintercept=1, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_precision.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=efficacy_precision, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_proportion.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_cost.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=cost_mean, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_tradeoff.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=-cost_mean, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_tradeoff_zoom.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=-cost_mean, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method, scales="free") +
    coord_cartesian(ylim=c(0.5, NA)) +
    labs(col=pp)
  print(pt)
}
dev.off()


pdf("graph_tradeoff_precision.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=-cost_mean, y=efficacy_precision, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method, scales="free") +
    coord_cartesian(ylim=c(0.5, NA)) +
    labs(col=pp)
  print(pt)
}
dev.off()


## Investigating bias
results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=NULL, output="extended")


results |>
  group_by(parasite, design, scenario, method, n_individ, mean_epg, reduction, weight, recovery, mu_pre) |>
  mutate(bias_efficacy = (1-reduction)-efficacy) |>
  mutate(bias_reduction = (reduction)-(1-efficacy)) |>
  summarise(pre_mean=mean(pre_mean, na.rm=TRUE), post_mean=mean(post_mean, na.rm=TRUE), pre_imean=mean(pre_imean, na.rm=TRUE), post_imean=mean(post_imean, na.rm=TRUE), efficacy = mean(efficacy, na.rm=TRUE), bias_efficacy = mean(bias_efficacy, na.rm=TRUE), bias_reduction = mean(bias_reduction, na.rm=TRUE), .groups="drop") ->
  sumstats

sumstats

pdf("graph_diff.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(sumstats |> filter(parasite==pp), aes(x=n_individ, y=efficacy-(1-post_mean/pre_mean), col=design, group=design)) +
    geom_hline(yintercept=0, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()


pdf("graph_bias_reduction.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(sumstats |> filter(parasite==pp), aes(x=n_individ, y=reduction-(post_mean/pre_mean), col=design, group=design)) +
    geom_hline(yintercept=0, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_bias_efficacy.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(sumstats |> filter(parasite==pp), aes(x=n_individ, y=(1-reduction)-(1-post_mean/pre_mean), col=design, group=design)) +
    geom_hline(yintercept=0, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()




# No bias in reduction estimates at individual mean level:
pdf("graph_bias_1.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(sumstats |> filter(parasite==pp), aes(x=n_individ, y=(reduction/(post_imean/pre_imean)), col=design, group=design)) +
    geom_hline(yintercept=1, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()

# There is no bias in efficacy at individual level:
pdf("graph_bias_1.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(sumstats |> filter(parasite==pp), aes(x=n_individ, y=(pre_imean*reduction)/post_imean, col=design, group=design)) +
    geom_hline(yintercept=1, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp)
  print(pt)
}
dev.off()
