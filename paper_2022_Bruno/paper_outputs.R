## Comparison of code bases for Bruno's paper

library("tidyverse")
library("eggSim")

n_individ_us <- seq(100,1000,by=5)
params <- survey_parameters()
scen <- survey_scenario()

iters <- 1e4
#iters <- 1e3

cl <- parallel::makePSOCKcluster(5)
system.time({
  results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=cl, output="summarised")
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

pdf("graph_bias.pdf")
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

