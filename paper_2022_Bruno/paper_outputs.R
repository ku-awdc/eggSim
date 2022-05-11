## Comparison of code bases for Bruno's paper

library("tidyverse")
library("eggSim")

# Get default parameter values for all scenarios:
n_individ_us <- seq(100,5000,by=5)
params <- survey_parameters()
scen <- survey_scenario()

# Set iterations (1e3 for testing, 1e4 for final results):
iters <- 1e4
#iters <- 1e3

# Get all results using 6 parallel cores:
results <- survey_sim(n_individ = n_individ_us, scenario=scen, parameters = params, iterations=iters, cl=6L, output="summarised")

# Remove estimates for means and efficacies where these are based on fewer than 50 successful iterations:
results <- results |>
  mutate(across(c(pre_mean, post_mean, efficacy_mean, efficacy_variance), ~ case_when((total_n-failure_n) >= 50L ~.x, TRUE ~ NA_real_)))

# Maximum N to show on graphs:
nmax <- 1000

# Max cost is defined by NS11 with sample size of 1000:
results |>
  filter(design=="NS_11", n_individ==nmax) |>
  group_by(parasite) |>
  summarise(maxcost = max(cost_mean), .groups="drop")
# i.e. around 10000
maxcost <- 10000

# Verify that all designs get over the max cost for all methods/scenarios/parasites
results |>
  group_by(design,method,scenario,parasite) |>
  summarise(maxcost_dmsp = max(cost_mean), .groups="drop") |>
  arrange(maxcost_dmsp) ->
  allmcost
allmcost
with(allmcost, stopifnot(all(maxcost_dmsp >= maxcost)))

# Change working directory if needed:
if(!grepl("paper_2022_Bruno", getwd())) setwd("paper_2022_Bruno")

# Make graphs:
theme_set(theme_light())

# Fig 1: Pipelines of analysis for one STH  (HK) and test (KK)  mean = 3.7; (failure rate; number of individuals vs. power; number of individuals vs. cost; cost vs. power). All the others in supplementary figures. Also use this example for the sensitivity analysis.

fig1data <- results |>
  filter(parasite=="hookworm", method=="kk", scenario==1L) |>
  mutate(nclass = case_when(n_individ > nmax ~ "over1k", TRUE ~ "under1k"))

fig1a <- ggplot(fig1data |> filter(n_individ <= nmax), aes(x=n_individ, y=failure_n/total_n, col=design)) +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  ylab("Failure Rate") +
  scale_x_continuous() +
  xlab("Total Individuals") +
  theme(legend.pos="none")
fig1a
fig1b <- ggplot(fig1data |> filter(n_individ <= nmax), aes(x=n_individ, y=below_cutoff/total_n, col=design)) +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  ylab("Proportion Below Threshold") +
  scale_x_continuous() +
  xlab("Total Individuals") +
  theme(legend.pos="none")
fig1b
fig1c <- ggplot(fig1data |> filter(n_individ <= nmax), aes(x=n_individ, y=cost_mean/1e3, col=design)) +
  geom_line() +
  scale_y_continuous() +
  ylab("Mean Cost ($k)") +
  scale_x_continuous() +
  xlab("Total Individuals") +
  theme(legend.pos="none")
fig1c
fig1d <- ggplot(fig1data |> filter(cost_mean <= maxcost), aes(x=cost_mean/1e3, y=below_cutoff/total_n, col=design, lty=nclass)) +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  ylab("Proportion Below Threshold") +
  scale_x_continuous(trans="reverse") +
  xlab("Mean Cost ($k)") +
  theme(legend.pos="none") +
  scale_linetype_manual(values=c(under1k="solid", over1k="dotted"))
fig1d

# From http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(fig1a + theme(legend.pos="bottom", legend.title=element_blank()))

fig1 <- grid.arrange(fig1a, fig1b, fig1c, fig1d, fig1e, ncol=2, nrow = 3,
             layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
             widths = c(2.7, 2.7), heights = c(2.5, 2.5, 0.5))
ggsave("Figure_1.pdf", fig1, height=7, width=6)

# Fig 2: KK best options to reduce the number (Hookworm: lowest value; 6 panels each representing one survey design and 3 methods (lines); cost. Vs power

fig2data <- results |>
  filter(parasite=="hookworm", scenario==1L) |>
  mutate(nclass = case_when(n_individ > nmax ~ "over1k", TRUE ~ "under1k")) |>
  mutate(design = factor(design, levels=c("SS_11","NS_11","SSR_11","SS_12","NS_12","SSR_12")))

fig2 <- ggplot(fig2data |> filter(cost_mean <= maxcost), aes(x=cost_mean/1e3, y=below_cutoff/total_n, col=method, lty=nclass)) +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  ylab("Proportion Below Threshold") +
  scale_x_continuous(trans="reverse") +
  xlab("Mean Cost ($k)") +
  facet_wrap(~design, scales="fixed") +
  theme(legend.pos="bottom") +
  scale_linetype_manual(values=c(under1k="solid", over1k="dotted")) +
  guides(linetype="none", col=guide_legend(""))
fig2
ggsave("Figure_2.pdf", fig2, height=6, width=7)

# Fig 3: Cost vs. power for KK only across all scenarios of survey design and endmicity. Figures complete range of values. Dashed lines from where it becomes not feasible

fig3data <- results |>
  filter(method=="kk") |>
  mutate(nclass = case_when(n_individ > nmax ~ "over1k", TRUE ~ "under1k")) |>
  mutate(scenario = str_c("Scenario ", scenario))

fig3 <- ggplot(fig3data |> filter(cost_mean <= maxcost), aes(x=cost_mean/1e3, y=below_cutoff/total_n, col=design, lty=nclass)) +
  geom_line() +
  scale_y_continuous(labels=scales::percent) +
  ylab("Proportion Below Threshold") +
  scale_x_continuous(trans="reverse") +
  xlab("Mean Cost ($k)") +
  facet_grid(scenario ~ parasite) +
  theme(legend.pos="bottom") +
  scale_linetype_manual(values=c(under1k="solid", over1k="dotted")) +
  guides(linetype="none", col=guide_legend(""))
fig3
ggsave("Figure_3.pdf", fig3, height=7, width=6)

## TODO: work out how to shade the area under the pareto front??



stop("Older graphs below here")

ggplot(results |> filter(parasite=="trichuris"), aes(x=-cost_mean, y=proportion_below, col=design, group=design)) +
  geom_line() +
  #  geom_point() +
  facet_grid(mean_epg ~ method) +
  labs(col=pp) +
  xlim(-maxcost, 0)

pdf("figure2.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=-cost_mean, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(-maxcost, 0)
  print(pt)
}
dev.off()


pdf("graph_failure_rate.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=(failure_n/total_n), col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()

pdf("graph_bias_1.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=(true_efficacy - efficacy_mean), col=design, group=design)) +
    geom_hline(yintercept=0, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()

pdf("graph_bias_2.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=(true_efficacy - (1 - post_mean/pre_mean)), col=design, group=design)) +
    geom_hline(yintercept=0, lty='dashed') +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()

pdf("graph_precision.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=1 / efficacy_variance, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()

pdf("graph_proportion.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()

pdf("graph_cost.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=cost_mean, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()


pdf("graph_cost_sd.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=n_individ, y=cost_variance^0.5, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(0, nmax)
  print(pt)
}
dev.off()

pdf("graph_tradeoff.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=-cost_mean, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method) +
    labs(col=pp) +
    xlim(-maxcost, 0)
  print(pt)
}
dev.off()

pdf("graph_tradeoff_zoom.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp), aes(x=-cost_mean, y=proportion_below, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method, scales="free") +
    coord_cartesian(ylim=c(0.5, NA), xlim=c(-maxcost, 0)) +
    labs(col=pp)
  print(pt)
}
dev.off()

pdf("graph_cost_vs_precision.pdf")
for(pp in unique(scen$parasite)){
  pt <- ggplot(results |> filter(parasite==pp, cost_mean <= maxcost), aes(x=cost_mean, y=1/efficacy_variance, col=design, group=design)) +
    geom_line() +
    #  geom_point() +
    facet_grid(mean_epg ~ method, scales="free_y") +
    labs(col=pp)
  print(pt)
}
dev.off()

stop()

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
