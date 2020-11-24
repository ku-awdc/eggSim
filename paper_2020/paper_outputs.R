library('tidyverse')
library('eggSim')

set.seed(2020-11-24)

theme_set(theme_light())

### Generate figures and tables for 2020 paper
basewd <- "paper_2020"

## Fixed parameters:
R <- 5e3
budget <- 1200
cv_between <- 1.5
cv_within <- 0.75
cv_slide <- 0.25
cv_reduction <- 0
true_prevalence <- 1
# Overall k is about 0.23
combined_k(cv_between, cv_within, cv_slide, cv_reduction)

## Variable parameters:
community_mean <- c(3.16, 23.5)  # In EPG assuming 1/24 grams
reduction <- c(0.05, 0.2)

## Get output:
output <- eggSim(reduction, budget=budget, community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, true_prevalence=true_prevalence, R=R, summarise=FALSE)

# Prettify:
output$Efficacy <- factor(str_c(output$TrueArithmetic, "% True Efficacy"), levels=str_c(rev(100*(1-reduction)), "% True Efficacy"))
output$OverallMean <- factor(str_c(output$OverallMean, " EPG"), levels=str_c(community_mean, " EPG"))
output$Prevalence <- factor(round(output$ObsPrev*100))
levels(output$Prevalence) <- unique(str_c(output$Prevalence, "% Apparent Prevalence (", as.character(output$OverallMean), ")"))

# Check that all simulations got a non-zero pre-tx FEC:
stopifnot(all(output$Communities > 0))

# For our interest:
ggplot(output, aes(x=Design, y=ScreenProp)) +
  geom_boxplot() +
  facet_grid(Efficacy ~ OverallMean)

# Figure 1:
ggplot(output, aes(x=Design, y=ArithmeticEfficacy)) +
  geom_hline(aes(yintercept=yint), output %>% distinct(Efficacy, Prevalence, yint=100*(1-Reduction)), lty='solid', col="grey50") +
  geom_boxplot() +
  facet_grid(Efficacy ~ Prevalence, scales='fixed') +
  ylab("Observed Efficacy (%)") +
  xlab(NULL) +
  theme(panel.spacing = unit(1.5, "lines"))
ggsave(file.path(basewd, "figure_1.pdf"), height=6, width=7)

# For our interest:
ggplot(output %>% filter(Communities > 0), aes(x=ArithmeticEfficacy, col=Design)) +
  stat_ecdf() +
  facet_grid(OverallMean ~ TrueArithmetic, scales='free_x')

# Useful table?
tab <- output %>%
  mutate(Success = 100 * sum(Communities > 0) / n()) %>%
  filter(Communities > 0) %>%
  mutate(ObsPrev = ObsPrev*100) %>%
  group_by(OverallMean, ObsPrev, Prevalence, Efficacy, Design) %>%
  summarise(MeanN = mean(N), TotalBudget = mean(ScreenBudget + SampleBudget), ScreenProp = mean(ScreenProp), Mean = mean(ArithmeticEfficacy), Median = median(ArithmeticEfficacy), SD = sd(ArithmeticEfficacy), .groups='drop')
write_excel_csv(tab, file.path(basewd, "table_2.csv"))


## Calibrate mu to get desired prevalence levels:
community_mean <- c(pc10=3.16,pc15=5.45,pc20=8.45,pc25=12,pc35=23.5)
sumout <- eggSim(0.05, budget=budget, community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, true_prevalence=true_prevalence, R=R, summarise=TRUE, parallelise=length(community_mean))
sumout %>%
  filter(Design=="NS") %>%
  select(OverallMean, ObsPrev, Success)


## Figure 2:
reduction <- c(0.01,seq(0.05,0.95,by=0.05))
fig2out <- eggSim(reduction, budget=budget, community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, true_prevalence=true_prevalence, R=R, summarise=TRUE, parallelise=parallel::detectCores()/2)  # Limit to physical cores as the RAM requirements are quite high

fig2out$Prevalence <- factor(str_c(round(fig2out$ObsPrev*100), "%"))

bias <- fig2out %>%
  select(Design, Prevalence, Target, Estimate=MedianBias) %>%
  mutate(Type = "Median Bias")
variance <- fig2out %>%
  select(Design, Prevalence, Target, Estimate=Variance) %>%
  mutate(Estimate = sqrt(Estimate), Type = "Std. Dev.")
relvariance <- variance %>%
  filter(Design=="NS") %>%
  select(Prevalence, Target, Reference=Estimate) %>%
  full_join(variance, by = c("Prevalence", "Target")) %>%
  mutate(Estimate = (Estimate / Reference), Type = "Relative SD")

plotdata <- bind_rows(bias, variance, relvariance) %>%
  mutate(Type = factor(Type, levels=c("Median Bias", "Std. Dev.", "Relative SD")))

yintdat <- plotdata %>%
  distinct(Design, Type) %>%
  mutate(yint = case_when(
    Type == "Relative Variance" ~ 1,
    Type == "Relative SD" ~ 1,
    TRUE ~ 0
  ))
ggplot(plotdata, aes(x=Target, y=Estimate, col=Prevalence, group=Prevalence)) +
  geom_hline(aes(yintercept=yint), yintdat, lty='solid', col="grey50") +
  #  stat_smooth(method='lm', formula = y ~ poly(x,3), se=FALSE) +
  #  geom_point(size=0.1) +
  geom_line() +
  facet_grid(Type ~ Design, scales='free_y') +
  ylab(NULL) +
  xlab("\nTrue Efficacy (%)") +
  scale_x_continuous(limits=c(1,99), breaks=c(10,30,50,70,90)) +
  guides(col = guide_legend("Apparent Prevalence:")) +
  theme(panel.spacing = unit(0.75, "lines"), legend.pos='bottom') #, legend.title=element_blank())
ggsave(file.path(basewd, "figure_2.pdf"), width=10, height=6)


save(output, plotdata, file=file.path(basewd, "simres.Rdata"))
