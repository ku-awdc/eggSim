### This file recreates the analyses presented in the 2020 paper:
# Coffeng LE, Levecke B, Hattendorf J, Walker M, Denwood MJ (under review). Survey design to monitor drug efficacy for the control of soil-transmitted helminthiasis and schistosomiasis. Clin Infect Dis.

## Code last updated 2020-11-30
## Intended for use with eggSim version 0.9.0

library('tidyverse')
library('eggSim')

set.seed(2020-11-24)

theme_set(theme_light())

xlabs <- c(
  NS = '"NS"["1x1"]',
  NS2 = '"NS"["2x1"]',
  NS3 = '"NS"["1x2"]',
  SS = '"SS  "',
  SSR1 = '"SSR"["1x1"]',
  SSR2 = '"SSR"["2x1"]',
  SSR3 = '"SSR"["1x2"]'
)

## Set a working directory subfolder within the package:
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
levels(output$Prevalence) <- unique(str_c(output$Prevalence, "% Apparent Baseline Prevalence (", as.character(output$OverallMean), ")"))

# Check that all simulations got a non-zero pre-tx FEC:
stopifnot(all(output$Communities > 0))

# For our interest:
ggplot(output, aes(x=Design, y=ScreenProp)) +
  geom_boxplot() +
  facet_grid(Efficacy ~ OverallMean)

output$DesignTxt <- factor(output$Design, levels=names(xlabs), labels=xlabs)

# Figure 1:
ggplot(output, aes(x=DesignTxt, y=ArithmeticEfficacy)) +
  geom_hline(aes(yintercept=yint), output %>% distinct(Efficacy, Prevalence, yint=100*(1-Reduction)), lty='dashed', col="black") +
  geom_boxplot(outlier.shape=1) +
  facet_grid(Efficacy ~ Prevalence, scales='fixed') +
  ylab("Observed Efficacy (%)") +
  xlab(NULL) +
  theme(panel.spacing = unit(1.5, "lines"), axis.text.x = element_text(angle=-45, vjust=0.5, hjust=0.25)) +
  scale_x_discrete(labels=function(x) parse(text=x))
ggsave(file.path(basewd, "figure_1.pdf"), height=6, width=7)

# For our interest:
ggplot(output %>% filter(Communities > 0), aes(x=ArithmeticEfficacy, col=Design)) +
  stat_ecdf() +
  facet_grid(OverallMean ~ TrueArithmetic, scales='free_x') +
  labs(y = "Empirical CDF")

# Useful table:
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

plotdata$Design2 <- factor(plotdata$Design, levels=names(xlabs), labels=xlabs)

ggplot(plotdata, aes(x=Target, y=Estimate, col=Prevalence, group=Prevalence)) +
  geom_hline(aes(yintercept=yint), yintdat, lty='dashed', col="black") +
  #  stat_smooth(method='lm', formula = y ~ poly(x,3), se=FALSE) +
  #  geom_point(size=0.1) +
  geom_line() +
  facet_grid(Type ~ Design2, scales='free_y', labeller=labeller(Design2=label_parsed, Type=label_value)) +
  ylab(NULL) +
  xlab("\nTrue Efficacy (%)") +
  scale_x_continuous(limits=c(1,99), breaks=c(10,30,50,70,90)) +
  guides(col = guide_legend("Apparent Baseline Prevalence:")) +
  theme(panel.spacing = unit(0.75, "lines"), legend.pos='bottom') #, legend.title=element_blank())
ggsave(file.path(basewd, "figure_2.pdf"), width=10, height=6)


save(output, plotdata, file=file.path(basewd, "simres.Rdata"))


## Appendix:

# NS strategy using two pre-treatment samples (and a single post-treatment sample) is worse than default NS:

output <- eggSim(reduction, budget=budget, community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, true_prevalence=true_prevalence, R=R, summarise=TRUE, design=c("NS1","NS2","NS4"))

tab <- output %>%
  select(Design, ObsPrev, Target, MedianBias, Variance) %>%
  arrange(ObsPrev, Target, Design)
write_excel_csv(tab, file.path(basewd, "table_S1.csv"))

# Look at relative variance for different SS and NS designs vs slide cost:
reduction <- c(0.01,seq(0.05,0.95,by=0.05))
community_mean <- c(pc10=3.16,pc15=5.45,pc20=8.45,pc25=12,pc35=23.5)
second_slide_cost <- c(0.1, 0.621, 1)

fig3out <- eggSim(reduction, budget=budget, second_slide_cost=second_slide_cost, design=c("NS","NS3","SSR3"), community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, true_prevalence=true_prevalence, R=R, summarise=TRUE, parallelise=parallel::detectCores()/2)  # Limit to physical cores as the RAM requirements are quite high
save(fig3out, file=file.path(basewd, "simres_suppl_fig_S1.Rdata"))

fig3out$Prevalence <- factor(str_c(round(fig3out$ObsPrev*100), "%"))
fig3out$DesignTxt <- with(fig3out, case_when(
  Design %in% c("NS3","SSR3") ~ str_c(Design, " - ", SecondSlideCost, " unit", ifelse(SecondSlideCost==1, "", "s")),
  TRUE ~ Design
))
fig3out$CostTxt <- with(fig3out, case_when(
  Design %in% c("NS3","SSR3") ~ str_c(SecondSlideCost, " unit", ifelse(SecondSlideCost==1, "", "s")),
  TRUE ~ ""
))

bias <- fig3out %>%
  select(Design, DesignTxt, CostTxt, Prevalence, Target, Estimate=MedianBias) %>%
  mutate(Type = "Median Bias")
variance <- fig3out %>%
  select(Design, DesignTxt, CostTxt, Prevalence, Target, Estimate=Variance) %>%
  mutate(Estimate = sqrt(Estimate), Type = "Std. Dev.")
relvariance <- variance %>%
  filter(Design=="NS") %>%
  select(Prevalence, Target, Reference=Estimate) %>%
  full_join(variance, by = c("Prevalence", "Target")) %>%
  mutate(Estimate = (Estimate / Reference), Type = "Relative SD")

plot3data <- bind_rows(bias, variance, relvariance) %>%
  mutate(Type = factor(Type, levels=c("Median Bias", "Std. Dev.", "Relative SD")))

yintdat <- plot3data %>%
  distinct(Design, DesignTxt, Type) %>%
  mutate(yint = case_when(
    Type == "Relative Variance" ~ 1,
    Type == "Relative SD" ~ 1,
    TRUE ~ 0
  ))

plot3data$Design2 <- factor(plot3data$Design, levels=names(xlabs), labels=xlabs)

ggplot(plot3data %>% filter(Type == "Relative SD", Design%in%c("NS3","SSR3")), aes(x=Target, y=Estimate, col=Prevalence, group=Prevalence)) +
  geom_hline(aes(yintercept=yint), yintdat %>% filter(Type == "Relative SD"), lty='dashed', col="black") +
  #  stat_smooth(method='lm', formula = y ~ poly(x,3), se=FALSE) +
  #  geom_point(size=0.1) +
  geom_line() +
  facet_grid(Design2 ~ CostTxt, scales='fixed', labeller=labeller(Design2=label_parsed, CostTxt=label_value)) +
  ylab(NULL) +
  xlab("\nTrue Efficacy (%)") +
  ylab("Relative standard deviation\n") +
  scale_x_continuous(limits=c(1,99), breaks=c(10,30,50,70,90)) +
  guides(col = guide_legend("Apparent Baseline Prevalence:")) +
  theme(panel.spacing = unit(0.75, "lines"), legend.pos='bottom') #, legend.title=element_blank())
ggsave(file.path(basewd, "figure_S1.pdf"), width=6, height=4)

