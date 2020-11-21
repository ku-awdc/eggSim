library('tidyverse')
library('eggSim')

theme_set(theme_light())

### Generate figures and tables for 2020 paper

## Fixed parameters:
budget <- 1200
cv_between <- 1.5
cv_within <- 0.75
cv_slide <- 0.25
cv_reduction <- 0
# Overall k is about 0.23
combined_k(cv_between, cv_within, cv_slide, cv_reduction)

## Variable parameters:
community_mean <- c(6, 24)  # In EPG assuming 1/24 grams
reduction <- c(0.05, 0.2)

## Get output:
output <- eggSim(reduction, budget=budget, community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, R=10^3, summarise=FALSE)

# Prettify:
output$Efficacy <- factor(str_c(output$TrueArithmetic, "% True Efficacy"), levels=str_c(rev(100*(1-reduction)), "% True Efficacy"))
output$OverallMean <- factor(str_c(output$OverallMean, " EPG"), levels=str_c(community_mean, " EPG"))

# Check that all simulations got a non-zero pre-tx FEC:
stopifnot(all(output$Communities > 0))

# For our interest:
ggplot(output, aes(x=Design, y=ScreenProp)) +
  geom_boxplot() +
  facet_grid(Efficacy ~ OverallMean)

# Figure 1:
ggplot(output %>% filter(Communities > 0), aes(x=Design, y=ArithmeticEfficacy)) +
  geom_boxplot() +
  facet_grid(Efficacy ~ OverallMean, scales='free_y') +
  ylab("Observed Efficacy") +
  xlab(NULL)
ggsave("figure_1.pdf")

# For our interest:
ggplot(output %>% filter(Communities > 0), aes(x=ArithmeticEfficacy, col=Design)) +
  stat_ecdf() +
  facet_grid(OverallMean ~ TrueArithmetic, scales='free_x')

# Useful table?
tab <- output %>%
  mutate(Success = 100 * sum(Communities > 0) / n()) %>%
  filter(Communities > 0) %>%
  group_by(OverallMean, Efficacy, Design) %>%
  summarise(Prevalence = mean(Prevalence), MeanN = mean(N), ScreenProp = mean(ScreenProp), Mean = mean(ArithmeticEfficacy), Median = median(ArithmeticEfficacy), Variance = var(ArithmeticEfficacy), .groups='drop')
write_excel_csv(tab, "table.csv")


sumout <- eggSim(reduction, budget=budget, community_mean=community_mean, cv_between=cv_between, cv_within=cv_within, cv_slide=cv_slide, cv_reduction=cv_reduction, R=10^3, summarise=TRUE)
sumout
