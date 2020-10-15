library(eggSim)

R <- 100
# Reduction (either geometric or arithmetic mean, depending on LP vs GP) - can be vectorised:
reduction <- 0.05
# Community means (either geometric or arithmetic mean, depending on LP vs GP):
community_mean <- 240 #c(100, 100, 100)

cv_between <- 0.6 #c(0.6,0.8,1)
# CV within individuals (day-to-day variation):
cv_within <- 0.4
# CV within individual and day (slide-to-slide variation):
cv_slide <- 0.75
# Variation in efficacy between individuals:
cv_reduction <- 0.5

dd <- cgpDataSim(1, 10^5, reduction, community_mean, cv_between, cv_within, cv_slide, cv_reduction, overall_mean = mean(community_mean), edt=24)

## Ensure that the mathematically calculated k matches the empirical k:
combined_k(cv_between, cv_within, cv_slide, cv_reduction)
mean(dd$PreCount)^2 / (var(dd$PreCount) - mean(dd$PreCount))
# == mean(dd$PreSlide)^2 / (var(dd$PreSlide))
mean(dd$PostCount1a)^2 / (var(dd$PostCount1a) - mean(dd$PostCount1a))
# == mean(dd$PostCount1b)^2 / (var(dd$PostCount1b) - mean(dd$PostCount1b))
# == mean(dd$PostCount2)^2 / (var(dd$PostCount2) - mean(dd$PostCount2))
# == mean(dd$PostSlide1a)^2 / (var(dd$PostSlide1a))

# The lengths of these must match:
stopifnot(length(cv_between)==length(community_mean))

# The number of individuals used for each design type:
N = c(NS=600, SS=109, SSR1=100, SSR2=92, SSR3=95)


reduction <- seq(0,1,by=0.05)

t <- Sys.time()
results <- eggSim(reduction, N, community_mean, cv_between, cv_within, cv_slide, cv_reduction, count=FALSE, log_constant=if(count) 1 else 0, screen_threshold = 0, edt=24, R=R)
difftime(Sys.time(), t)

ggplot(results %>% filter(!is.na(Bias)), aes(x=Target, y=Bias, ymin=LCI, ymax=UCI, col=Design)) +
  geom_hline(yintercept=0) +
  geom_errorbar(position='dodge') +
  geom_line() +
  geom_point() +
  facet_wrap(~ Set, scales='free_y') +
  xlab('True efficacy (arithmetic or geometric mean as indicated)')
