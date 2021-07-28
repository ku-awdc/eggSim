# Test Stan model for variance decomposition of egg counts
# Author: Luc Coffeng
# Date created: 23 July 2021

# Set up R Session ----
rm(list = ls())

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(data.table)


# Set paths ----
base_dir <- "/Users/luc/Documents/Research/Ghent/06_Quantifying_egg_count_CV"
data_dir <- file.path(base_dir, "01_Data")
code_dir <- file.path(base_dir, "02_Code")
output_dir <- file.path(base_dir, "03_Output")
NB_model_file <- file.path(code_dir, "00_Source/KK_variance_structure_gamma_ZI.stan")
lognormal_model_file <- file.path(code_dir, "00_Source/KK_variance_structure_lognormal_ZI.stan")

# Simulation setup ----
n_clus <- 2            # Number of clusters
n_clus <- 1            # Number of clusters
n_ind <- 2e2           # Number of individuals per cluster
n_days <- 2            # Number of days per individual
n_slides <- 2          # Number of slides per day

mu_clus <- c(4, 10)       # Mean egg count per cluster (must be of length 'n_clus')
mu_clus <- 10             # Mean egg count per cluster (must be of length 'n_clus')
                          # N.B.: used as arithm. mean in NB model and used as geom mean in lnorm model
p_zero_clus <- c(.1, .3)  # Proportion not infected per cluster (must be of length 'n_clus')
p_zero_clus <- .0         # Proportion not infected per cluster (must be of length 'n_clus')
cv_between <- sqrt(4)     # CV for variation between individuals
cv_within <- sqrt(2)      # CV for variation between days within individuals
cv_slide <- sqrt(1)       # CV for variation between slides within individuals and days
cv_total <- sqrt((cv_between^2 + 1) * (cv_within^2 + 1) * (cv_slide^2 + 1) - 1)

k_between <- cv_between^-2
k_within <- cv_within^-2
k_slide <- cv_slide^-2

sigma_between <- sqrt(log(cv_between^2 + 1))
sigma_within <- sqrt(log(cv_within^2 + 1))
sigma_slide <- sqrt(log(cv_slide^2 + 1))

data_faux <- data.table(expand.grid(clus = 1:n_clus,
                                    ID = 1:n_ind,
                                    day = 1:n_days,
                                    slide = 1:n_slides))
data_faux[, mu := mu_clus[clus]]
data_faux[, p_zero := p_zero_clus[clus]]
data_faux[, ID := as.integer(interaction(ID, clus))]
setkey(data_faux, clus, ID, day, slide)


# Prior distribution ----
prior_nb <- list(
  # Log of arithmetric mean egg count (N.B. not the geometric mean)
  mu_log_mu = log(mean(mu_clus)) - 1^2 / 2, mu_log_sd = 1,
  
  # Log of total CV
  cv_log_mu = log(cv_total) - 1^2 / 2, cv_log_sd = 1,
  
  # Variance decomposition
  alpha = rep(1, 3),  # uniform prior
  
  # Zero-inflation probability
  p_zero_a = 1, p_zero_b = 1)  # uniform prior
              
prior_lnorm <- prior_nb
prior_lnorm$mu_log_mu <- mean(log(mu_clus))
prior_lnorm$mu_log_sd <- max(1, sd(log(mu_clus)) * 2, na.rm = TRUE)


# Generate and prep data for the gamma model ----
data_NB <- copy(data_faux)

set.seed(12345)
data_NB[, mu_ind := rgamma(n = 1, shape = k_between, rate = k_between / mu[1]), by = .(clus, ID)]
data_NB[, mu_day := rgamma(n = 1, shape = k_within, rate = k_within / mu_ind[1]), by = .(clus, ID, day)]
data_NB[, mu_slide := rgamma(n = .N, shape = k_slide, rate = k_slide / mu_day)]
data_NB[, count := rpois(n = .N, mu_slide)]

if (any(p_zero_clus > 0)) {
  ID_index <- data_NB[, sample(x = unique(ID), size = round(p_zero[1] * n_ind)),
                   by = clus][, V1]
  data_NB[ID %in% ID_index, count := 0]
  data_NB[, all_zero := sum(count) == 0, by = ID]
  data_NB[, ID := match(ID, unique(ID)),
          by = all_zero]
  setkey(data_NB, all_zero, clus, ID, day, slide)
}

if (all(p_zero_clus == 0)) {
  data_NB_prepped <- list(N_clus = n_clus,
                          N_ind = data_NB[, length(unique(ID))],
                          N_days = n_days,
                          N_obs = data_NB[, .N],
                          
                          count = data_NB[, count],
                          clus = data_NB[, clus],
                          ID = data_NB[, ID],
                          day = data_NB[, day],
                          
                          N_ind0 = 0,
                          clus0 = integer(),
                          N_count0 = array(0L, dim = c(0, n_days)))
} else {
  data_NB_prepped <- list(N_clus = n_clus,
                          N_ind = data_NB[!(all_zero), length(unique(ID))],
                          N_days = n_days,
                          N_obs = data_NB[!(all_zero), .N],
                          
                          count = data_NB[!(all_zero), count],
                          clus = data_NB[!(all_zero), clus],
                          ID = data_NB[!(all_zero), ID],
                          day = data_NB[!(all_zero), day],
                          
                          N_ind0 = data_NB[(all_zero), length(unique(ID))],
                          clus0 = data_NB[(all_zero), clus[1], by = ID][, V1],
                          N_count0 = as.matrix(dcast(data_NB[(all_zero),
                                                             .N, by = .(ID, day)],
                                                     formula = ID ~ day,
                                                     value.var = "N"))[, -1])
}

# Test gamma model ----
fit_NB <- stan(file = NB_model_file,
               model_name = "NB_count_model",
               data = c(data_NB_prepped, prior_nb),
               pars = c("mu", "p_zero", "cv", "kappa",
                        "cv_between", "cv_within", "cv_slide",
                        "k_between", "k_within", "k_slide"),
               control = list(adapt_delta = .95,
                              max_treedepth = 12),
               chains = 4,
               iter = 2000)

print(fit_NB,
      probs = c(0.025, 0.975),
      digits_summary = 3)

pairs(fit_NB, pars = "kappa")
pairs(fit_NB, pars = c("cv_between", "cv_within", "cv_slide"))
pairs(fit_NB, pars = c("k_between", "k_within", "k_slide"))

# Generate data for the log-normal model ----
data_lnorm <- copy(data_faux)

set.seed(12345)
data_lnorm[, mu_geom_log_ind := rnorm(n = 1, mean = log(mu[1]), sd = sigma_between), by = .(clus, ID)]
data_lnorm[, mu_geom_log_day := rnorm(n = 1, mean = mu_geom_log_ind[1], sd = sigma_within), by = .(clus, ID, day)]
data_lnorm[, mu_geom_log_slide := rnorm(n = .N, mean = mu_geom_log_day, sd = sigma_slide)]
data_lnorm[, count := rpois(n = .N, exp(mu_geom_log_slide))]

if (any(p_zero_clus > 0)) {
  ID_index <- data_lnorm[, sample(x = unique(ID), size = round(p_zero[1] * n_ind)),
                         by = clus][, V1]
  data_lnorm[ID %in% ID_index, count := 0]
  data_lnorm[, all_zero := sum(count) == 0, by = ID]
  data_lnorm[, ID := match(ID, unique(ID)),
             by = all_zero]
  setkey(data_lnorm, all_zero, clus, ID, day, slide)
}

if (all(p_zero_clus == 0)) {
  data_lnorm_prepped <- list(N_clus = n_clus,
                             N_ind = data_NB[, length(unique(ID))],
                             N_days = n_days,
                             N_obs = data_NB[, .N],
                             
                             count = data_NB[, count],
                             clus = data_NB[, clus],
                             ID = data_NB[, ID],
                             day = data_NB[, day],
                             
                             N_obs0 = 0,
                             N_ind0 = 0,
                             clus0 = integer(),
                             ID0 = integer(),
                             day0 = integer())
} else {
  data_lnorm_prepped <- list(N_clus = n_clus,
                          N_ind = data_NB[!(all_zero), length(unique(ID))],
                          N_days = n_days,
                          N_obs = data_NB[!(all_zero), .N],
                          
                          count = data_NB[!(all_zero), count],
                          clus = data_NB[!(all_zero), clus],
                          ID = data_NB[!(all_zero), ID],
                          day = data_NB[!(all_zero), day],
                          
                          N_obs0 = data_NB[(all_zero), .N],
                          N_ind0 = data_NB[(all_zero), length(unique(ID))],
                          clus0 = data_NB[(all_zero), clus[1], by = ID][, V1],
                          ID0 = data_NB[(all_zero), ID],
                          day0 = data_NB[(all_zero), day])
}


# Test log-normal model ----
fit_lnorm <- stan(file = lognormal_model_file,
                  model_name = "lnorm_pois_count_model",
                  data = c(data_lnorm_prepped, prior_lnorm),
                  pars = c("mu_arith", "mu_geom", "p_zero", "cv", "kappa",
                           "cv_between", "cv_within", "cv_slide",
                           "sigma_between", "sigma_within", "sigma_slide"),
                  control = list(adapt_delta = .95,
                                 max_treedepth = 13),
                  chains = 4,
                  iter = 2e2)

print(fit_lnorm,
      probs = c(0.025, 0.975),
      digits_summary = 3)

pairs(fit_lnorm, pars = "kappa")
pairs(fit_lnorm, pars = c("cv_between", "cv_within", "cv_slide"))
pairs(fit_lnorm, pars = c("sigma_between", "sigma_within", "sigma_slide"))

### END OF CODE ### ----
