library(brms)
library(tidyverse)
library(posterior)

# model parameters
n_iter <- 50000
n_thin <- 25

# relatedness matrix:
snp.matrix <- read.table('Samples_snpMatrix.txt', h = T, row.names = 1)

# Whole dataset
data <- read.table('PSdataset.txt', h = T)
data$animal <- data$SampleID

# NA Age first repr for PP
data_pp <- read.table('PPdataset.txt', h = T)
data_pp$animal <- data_pp$SampleID

# models
########
PP_mod <- brms::brm(Pup_produced ~ Birth_year_SAM_ANN_scaled + Repr_year_SAM_ANN_scaled + Mismatch_ANN_scaled_abs + Age_new_scaled + I(Age_new_scaled^2) + Froh_scaled + JcatAll + Froh_scaled*Repr_year_SAM_ANN_scaled + Age_new_scaled*Repr_year_SAM_ANN_scaled + I(Age_new_scaled^2)*Repr_year_SAM_ANN_scaled + Females_count_scaled + Froh_scaled*JcatAll + Froh_scaled*Age_new_scaled + I(Age_new_scaled^2)*Froh_scaled + (1|SampleID) + (1|gr(animal, cov = snp.matrix)), 
				  data = data_pp, 
				  data2 = list(snp.matrix = snp.matrix),
				  family = bernoulli, 
				  chains = 4, iter = n_iter, thin = n_thin, 
				  set_prior("normal(0,1)", class = "b"),
				  cores = 4,
				  control = list(adapt_delta = 0.99, max_treedepth = 15))

plot(PP_mod)

#####################

PS_mod <- brms::brm(Pup_survived ~ Birth_year_SAM_ANN_scaled + Pup_sex_bin + Repr_year_SAM_ANN_scaled + Mismatch_ANN_scaled_abs + Age_scaled + I(Age_scaled^2) + Froh_scaled + JcatAll + Froh_scaled*Repr_year_SAM_ANN_scaled + Age_scaled*Repr_year_SAM_ANN_scaled + I(Age_scaled^2)*Repr_year_SAM_ANN_scaled + Females_count_scaled + Froh_scaled*JcatAll + Froh_scaled*Age_scaled + I(Age_scaled^2)*Froh_scaled + (1|SampleID) + (1|gr(animal, cov = snp.matrix)), 
				  data = data, 
				  data2 = list(snp.matrix = snp.matrix),
				  family = bernoulli, 
				  chains = 4, iter = n_iter, thin = n_thin, 
				  set_prior("normal(0,1)", class = "b"),
				  cores = 4,
				  control = list(adapt_delta = 0.99, max_treedepth = 15))

plot(PS_mod)

##################################

# make pop grids
newdata_pp <- expand.grid(
  Repr_year_SAM_ANN_scaled = seq(-2, 2, length.out = 50), # increase max if you want to include harsher environment
  Froh_scaled = c(-2, 2), # low, high inbreeding
  
  # hold others constant (adjust as needed)
  Birth_year_SAM_ANN_scaled = 0,
  Mismatch_ANN_scaled_abs = 0,
  Age_new_scaled = 0,
  JcatAll = 1,
  Females_count_scaled = 0
)

newdata_ps <- expand.grid(
  Repr_year_SAM_ANN_scaled = seq(-2, 2, length.out = 50), # increase max if you want to include harsher environment
  Froh_scaled = c(-2, 2), # low, high inbreeding
  
  # hold others constant (adjust as needed)
  Birth_year_SAM_ANN_scaled = 0,
  Mismatch_ANN_scaled_abs = 0,
  Age_scaled = 0,
  JcatAll = 1,
  Females_count_scaled = 0,
  Pup_sex_bin = 0
)

# posterior draws
pp_draws <- posterior_epred(PP_mod, newdata = newdata_pp, re_formula = NA)
ps_draws <- posterior_epred(PS_mod, newdata = newdata_ps, re_formula = NA)

# Compute early fitness (PP * PS)
w_draws <- pp_draws * ps_draws

# summarize
w_summary <- summarise_draws(
  as_draws_matrix(w_draws),
  mean = mean,
  q2.5 = ~quantile(.x, 0.025),
  q97.5 = ~quantile(.x, 0.975)
)

results <- cbind(newdata_ps, w_summary)
# write.table(results, 'results_fig4.txt', quote = F, col.names = T, row.names = F, sep = '\t')

# Define function to compute beta for high F vs low F and calculate it for a given SAM value:
compute_beta_env <- function(draws, newdata, trait_name, f_low, f_high, env_value) {
  
  idx_low <- which(newdata[[trait_name]] == f_low &
                     newdata$Repr_year_SAM_ANN_scaled == env_value)
  
  idx_high <- which(newdata[[trait_name]] == f_high &
                      newdata$Repr_year_SAM_ANN_scaled == env_value)
  
  delta_f <- f_high - f_low
  
  beta <- (draws[, idx_high] - draws[, idx_low]) / delta_f
  
  return(beta)
}

# BAD years (e.g. SAM = 2)
beta_pp_bad <- compute_beta_env(pp_draws, newdata_pp, "Froh_scaled", -2, 2, env_value = 2)
beta_ps_bad <- compute_beta_env(ps_draws, newdata_ps, "Froh_scaled", -2, 2, env_value = 2)
beta_w_bad  <- compute_beta_env(w_draws,  newdata_ps, "Froh_scaled", -2, 2, env_value = 2)

# GOOD years (e.g. SAM = 2)
beta_pp_good <- compute_beta_env(pp_draws, newdata_pp, "Froh_scaled", -2, 2, env_value = -2)
beta_ps_good <- compute_beta_env(ps_draws, newdata_ps, "Froh_scaled", -2, 2, env_value = -2)
beta_w_good  <- compute_beta_env(w_draws,  newdata_ps, "Froh_scaled", -2, 2, env_value = -2)

# Summarize
# summary function
summarize_beta <- function(beta) {
  data.frame(
    mean = mean(beta),
    sd = sd(beta),
    q2.5 = quantile(beta, 0.025),
    q97.5 = quantile(beta, 0.975)
  )
}

summary_bad <- rbind(
  cbind(stage = "PP", summarize_beta(beta_pp_bad)),
  cbind(stage = "PS", summarize_beta(beta_ps_bad)),
  cbind(stage = "W", summarize_beta(beta_w_bad))
)

summary_bad$env <- "High SAM"

summary_good <- rbind(
  cbind(stage = "PP", summarize_beta(beta_pp_good)),
  cbind(stage = "PS", summarize_beta(beta_ps_good)),
  cbind(stage = "W", summarize_beta(beta_w_good))
)

summary_good$env <- "Low SAM"

selection_env_df <- rbind(summary_bad, summary_good)

# reformat data
df_fig4 <- bind_rows(
  tibble(value = beta_pp_good, type = "PP", env = "Low SAM"),
  tibble(value = beta_ps_good, type = "PS", env = "Low SAM"),
  tibble(value = beta_w_good,  type = "PW",  env = "Low SAM"),
  
  tibble(value = beta_pp_bad, type = "PP", env = "High SAM"),
  tibble(value = beta_ps_bad, type = "PS", env = "High SAM"),
  tibble(value = beta_w_bad,  type = "PW",  env = "High SAM")
)

df_fig4$env <- factor(df_fig4$env, levels = c('Low SAM', 'High SAM'))

# fix type for bad values
df_fig4$type[(length(beta_pp_good) + length(beta_ps_good) + length(beta_w_good) + 1):
          nrow(df_fig4)] <- rep(c("PP","PS","PW"),
                            times = c(length(beta_pp_bad),
                                      length(beta_ps_bad),
                                      length(beta_w_bad)))

# write.table(df_fig4, 'df_fig4.txt', quote = T, col.names = T, row.names = F, sep = '\t')
