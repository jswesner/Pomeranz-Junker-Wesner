# Bayesian models

# fit candidate models for ISD exponent and total log10 community biomass. 
# compare models using bayesian stacking
# make plots for main text

library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)
library(ggdist)
library(scales)

ID_key <- readRDS("data/sample_ID_key.RDS")
# read in explanatory variables
abiotic <- read.csv("data/abiotic.csv")
# standardize variables
abiotic_s <- scale(abiotic[,2:8]) 
# remove scaled attributes (mean and SD for each column)
scaled_attr <- attributes(abiotic_s)
# make scaled values into data frame and add siteID
abiotic_s <- as.data.frame(abiotic_s)
abiotic_s$siteID <- abiotic$Site

# # read in b-exponent and biomass data and wrangle objects
# MLEbins <- readRDS("data/MLEbins.RDS")
# MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

biomass <- readRDS("data/mean_biomass.RDS")
biomass <- biomass[,c("siteID", "ID", "year",
                      "sampleEvent", "u_biomass")]

# join data sets
#d <- left_join(MLEbins, biomass)
d <- left_join(biomass, abiotic_s)
#d <- left_join(biomass, ID_key)

# log 10 mean dry weight estimate
d$log_mg <- log10(d$u_biomass)
d$scale_mg <- d$u_biomass/max(d$u_biomass)
d$biomass_g <- d$u_biomass/1000

# community biomass models ------------------------------------------------

# global model using standardized variables
# temperature + nutrients + canopy
mod1 <- brm(data = d,
            biomass_g ~ mat.c + map.mm + tdn + tdp + canopy +
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              # 95% of slope prior between -1 and 1
              # response is on log10 scale
              # log-link exponentiates
              c(prior(normal(0, 1),
                      class = "b"),
                # 95% of intercept prior between 
                # -3.5 and 0.5
                # i.e., exponent value at mean values of
                # all standardized variables 
                prior(normal(0,2), 
                      class = "Intercept"),
                prior(exponential(1),
                      class="sd")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 2000,
            cores = 4)


# temp + nutrients
mod2 <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 8)
# "Climate model" temp + precipitation
mod3 <- update(mod1, formula. = . ~ . -canopy -tdn -tdp,
               cores = 8)
# temp + canopy
mod4 <- update(mod1,
               formula. = . ~ . -map.mm -tdn - tdp,
               cores = 8)
# just temperature
mod5 <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy)
# resources - autochthonous resources = TDN +TDP
# Allochthonous resources = canopy
mod6 <- update(mod1,
               formula. = . ~ . -map.mm,
               cores = 8)

dir.create("models_jsw")
saveRDS(mod1, file = "models_jsw/mod1.rds")
saveRDS(mod2, file = "models_jsw/mod2.rds")
saveRDS(mod3, file = "models_jsw/mod3.rds")
saveRDS(mod4, file = "models_jsw/mod4.rds")
saveRDS(mod5, file = "models_jsw/mod5.rds")
saveRDS(mod6, file = "models_jsw/mod6.rds")

# leave-one-out cross validation with bayesian stacking weights
loo_1 <- loo(mod1,
             reloo = TRUE,
             cores = 6)
loo_2 <- loo(mod2,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_3 <- loo(mod3,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_4 <- loo(mod4,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_5 <- loo(mod5,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)
loo_6 <- loo(mod6,
             reloo = TRUE,
             seed = TRUE,
             cores = 6)


saveRDS(loo_1, file = "models_jsw/loo_1.rds")
saveRDS(loo_2, file = "models_jsw/loo_2.rds")
saveRDS(loo_3, file = "models_jsw/loo_3.rds")
saveRDS(loo_4, file = "models_jsw/loo_4.rds")
saveRDS(loo_5, file = "models_jsw/loo_5.rds")
saveRDS(loo_6, file = "models_jsw/loo_6.rds")


# note on warnings:
# reloo uses the future package internally, and the future package recently updated to throw a warning when random numbers are generated.
# It seems like the warnings can be ignored based on this post: https://discourse.mc-stan.org/t/psis-loo-warns-about-an-unreliable-value/19437/2

# Allegedly you should be able to set seed=TRUE for the future package, but can't figure out how to do it. 

# The results I got previously (before the update) were similar. for now, I'm using as-is

# summarize model coefficients for SI ####
coef_mods_list <- list(global = as.data.frame(fixef(mod1)),
                       nutrient = as.data.frame(fixef(mod2)),
                       Climate  = as.data.frame(fixef(mod3)),
                       canopy = as.data.frame(fixef(mod4)),
                       T_only = as.data.frame(fixef(mod5)),
                       resources = as.data.frame(fixef(mod6)))
coef_names_list <- list(row.names(fixef(mod1)),
                        row.names(fixef(mod2)),
                        row.names(fixef(mod3)),
                        row.names(fixef(mod4)),
                        row.names(fixef(mod5)),
                        row.names(fixef(mod6)))
coef_mods_table <- map2(coef_mods_list,
                        coef_names_list,
                        ~cbind(.x, coef_name = .y))

coef_mods_table <- bind_rows(coef_mods_table, .id = "MOD")
row.names(coef_mods_table) <- NULL
coef_mods_table[,2:5] <- round(coef_mods_table[,2:5], 3)
coef_mods_table <- coef_mods_table[,c(1, 6, 2, 4, 5)]
write_csv(coef_mods_table, "results/all_biomass_model_coef.csv")

# model weights ####

(b_weights <- loo_model_weights(
  list(global1 = loo_1,
       T_nutr2 = loo_2,
       Clim = loo_3,
       T_can = loo_4,
       T_only = loo_5,
       resource = loo_6)))

# Rerun best model with prior samples and save
mod_best <- update(mod6,
                   sample_prior = TRUE)
saveRDS(mod_best, "results/log_mg_mat_can_brms.RDS")
#mod_best <- readRDS("results/log_mg_mat_can_brms.RDS")



# summary out puts --------------------------------------------------------

beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}

temp_beta_log_mg <- bind_rows(fixef(mod1)["mat.c",],
                              fixef(mod2)["mat.c",],
                              fixef(mod3)["mat.c",],
                              fixef(mod4)["mat.c",],
                              fixef(mod5)["mat.c",],
                              fixef(mod6)["mat.c",])

temp_beta_log_mg$prob_less_0 <- 
  c(beta_0(mod1, "b_mat.c")$less,
    beta_0(mod2, "b_mat.c")$less,
    beta_0(mod3, "b_mat.c")$less,
    beta_0(mod4, "b_mat.c")$less,
    beta_0(mod5, "b_mat.c")$less,
    beta_0(mod6, "b_mat.c")$less)

temp_beta_log_mg$mod <- c("Global",
                          "Nutrient",
                          "Climate",
                          "Canopy",
                          "MAT.c",
                          "Resources")
temp_beta_log_mg$weight <- round(b_weights, 3)
saveRDS(temp_beta_log_mg, "results/log_mg_mat_coef_table.RDS")

# "best" model coefficients, probability < 0
fixef(mod6)
beta_0(mod6, "b_mat.c")$less

# conditional effects plot
plot(conditional_effects(mod6), points = TRUE)

# summary of best model ####

# probability that ISD exponent is negatively related to mat.c
beta_0(mod_best, "b_mat.c")$less

# Range of site-specific median ISD exponents
mod_best %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c,
               b_canopy) %>%
  filter(.iteration < 1000) %>%
  left_join(abiotic_s[,c("mat.c", "tdn", "tdp", "canopy", "siteID")]) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID + r_year) + #intercept
      b_mat.c * mat.c + b_canopy*canopy) %>%
  group_by(siteID) %>%
  summarize(med_biom = median(fitted_b)) %>%
  arrange(med_biom) %>%
  slice(c(1, n()))

# Save data for plots -----------------------------------------------------

# create plots directory, if it doesn't already exist.
# this should have already been run in script 4, but including here in case people run scripts out of order. 
if(!dir.exists("plots")){
  dir.create("plots")
}

# plot fitted log 10 mg mat.c coefficient
(mg_dist_plot <- mod6 %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c,
               b_canopy,
               b_tdn,
               b_tdp) %>%
  left_join(abiotic_s) %>%
  mutate(
    fitted_g =
      (b_Intercept + r_siteID) + #intercept
      b_mat.c * mat.c + b_canopy*canopy + b_tdn*tdn + b_tdp*tdp) %>%
  ggplot(aes(x = exp(fitted_g),
             fill = mat.c,
             y = reorder(siteID, mat.c))) +
  stat_halfeye() +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw() +
  scale_x_log10() +
  labs(y = "Site",
       x = "Grams dry weight per m2") +
  NULL)


saveRDS(mg_dist_plot, "plots/log_mg_mat_canopy_post_dist.RDS")

posts_6 <- posterior_samples(mod6) %>% mutate(iter = 1:nrow(.)) %>% 
  expand_grid(mat.c = seq(min(d$mat.c),
                          max(d$mat.c), 
                          length.out = 10)) %>% 
  mutate(fitted_g = b_Intercept + b_mat.c*mat.c) %>% 
  group_by(mat.c) %>% 
  summarize(median = median(exp(fitted_g)),
            upper = quantile(exp(fitted_g), probs = 0.975),
            lower = quantile(exp(fitted_g), probs = 0.025))

mg_fit_mat <- posts_6 %>% 
  ggplot(aes(
    x = mat.c, 
    y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.7) +
  scale_fill_manual(values = c("gray80")) +
  # scale_fill_grey()+#start = 0.2, end = 0.7) +
  geom_point(data = d, aes(y = biomass_g, color = mat.c),
             size = 2.5, 
             alpha = 0.8,
             data = d) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw()

saveRDS(mg_fit_mat, "plots/log_mg_mat.RDS")


mg_fit_can <- d %>%
  modelr::data_grid(
    siteID = siteID,
    year = year,
    mat.c = 0,
    canopy = seq(min(canopy),
                max(canopy), 
                length.out = 10)) %>%
  add_fitted_draws(mod6) %>%
  ggplot(aes(
    x = canopy, 
    y = log_mg)) +
  stat_lineribbon(aes(y = .value),
                  .width = c(0.95),
                  alpha = 0.7
  ) +
  scale_fill_manual(values = c("gray80")) +
  # scale_fill_grey()+#start = 0.2, end = 0.7) +
  geom_point(aes(color = canopy),
             size = 2.5, 
             alpha = 0.8,
             data = d) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw()

saveRDS(mg_fit_can, "plots/log_mg_can.RDS")
