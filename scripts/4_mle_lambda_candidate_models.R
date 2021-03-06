# 4 Bayesian models

# fit candidate models for ISD lambda exponent 
# compare models using bayesian stacking
# make plots for main text

library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)

# read in explanatory variables
abiotic <- read.csv("data/abiotic.csv")
# standardize variables
abiotic_s <- scale(abiotic[,2:9]) 
# remove scaled attributes (mean and SD for each column)
scaled_attr <- attributes(abiotic_s)
# make scaled values into data frame and add siteID
abiotic_s <- as.data.frame(abiotic_s)
abiotic_s$siteID <- abiotic$Site

# read in b-exponent and biomass data and wrangle objects
MLEbins <- readRDS("data/MLEbins.RDS")
MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

# biomass <- readRDS("data/mean_biomass_latitude.RDS")
# biomass <- biomass[,c("siteID", "ID",
#                       "sampleEvent", "u_biomass")]
# 
# # join data sets
# d <- left_join(MLEbins, biomass)
d <- left_join(MLEbins, abiotic_s)

# log 10 mean dry weight estimate
#d$log_mg <- log10(d$u_biomass)

# preliminary stuff ---------------------------------------------------


cor(abiotic$latitude, abiotic[,c(-1, -9)])
ab_cor <- cor(abiotic[,c(2,3,4,5,7,8)])
ifelse(abs(ab_cor)<0.5, 0, ab_cor)
max(abs(ifelse(abs(ab_cor)==1, 0, ab_cor)))

sapply(abiotic[,-1], range)
sapply(abiotic[,-1], quantile, probs = c(0.05, 0.95))



# MLEbins exponent models -------------------------------------------------


# global model using standardized variables
mod1 <- brm(data = d,
            b ~ mat.c + map.mm + tdn + tdp + canopy +
              (1 |siteID) + (1|year), 
            family = gaussian(),
            prior =
              # 95% of slope prior between -0.5 and 0.5
              c(prior(normal(0,0.25),
                      class = "b"),
                # 95% of intercept prior between 
                # -3.5 and 0.5
                # i.e., exponent value at mean values of
                # all standardized variables 
                prior(normal(-1.5, 1), 
                      class = "Intercept"),
                prior(exponential(2),
                      class="sd"),
                prior(exponential(2),
                      class="sigma")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99)
            )
# temp + nutrients
mod2 <- update(mod1, formula. = . ~ . -canopy -map.mm,
               cores = 4)
# "Climate model" temp + precipitation
mod3 <- update(mod1, formula. = . ~ . -canopy -tdn -tdp,
               cores = 4)
# temp + canopy
mod4 <- update(mod1,
               formula. = . ~ . -map.mm -tdn - tdp,
               cores = 4)
# just temperature
mod5 <- update(mod1,
               formula. = . ~ . -map.mm -tdn -tdp -canopy,
               cores = 4)
# resources - autochthonous resources = TDN +TDP
# Allochthonous resources = canopy
mod6 <- update(mod1,
               formula. = . ~ . -map.mm,
               cores = 4)


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
write_csv(coef_mods_table, "results/all_b_model_coef.csv")

# model weights ####
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

# note on warnings:
# reloo uses the future package internally, and the future package recently updated to throw a warning when random numbers are generated.
# It seems like the warnings can be ignored based on this post: https://discourse.mc-stan.org/t/psis-loo-warns-about-an-unreliable-value/19437/2

# Allegedly you should be able to set seed=TRUE for the future package, but can't figure out how to do it. 

# The results I got previously (before the update) were similar for now, I'm using as-is

(b_weights <- loo_model_weights(
  list(global1 = loo_1,
       T_nutr2 = loo_2,
       Clim = loo_3,
       T_can = loo_4,
       T_only = loo_5,
       resource = loo_6)))

round(fixef(mod1), 3)
round(fixef(mod2), 3)
round(fixef(mod3), 3)
round(fixef(mod4), 3)
round(fixef(mod5), 3)
round(fixef(mod6), 3)

# Rerun best model with prior samples
# "best" model is mod
mod_best <- update(mod5,
                   sample_prior = TRUE)


saveRDS(mod_best, "results/b_mat_c_brms.RDS")
mod_best <- readRDS("results/b_mat_c_brms.RDS")

# function to calculate probability that coef is < or > 0
beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}

# summary of candidate models ####

temp_beta_b <- bind_rows(fixef(mod1)["mat.c",],
                         fixef(mod2)["mat.c",],
                         fixef(mod3)["mat.c",],
                         fixef(mod4)["mat.c",],
                         fixef(mod5)["mat.c",],
                         fixef(mod6)["mat.c",])
temp_beta_b$prob_less_0 <- c(
  beta_0(mod1, "b_mat.c")$less,
  beta_0(mod2, "b_mat.c")$less,
  beta_0(mod3, "b_mat.c")$less,
  beta_0(mod4, "b_mat.c")$less,
  beta_0(mod5, "b_mat.c")$less,
  beta_0(mod6, "b_mat.c")$less)

temp_beta_b$mod <- c("Global",
                     "Nutrient",
                     "Climate",
                     "Canopy",
                     "MAT.c",
                     "Resources")
temp_beta_b$weight <- round(b_weights, 3)
saveRDS(temp_beta_b, "results/b_mat_coef_table.RDS")
temp_beta_b <- readRDS("results/b_mat_coef_table.RDS")

# summary of best model ####

# probability that ISD exponent is negatively related to mat.c
beta_0(mod_best, "b_mat.c")$less

# Range of site-specific median ISD exponents
mod_best %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c) %>%
  filter(.iteration < 1000) %>%
  left_join(abiotic_s[,c("mat.c", "siteID")]) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID + r_year) + #intercept
      b_mat.c * mat.c) %>%
  group_by(siteID) %>%
  summarize(med_isd = median(fitted_b)) %>%
  arrange(med_isd) %>%
  slice(c(1, n()))

# plots with tidybayes ----------------------------------------------------

# create plots directory, if it doesn't already exist. 
if(!dir.exists("plots")){
  dir.create("plots")
}

# plots with tidybayes
get_variables(mod_best)

# plot fitted b
mat.c <- data.frame(mat.c =unique(d$mat.c))

mat_raw <- abiotic %>% select(Site, mat.c) %>% rename(siteID = Site,mat_raw = mat.c) %>% left_join(abiotic_s %>% select(siteID, mat.c)) %>%
  select(-siteID) %>% 
  right_join(mat.c) %>% 
  distinct(mat.c, mat_raw)

(b_post_plot <- mod_best %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c) %>%
  filter(.iteration < 1000) %>%
  left_join(abiotic_s[,c("mat.c", "siteID")]) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID + r_year) + #intercept
      b_mat.c * mat.c) %>%
  left_join(mat_raw) %>% 
  ggplot(aes(x = fitted_b,
             fill = mat_raw,
             y = reorder(siteID, mat_raw))) +
  stat_halfeye() +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Site",
       x = "ISD exponent") +
  NULL)

#ggsave("plots/b_mat_c_post_dist.jpg")
saveRDS(b_post_plot, "plots/b_mat_c_post_dist.RDS")

# global median of exponent
mod_best %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c) %>%
  left_join(abiotic_s[,c("mat.c", "siteID")]) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID + r_year) +
      b_mat.c * mat.c) %>%
  ungroup() %>%
  summarize(est = median(fitted_b),
            lower = quantile(fitted_b, probs = 0.05),
            upper = quantile(fitted_b, probs = 0.95))

# Site-specific medians
mod_best %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c) %>%
  left_join(abiotic_s[,c("mat.c", "siteID")]) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID) +
      b_mat.c * mat.c) %>%
  ungroup() %>%
  group_by(siteID) %>%
  summarize(est = median(fitted_b),
            lower = quantile(fitted_b, probs = 0.05),
            upper = quantile(fitted_b, probs = 0.95)) %>%
  arrange(est) %>%
  slice(1, n())


post_fits <- fitted(mod_best, newdata = data.frame(mat.c = unique(d$mat.c)),
                    re_formula = NA) %>% as_tibble() %>% clean_names() %>% 
  mutate(mat.c = unique(d$mat.c)) %>% 
  left_join(mat_raw)

(b_fit_plot <- post_fits %>% 
    ggplot(aes(x = mat_raw)) + 
    geom_ribbon(aes(ymin = q2_5, ymax = q97_5), alpha = 0.2) + 
    geom_line(aes(y = estimate)) + 
    geom_point(data = d %>% left_join(mat_raw), aes(color = mat_raw, y = b),
             size = 2.5, 
             alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
    guides(color = F) +
  labs(y = "ISD exponent",
       x = expression("Mean Annual Temperature " ( degree*C)),
       color = "Temperature") +
  theme_bw() )

#ggsave("plots/b_mat_c.jpg")
saveRDS(b_fit_plot, "plots/b_mat_c.RDS")
