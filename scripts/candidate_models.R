# Bayesian models

# fit candidate models for ISD exponent and total log10 community biomass. 
# compare models using bayesian stacking
# make plots for main text

library(tidyverse)
library(brms)
library(tidybayes)
library(viridis)

# read in explanatory variables
abiotic <- read.csv("data/abiotic.csv")
# standardize variables
abiotic_s <- scale(abiotic[,2:8]) 
# remove scaled attributes (mean and SD for each column)
scaled_attr <- attributes(abiotic_s)
# make scaled values into data frame and add siteID
abiotic_s <- as.data.frame(abiotic_s)
abiotic_s$siteID <- abiotic$Site

# read in b-exponent and biomass data and wrangle objects
MLEbins <- readRDS("data/MLEbins.RDS")
MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

biomass <- readRDS("data/mean_biomass_latitude.RDS")
biomass <- biomass[,c("siteID", "ID",
                      "sampleEvent", "u_biomass")]

# join data sets
d <- left_join(MLEbins, biomass)
d <- left_join(d, abiotic_s)

# log 10 mean dry weight estimate
d$log_mg <- log10(d$u_biomass)

# preliminary stuff ---------------------------------------------------


cor(abiotic$latitude, abiotic[,-1])
ab_cor <- cor(abiotic[,c(2,4,6,7)])
ifelse(abs(ab_cor)<0.5, 0, ab_cor)
max(abs(ifelse(abs(ab_cor)==1, 0, ab_cor)))

sapply(abiotic[,-1], range)
sapply(abiotic[,-1], quantile, probs = c(0.05, 0.95))



# MLEbins exponent models -------------------------------------------------


# global model using standardized variables
mod1 <- brm(data = d,
            b ~ mat.c + tdn + tdp + canopy +
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
mod2 <- update(mod1, formula. = . ~ . -canopy,
               cores = 4)
# just temperature
mod5 <- update(mod1,
               formula. = . ~ .- tdn - tdp - canopy,
               cores = 4)
# temp + canopy
mod8 <- update(mod1,
               formula. = . ~ . -tdn - tdp,
               cores = 4)


loo_1 <- loo(mod1,
             reloo = TRUE,
             cores = 3)
loo_2 <- loo(mod2, reloo = TRUE,
             cores = 3)
loo_5 <- loo(mod5, reloo = TRUE,
             cores = 3)
loo_8 <- loo(mod8, reloo = TRUE,
             cores = 3)

b_weights <- loo_model_weights(list(global1 = loo_1,
                       T_nutr2 = loo_2,
                       Temp5 = loo_5,
                       T_c8 = loo_8))
round(fixef(mod1), 3)
round(fixef(mod2), 3)
round(fixef(mod8), 3)
round(fixef(mod5), 3)

# save best model
# model 5 with prior samples
mod5 <- brm(data = d,
            b ~ mat.c +
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
            sample_prior = TRUE,
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99)
)
saveRDS(mod5, "results/b_mat_c_brms.RDS")


beta_0 <- function(model, b_est){
  post <- posterior_samples(model)
  less <- sum(post[[b_est]] < 0)/ nrow(post)
  more <- sum(post[[b_est]] > 0)/ nrow(post)
  list(less = less, more = more)
}


temp_beta_b <- bind_rows(fixef(mod1)["mat.c",],
                              fixef(mod2)["mat.c",],
                              fixef(mod5)["mat.c",],
                              fixef(mod8)["mat.c",])
temp_beta_b$prob_less_0 <- c(
  beta_0(mod1, "b_mat.c")$less,
  beta_0(mod2, "b_mat.c")$less,
  beta_0(mod5, "b_mat.c")$less,
  beta_0(mod8, "b_mat.c")$less)
temp_beta_b$mod <- c("Resources", 
                     "Nutrient",
                     "MAT.c",
                     "Canopy")
temp_beta_b$weight <- round(b_weights[1:4], 3)
saveRDS(temp_beta_b, "results/b_mat_coef_table.RDS")

# plots with tidybayes ----------------------------------------------------

# plots with tidybayes
get_variables(mod5)

# plot densities of fixed effects
mod5 %>%
  gather_draws(#b_Intercept,
    b_mat.c) %>%
  #median_qi() %>%
  ggplot(aes(x = .value,
             y = .variable))+
  stat_halfeye()

# plot site specific intercepts
mod8 %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term]) %>%
  left_join(abiotic_s) %>%
  mutate(site_intercept = b_Intercept +
              r_siteID + r_year) %>%
  ggplot(
    aes(x = site_intercept,
        color = mat.c,
        y = interaction(
          year, reorder(siteID, mat.c))))+
  stat_pointinterval() +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  NULL

# plot fitted b exponent
mod8 %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c,
               b_canopy) %>%
  left_join(abiotic_s[,c("mat.c", "canopy", "siteID")]) %>%
  #mutate(site_intercept = b_Intercept + r_siteID) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID + r_year) + #intercept
      b_mat.c * mat.c + b_canopy * canopy) %>%
  ggplot(aes(x = fitted_b,
             #xmin = .lower, 
             #xmax = .upper,
             color = mat.c,
             y = interaction(year, reorder(siteID,
                                           mat.c))))+
  #stat_halfeye(adjust = 0.95)
  stat_pointinterval() +
  #facet_wrap(.~year, nrow = 3)
  NULL

# plot fitted b
b_post_plot <- mod5 %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c) %>%
  left_join(abiotic_s[,c("mat.c", "siteID")]) %>%
  mutate(
    fitted_b =
      (b_Intercept + r_siteID + r_year) + #intercept
      b_mat.c * mat.c) %>%
  ggplot(aes(x = fitted_b,
             fill = mat.c,
             y = reorder(siteID, mat.c)))+
  stat_halfeye(interval_color = "black",
               point_color = "black") +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Site",
       x = "ISD exponent") +
  NULL
ggsave("plots/b_mat_c_post_dist.jpg")
saveRDS(b_post_plot, "plots/b_mat_c_post_dist.RDS")

# global median of exponent
mod5 %>%
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
mod5 %>%
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
  group_by(siteID) %>%
  summarize(est = median(fitted_b),
            lower = quantile(fitted_b, probs = 0.05),
            upper = quantile(fitted_b, probs = 0.95)) %>%
  arrange(est) %>%
  slice(1, n())


b_fit_plot <- d %>%
  modelr::data_grid(
    siteID = siteID,
    year = year,
    #log_mg = log_mg,
    mat.c = seq(min(mat.c),
                max(mat.c), 
                length.out = 10)) %>%
  add_fitted_draws(mod5) %>%
  ggplot(aes(
    x = mat.c, 
    y = b)) +
  stat_lineribbon(aes(y = .value),
                  .width = c(0.5, 0.89),
                  alpha = 0.7
  ) +
  scale_fill_manual(values = c("gray80", "gray55")) +
  geom_point(aes(color = mat.c),
             size = 2.5, 
             alpha = 0.8,
             data = d) +
  scale_color_viridis_c(option = "plasma") +
  labs(y = "ISD exponent",
       x = "Standardized Mean Annual Temperature",
       color = "Temperature") +
  theme_bw() 

ggsave("plots/b_mat_c.jpg")
saveRDS(b_fit_plot, "plots/b_mat_c.RDS")

# biomass in grams --------------------------------------------------------



# community biomass models ------------------------------------------------

d$log_mg <- log10(d$u_biomass)

# global model using standardized variables
# temperature + nutrients + canopy
mod1 <- brm(data = d,
            log_mg ~ mat.c + tdn + tdp + canopy +
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              # 95% of slope prior between -1 and 1
              # response is on log10 scale
              # log-link exponentiates
              c(prior(normal(0,0.5),
                      class = "b"),
                # 95% of intercept prior between 
                # -3.5 and 0.5
                # i.e., exponent value at mean values of
                # all standardized variables 
                prior(normal(1.1, 0.25), 
                      class = "Intercept"),
                prior(exponential(2),
                      class="sd")),
            chains = 4, 
            sample_prior = FALSE,
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99)
)

# temp + nutrients
mod2 <- update(mod1, formula. = . ~ . -canopy)

# just temperature
mod5 <- update(mod1, formula. = . ~ .- tdn -
                 tdp - canopy)

# temp + canopy
mod8 <- update(mod1, formula. = . ~ . -tdn - tdp)

loo_1 <- loo(mod1,
             reloo = TRUE,
             cores = 3)
loo_2 <- loo(mod2, reloo = TRUE,
             cores = 3)

loo_5 <- loo(mod5, reloo = TRUE,
             cores = 3)

loo_8 <- loo(mod8, reloo = TRUE,
             cores = 3)

loo_model_weights(list(global1 = loo_1,
                       T_nutr2 = loo_2,
                       Temp5 = loo_5,
                       T_c8 = loo_8))

# update model 5 with prior samples
mgmod5 <- brm(data = d,
            log_mg ~ mat.c +
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              # 95% of slope prior between -1 and 1
              # response is on log10 scale
              # log-link exponentiates
              c(prior(normal(0,0.5),
                      class = "b"),
                # 95% of intercept prior between 
                # -3.5 and 0.5
                # i.e., exponent value at mean values of
                # all standardized variables 
                prior(normal(1.1, 0.25), 
                      class = "Intercept"),
                prior(exponential(2),
                      class="sd")),
            chains = 4, 
            sample_prior = TRUE,
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99)
)
# # save best model
# mod5 <- update(mod5,
#                chains = 4, 
#                iter = 6000,
#                cores = 4)
saveRDS(mgmod5, "results/log_mg_mat_c_brms.RDS")
#mod5 <- readRDS("results/log_mg_mat_c_brms.RDS")
temp_beta_log_mg <- bind_rows(fixef(mod1)["mat.c",],
                              fixef(mod2)["mat.c",],
                              fixef(mod5)["mat.c",],
                              fixef(mod8)["mat.c",])
temp_beta_log_mg$prob_less_0 <- 
  c(
  beta_0(mod1, "b_mat.c")$less,
  beta_0(mod2, "b_mat.c")$less,
  beta_0(mod5, "b_mat.c")$less,
  beta_0(mod8, "b_mat.c")$less)
temp_beta_log_mg$mod <- c(
  "Resources", 
  "Nutrient",
  "MAT.c",
  "Canopy")
saveRDS(temp_beta_log_mg, "results/log_mg_mat_coef_table.RDS")


plot(conditional_effects(mod5), points = TRUE)

# plot densities of fixed effect (b_mat.c)
mod5 %>%
  gather_draws(b_mat.c) %>%
  ggplot(aes(x = .value,
             y = .variable))+
  stat_halfeye() +
  theme_bw()+
  geom_vline(xintercept = 0,
             color = "red",
             lty = "dashed") +
  labs(y = "")

# plot fitted log 10 mg
mg_dist_plot <- mod5 %>%
  spread_draws(b_Intercept,
               r_siteID[siteID, term],
               r_year[year, term],
               b_mat.c) %>%
  left_join(abiotic_s[,c("mat.c", "siteID")]) %>%
  mutate(
    fitted_log_mg =
      (b_Intercept + r_siteID + r_year) + #intercept
      b_mat.c * mat.c) %>%
  ggplot(aes(x = fitted_log_mg,
             fill = mat.c,
             y = reorder(siteID, mat.c)))+
  stat_halfeye(interval_color = "black",
               point_color = "black") +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Site",
       x = "Log10 mg dry weight per m2") +
  NULL
ggsave("plots/log_mg_mat_c_post_dist.jpg")
saveRDS(mg_dist_plot, "plots/log_mg_mat_c_post_dist.RDS")

mg_fit_plot <- d %>%
  modelr::data_grid(
    siteID = siteID,
    year = year,
    #log_mg = log_mg,
    mat.c = seq(min(mat.c),
                max(mat.c), 
                length.out = 10)) %>%
  add_fitted_draws(mod5) %>%
  ggplot(aes(
    x = mat.c, 
    y = log_mg)) +
  stat_lineribbon(aes(y = .value),
                  .width = c(0.5, 0.89),
                  alpha = 0.7
                  ) +
  scale_fill_manual(values = c("gray80", "gray55")) +
  # scale_fill_grey()+#start = 0.2, end = 0.7) +
  geom_point(aes(color = mat.c),
             size = 2.5, 
             alpha = 0.8,
             data = d) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw()
saveRDS(mg_fit_plot, "plots/log_mg_mat_c.RDS")
ggsave("plots/log_mg_mat_c.jpg")

