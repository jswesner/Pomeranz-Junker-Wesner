# 3_community_biomass
# Mean annual temperature Celsius

library(tidyverse)
library(brms)
library(janitor)
library(ggridges)

# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS")

# read in site info, includes mean annual temp (mat.c)
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,6)]
names(site.info) <- c("siteID","latitude")
ID_key <- readRDS("data/sample_ID_key.RDS")

# # wrangle dw to match site.info object
# dw <- dw %>%
#   separate(
#     collectDate,
#     c("year", "month", "day",
#       "hour", "min", "sec")) %>%
#   select(-hour, -min, -sec)


# total biomass within a samples
biomass <- dw %>%
  select(siteID, collectDate,
         sampleID, dw, no_m2, ID) %>%
  group_by(siteID, collectDate, sampleID, ID) %>%
  mutate(dw_density = dw * no_m2) %>%
  summarize(
    sample_biomass = sum(dw_density,
                         na.rm = TRUE)) %>%
  ungroup()

# add latitude
biomass <- left_join(biomass,
                      site.info)
biomass <- left_join(biomass,
                     ID_key)
saveRDS(biomass, "data/sample_biomass_latitude.RDS")

# fit bayes
gamma.mod <- brm(data = biomass,
                 sample_biomass ~ latitude +
                   (1 |siteID) + (1|year), 
                 family = Gamma(link = "log"),
                 prior =
                   c(prior(normal(0,0.1),
                           class = "b"),
                     prior(normal(7, 1),
                           class = "Intercept"),
                     prior(exponential(2),
                           class="sd")),
                 chains = 4, 
                 sample_prior = TRUE,
                 iter = 6000,
                 cores = 4,
                 control = list(adapt_delta = 0.99))

# save gamma model
saveRDS(gamma.mod, "results/biomass_latitude_gamma.RDS")

# Deciding to make figures based off of mean biomass model. # if we confirm this, replace code specifying gamma model above and get rid of everyhting else. 
gamma.mod <- readRDS("results/u_biomass_lat_brm.RDS")
# gamma model output and summaries
gamma.mod
fixef(gamma.mod)
plot(conditional_effects(gamma.mod), points = TRUE)

# posterior predictive checks
# Supplemental figure S5 (saved in script 5)
pp_check(gamma.mod, type = "boxplot")

# slope positive or negative?
post.mod <- posterior_samples(gamma.mod)
# probability that slope is less than specified value
sum(post.mod$b_latitude < 0)/ nrow(post.mod)
# probability that slope is positive
sum(post.mod$b_latitude > 0)/ nrow(post.mod)

# main figure -------------------------------------------------------------

# main figure ####
# data to condition on
newdat <- tibble(
  latitude = seq(min(biomass$latitude),
              max(biomass$latitude),
              length.out = 100)) 

# posterior samples to estimate from
posts <- posterior_samples(gamma.mod) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  sample_n(1000) %>% 
  expand_grid(newdat) %>% 
  mutate(b = exp(b_intercept + 
                   b_latitude * latitude)) %>% 
  group_by(latitude) %>% 
  summarize(median = median(b),
            lower = quantile(b, probs = 0.025),
            upper = quantile(b, probs = 0.975))

# plot
(plot_biomass <- posts %>% 
    ggplot(aes(x = latitude, y = median)) +
    geom_ribbon(alpha = 0.2,
                aes(ymax = upper,
                    ymin = lower)) + 
    geom_line() +
    geom_point(data = biomass,
               aes(y = sample_biomass),
               size = 1.7) +
    labs(y = "b",
         x = "Mean annual temperature C") +
    # scale_color_brewer(type = "qual") +
    theme_classic() +
    # ylim(-5,5) +
    # scale_y_log10() +
    scale_y_log10(
      breaks =
        c(0.1, 1, 10, 100, 1000,
          10000, 100000, 1000000, 10000000),
      labels = c("0.1", "1","10","100","1,000",
                 "10,000","100,000","1,000,000",
                 "10,000,000")) +
    geom_text(aes(x = -Inf,
                  y = Inf,
                  label = "b)"),
              hjust = -.2,
              vjust = 1.2)  +
    NULL)

# save panel for use in script 4 to make the plots
saveRDS(plot_biomass, file = "plots/latitude_biomass.rds")


# posterior samples -------------------------------------------------------


# extract model posterior samples
site_biomass_est <- posterior_samples(gamma.mod) %>%
  clean_names() %>% # janitor package
  mutate(iter = 1:nrow(.)) %>% # add iteration number
  pivot_longer( # extract random effect of site
    contains("r_site"),
    names_to = "site_ranef", values_to = "site_offset") %>%
  pivot_longer(# extract random effect of year
    contains("r_year"), names_to = "year_ranef", values_to = "year_offset") %>%
  mutate( # add siteID and year column
    siteID = str_sub(site_ranef, 11,14),
    siteID = toupper(siteID),
    year = str_sub(year_ranef, 8, 11))
# add latitude 
site_biomass_est <- site_biomass_est %>%
  inner_join(site.info) 
# calculate model estimated biomass
site_biomass_est <- site_biomass_est %>%
  mutate(
    bm_est = 
      b_intercept + site_offset + year_offset +
      b_latitude * latitude,
    bm_raw = exp(bm_est)) %>%
  group_by(siteID) %>%
  mutate(median_b = median(bm_est,
                           na.rm = TRUE)) %>%
  ungroup()


# raw biomass values on log 10 scale
(biomass_post_B_panel <- site_biomass_est %>% 
  #filter(iter <= 10) %>%
  ggplot() +
  geom_density_ridges_gradient(
    aes(
      x = bm_raw, 
      y = reorder(siteID, latitude),
      group = interaction(siteID, latitude),
      fill = latitude),
    scale = 2,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = 2) +
  scale_fill_viridis_c(alpha = 0.5,
                       option = "plasma",
                       direction = -1) +
  theme_bw() +
  labs(y = "Site", x = "Biomass") +
  scale_x_log10() +
  #facet_wrap(.~year) + 
    geom_text(aes(x = 0,
                  y = Inf,
                  label = "b)"),
              hjust = -1.2,
              vjust = 1.2)  +
    NULL)

# save panel for use in script 4 to make the plots
saveRDS(biomass_post_B_panel,
        file =
          "plots/biomass_latitiude_post_panel_B.rds")

# summarize posterior samples ####
# global median
site_biomass_est %>% 
  summarize(
    med_raw = median(bm_raw, na.rm = TRUE),
    l95_raw = quantile(bm_raw, probs = 0.025,
                       na.rm = TRUE),
    u95_raw = quantile(bm_raw, probs = 0.975,
                       na.rm = TRUE))

# site:degree day estimates
site_biomass_est %>% 
  group_by(siteID, latitude) %>%
  summarize(
    med_raw = median(bm_raw, na.rm = TRUE),
    l95_raw = quantile(bm_raw, probs = 0.025,
                       na.rm = TRUE),
    u95_raw = quantile(bm_raw, probs = 0.975,
                       na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(med_raw) %>%
  slice(c(1, n())) # select first and last row



# models with sample event ------------------------------------------------

#gamma.mod <- readRDS("results/biomass_latitude_gamma.RDS")
# fit bayes
mod2 <- brm(data = biomass,
                 sample_biomass ~ latitude + sampleEvent +
                   (1 |siteID) + (1|year), 
                 family = Gamma(link = "log"),
                 prior =
                   c(prior(normal(0,0.1),
                           class = "b"),
                     prior(normal(7, 1),
                           class = "Intercept"),
                     prior(exponential(2),
                           class="sd")),
                 chains = 4, 
                 sample_prior = TRUE,
                 iter = 6000,
                 cores = 4,
                 control = list(adapt_delta = 0.99))
# save model 2
saveRDS(mod2, "results/biomass_latitude_season_brm.RDS")
#mod2 <- readRDS("results/biomass_latitude_season_brm.RDS")
# interaction model
mod3 <- brm(data = biomass,
            sample_biomass ~ latitude * sampleEvent +
              (1 |siteID) + (1|year), 
            family = Gamma(link = "log"),
            prior =
              c(prior(normal(0,0.1),
                      class = "b"),
                prior(normal(7, 1),
                      class = "Intercept"),
                prior(exponential(2),
                      class="sd")),
            chains = 4, 
            sample_prior = TRUE,
            iter = 6000,
            cores = 4,
            control = list(adapt_delta = 0.99))

# save model 3
saveRDS(mod3, "results/biomass_interaction_brm.RDS")
#mod3 <- readRDS("results/biomass_interaction_brm.RDS")

gamma.mod <- add_criterion(gamma.mod, "waic")
mod2 <- add_criterion(mod2, "waic")
mod3 <- add_criterion(mod3, "waic")
loo_compare(gamma.mod, mod2, mod3, criterion = "waic")

loomod <- loo(gamma.mod, reloo = F)
loomod2 <- loo(mod2, reloo = F)
loomod3 <- loo(mod3, reloo = F)
loo_model_weights(list(loomod, loomod2, loomod3))

fixef(gamma.mod)
#fixef(mod2)
fixef(mod3)



plot(conditional_effects(gamma.mod), points = TRUE)
plot(conditional_effects(mod2), points = TRUE)
plot(conditional_effects(mod3), points = TRUE)

library(future)
library(loo)
plan(multiprocess)
g_mod1_k <- kfold(gamma.mod, k = 10, chains = 1)
g_mod2_k <- kfold(mod2, k = 10, chains = 1)
g_mod3_k <- kfold(mod3, k = 10, chains = 1)

g_compare <- loo_compare(g_mod1_k, g_mod2_k, g_mod3_k)
saveRDS(g_compare, "results/sample_biomass_compare_kfold.RDS")

#loo_model_weights(list(g_mod1_k, g_mod2_k, g_mod3_k))


g_kfold_point <- cbind(
  g_mod1_k$pointwise[,"elpd_kfold"],
  g_mod2_k$pointwise[,"elpd_kfold"],
  g_mod3_k$pointwise[,"elpd_kfold"]
  )

stacking_weights(g_kfold_point)

# average biomass ---------------------------------------------------------



# fit bayes
mean.mod <- brm(data = biomass_mean,
                 u_biomass ~ latitude +
                   (1 |siteID) + (1|year), 
                 family = Gamma(link = "log"),
                 prior =
                   c(prior(normal(0,0.1),
                           class = "b"),
                     prior(normal(7, 1),
                           class = "Intercept"),
                     prior(exponential(2),
                           class="sd")),
                 chains = 4, 
                 sample_prior = TRUE,
                 iter = 6000,
                 cores = 4,
                 control = list(adapt_delta = 0.99))

saveRDS(mean.mod, "results/u_biomass_lat_brm.RDS")
mean.mod <- readRDS("results/u_biomass_lat_brm.RDS")
# sampleEvent main effect
mean.mod2 <- brm(data = biomass_mean,
                u_biomass ~ latitude + sampleEvent +
                  (1 |siteID) + (1|year), 
                family = Gamma(link = "log"),
                prior =
                  c(prior(normal(0,0.1),
                          class = "b"),
                    prior(normal(7, 1),
                          class = "Intercept"),
                    prior(exponential(2),
                          class="sd")),
                chains = 4, 
                sample_prior = TRUE,
                iter = 6000,
                cores = 4,
                control = list(adapt_delta = 0.99))
saveRDS(mean.mod2, "results/u_biomass_lat_season_brm.RDS")
# sampleEvent interaction
mean.mod3 <- brm(data = biomass_mean,
                 u_biomass ~ latitude * sampleEvent +
                   (1 |siteID) + (1|year), 
                 family = Gamma(link = "log"),
                 prior =
                   c(prior(normal(0,0.1),
                           class = "b"),
                     prior(normal(7, 1),
                           class = "Intercept"),
                     prior(exponential(2),
                           class="sd")),
                 chains = 4, 
                 sample_prior = TRUE,
                 iter = 6000,
                 cores = 4,
                 control = list(adapt_delta = 0.99))
saveRDS(mean.mod3, "results/u_biomass_interaction_brm.RDS")

mean.mod <- add_criterion(mean.mod, "waic")
mean.mod2 <- add_criterion(mean.mod2, "waic")
mean.mod3 <- add_criterion(mean.mod3, "waic")
loo_compare(mean.mod, mean.mod2, mean.mod3, criterion = "waic")

loo_m_mod <- loo(mean.mod, reloo = F)
loo_m_mod2 <- loo(mean.mod2, reloo = F)
loo_m_mod3 <- loo(mean.mod3, reloo = F)
loo_model_weights(list(loo_m_mod, loo_m_mod2, loo_m_mod3))

fixef(mean.mod)
fixef(mean.mod2)
fixef(mean.mod3)

# plot(conditional_effects(mean.mod), points = T)
# plot(conditional_effects(mean.mod2), points = T)
# plot(conditional_effects(mean.mod3), points = T)

library(future)
plan(multiprocess)
mod1_k <- kfold(mean.mod, k = 10, chains = 1)
mod2_k <- kfold(mean.mod2, k = 10, chains = 1)
mod3_k <- kfold(mean.mod3, k = 10, chains = 1)

u_compare <- loo_compare(mod1_k, mod2_k, mod3_k)
saveRDS(u_compare, "results/u_biomass_compare_kfold.RDS")

#loo_model_weights(list(mod1_k, mod2_k, mod3_k))

u_kfold_point <- cbind(
  mod1_k$pointwise[,"elpd_kfold"],
  mod2_k$pointwise[,"elpd_kfold"],
  mod3_k$pointwise[,"elpd_kfold"]
)

stacking_weights(u_kfold_point)

