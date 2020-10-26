# 3_community_biomass
# Mean annual temperature Celsius

library(tidyverse)
library(brms)
library(janitor)
library(ggridges)

# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS")

# read in site info, includes mean annual temp (mat.c)
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,10)]
names(site.info) <- c("siteID", "mat.c")

# wrangle dw to match degree days object
dw <- dw %>%
  separate(collectDate,
           c("year", "month", "day", "hour", "min", "sec")) %>%
  select(-hour, -min, -sec)

# remove lake and large river sites
lakes <- c("BARC", "CRAM", "LIRO", "PRLA", "PRPO",
           "SUGG", "TOOK", "BLWA", "FLNT", "TOMB")

dw <- dw %>%
  filter(!siteID %in% lakes)

# total biomass within a samples
biomass <- dw %>%
  select(siteID, year, month, day, sampleID,
         dw, no_m2) %>%
  group_by(siteID, sampleID, year, month, day) %>%
  mutate(dw_density = dw * no_m2) %>%
  summarize(sample_biomass = sum(dw_density, na.rm = TRUE)) %>%
  ungroup()

# add mat.c data
biomass <- inner_join(biomass,
                      site.info)
saveRDS(biomass, "data/sample_biomass_data.RDS")

# fit bayes
gamma.mod <- brm(data = biomass,
                 sample_biomass ~ mat.c +
                   (1 |siteID) + (1|year), 
                 family = Gamma(link = "log"),
                 prior = c(
                   prior(normal(0,0.1), class = "b"),
                   prior(normal(7, 1), class = "Intercept"),
                   prior(exponential(2), class="sd")),
                 chains = 4, 
                 sample_prior = TRUE,
                 iter = 6000,
                 cores = 4,
                 control = list(adapt_delta = 0.99))

# save gamma model
saveRDS(gamma.mod, "data/biomass_gamma_mod.RDS")

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
sum(post.mod$b_mat.c < 0)/ nrow(post.mod)
# probability that slope is positive
sum(post.mod$b_mat.c > 0)/ nrow(post.mod)

# main figure -------------------------------------------------------------

# main figure ####
# data to condition on
newdat <- tibble(
  mat.c = seq(min(biomass$mat.c),
              max(biomass$mat.c),
              length.out = 100)) 

# posterior samples to estimate from
posts <- posterior_samples(gamma.mod) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  sample_n(1000) %>% 
  expand_grid(newdat) %>% 
  mutate(b = exp(b_intercept + b_mat_c*mat.c)) %>% 
  group_by(mat.c) %>% 
  summarize(median = median(b),
            lower = quantile(b, probs = 0.025),
            upper = quantile(b, probs = 0.975))

# plot
(plot_biomass <- posts %>% 
    ggplot(aes(x = mat.c, y = median)) +
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
    scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
                  labels = c("0.1", "1","10","100","1,000","10,000","100,000","1,000,000","10,000,000")) +
    geom_text(aes(x = -Inf, y = Inf, label = "b)"),
              hjust = -.2,
              vjust = 1.2)  +
    NULL)

# save panel for use in script 4 to make the plots
saveRDS(plot_biomass, file = "plots/mat_c_biomass.rds")


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
# add degree data
site_biomass_est <- site_biomass_est %>%
  inner_join(site.info) 
# calculate model estimated biomass
site_biomass_est <- site_biomass_est %>%
  mutate(
    bm_est = b_intercept + site_offset + year_offset +
      b_mat_c * mat.c,
    bm_raw = exp(bm_est)) %>%
  group_by(siteID) %>%
  mutate(median_b = median(bm_est, na.rm = TRUE)) %>%
  ungroup()


# raw biomass values on log 10 scale
(biomass_post_B_panel <- site_biomass_est %>% 
  filter(iter <= 2000) %>%
  ggplot() +
  geom_density_ridges_gradient(
    aes(
      x = bm_raw, 
      y = reorder(siteID, mat.c),
      group = interaction(siteID, mat.c),
      fill = mat.c),
    scale = 2,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = 2) +
  scale_fill_viridis_c(alpha = 0.5, option = "plasma") +
  theme_bw() +
  labs(y = "Site", x = "Biomass") +
  scale_x_log10() +
  #facet_wrap(.~year) + 
    geom_text(aes(x = 0, y = Inf, label = "b)"),
              hjust = -1.2,
              vjust = 1.2)  +
    NULL)

# save panel for use in script 4 to make the plots
saveRDS(biomass_post_B_panel,
        file = "plots/biomass_post_panel_B.rds")

# summarize posterior samples ####
# global median
site_biomass_est %>% 
  summarize(
    med_raw = median(bm_raw, na.rm = TRUE),
    l95_raw = quantile(bm_raw, probs = 0.025, na.rm = TRUE),
    u95_raw = quantile(bm_raw, probs = 0.975, na.rm = TRUE))

# site:degree day estimates
site_biomass_est %>% 
  group_by(siteID, mat.c) %>%
  summarize(
    med_raw = median(bm_raw, na.rm = TRUE),
    l95_raw = quantile(bm_raw, probs = 0.025, na.rm = TRUE),
    u95_raw = quantile(bm_raw, probs = 0.975, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(med_raw) %>%
  slice(c(1, n())) # select first and last row

