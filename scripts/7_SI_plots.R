# SI figures

# load libraries
library(brms)
library(tidyverse)
library(janitor)
library(cowplot)

# create plots directory, if it doesn't already exist.
# this should have already been run in script 4, but including here in case people run scripts out of order. 
if(!dir.exists("plots")){
  dir.create("plots")
}

# load data and models 
b.mod <- readRDS("results/b_mat_c_brms.RDS")
b.mod.data <- readRDS("data/MLEbins.RDS")
b.mod.data$mat.c <- scale(b.mod.data$mat.c)

ID_mat <- b.mod.data %>%
  select(ID, mat.c)


biomass.clim <- readRDS("models_jsw/mod3.rds")
biomass.resource <- readRDS("models_jsw/mod6.rds")
biomass.mod.data <- biomass.clim$data %>% left_join(biomass.resource$data)

abiotic <- read.csv("data/abiotic.csv")[,c(1,5)]
names(abiotic) <- c("siteID", "canopy")
abiotic$canopy <- scale(abiotic$canopy)



#prior v posterior ----------------------------

draws_b_long <- posterior_samples(b.mod) %>%
  as_tibble() %>%
  clean_names() %>% 
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_site", "r_year", "lp"))) %>% 
  pivot_longer(cols = -iter) %>% 
  mutate(parameter = case_when(
    name == "b_intercept" ~ "intercept",
    name == "b_mat_c" ~ "b", 
    name == "sd_site_id_intercept" ~ "sd_site_id",
    name == "sd_year_intercept" ~ "sd_year",
    grepl("sigma", name) ~ "sigma",
    name == "prior_intercept" ~ "intercept",
    name == "prior_b" ~ "b",
    name == "prior_sd_site_id" ~ "sd_site_id",
    name == "prior_sd_year" ~ "sd_year"),
    model = case_when(grepl("prior", name) ~ "prior", TRUE ~ "posterior")) %>% 
  select(-name)

prior_post_b <- draws_b_long %>% 
  mutate(parameter = case_when(
    parameter == "b" ~ "\u03b2", 
    parameter == "intercept" ~ "\u03b1",
    parameter == "sd_site_id" ~ "\u03c3_s",
    parameter == "sd_year" ~ "\u03c3_y",
    parameter == "sigma" ~ "\u03c3")) %>% 
  ggplot(aes(x = value, y = ..scaled.., fill = model)) + 
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("black", "NA")) +
  labs(y = "Scaled Density", 
       x = "Parameter Value") +
  theme_classic()


biomass.clim_post <- posterior_samples(biomass.clim) %>% mutate(model = "Climate Model",
                                                                prior_intercept = rnorm(nrow(.), 0, 2),
                                                                prior_b_mat.c = rnorm(nrow(.), 0, 1),
                                                                prior_b_map.mm = rnorm(nrow(.), 0, 1),
                                                                prior_sd_site = rexp(nrow(.), 1),
                                                                prior_sd_year = rexp(nrow(.), 1),
                                                                prior_shape = rgamma(nrow(.), 0.01, 0.01))

biomass.resource_post <- posterior_samples(biomass.resource) %>% mutate(model = "Resource Model",
                                                                prior_intercept = rnorm(nrow(.), 0, 2),
                                                                prior_b_mat.c = rnorm(nrow(.), 0, 1),
                                                                prior_b_tdn = rnorm(nrow(.), 0, 1),
                                                                prior_b_tdp = rnorm(nrow(.), 0, 1),
                                                                prior_b_canopy =rnorm(nrow(.), 0, 1),
                                                                prior_sd_site = rexp(nrow(.), 1),
                                                                prior_sd_year = rexp(nrow(.), 1),
                                                                prior_shape = rgamma(nrow(.), 0.01, 0.01))


biomass_posts <- bind_rows(biomass.clim_post, biomass.resource_post) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_site", "r_year", "lp"))) %>% 
  pivot_longer(cols = c(-model, -iter))

draws_dm_long <-  biomass_posts %>% 
  mutate(parameter = case_when(
    name == "sd_site_id_intercept" ~ "\u03c3_s",
    name == "sd_year_intercept" ~ "\u03c3_y",
    grepl("mat.c", name) ~ "\u03b2_mat_c",
    grepl("intercept", name) ~ "\u03b1",
    grepl("mm", name) ~ "\u03b2_map.mm",
    grepl("tdn", name) ~ "\u03b2_tdn",
    grepl("tdp", name) ~ "\u03b2_tdp",
    grepl("canopy", name) ~ "\u03b2_canopy",
    grepl("shape", name) ~ "shape",
    name == "prior_sd_site" ~ "\u03c3_s",
    name == "prior_sd_year" ~ "\u03c3_y",
    TRUE ~ name),
    post_prior = case_when(grepl("prior", name) ~ "prior", TRUE ~ "posterior")) %>% 
  select(-name)



prior_post_dm <- draws_dm_long %>% 
  ggplot(aes(x = value, y = ..scaled.., fill = post_prior)) + 
  geom_density() +
  facet_grid(parameter ~ model, scales = "free") +
  scale_fill_manual(values = c("black", "NA")) +
  labs(y = "Scaled Density", 
       x = "Parameter Value") +
  theme_classic() +
  xlim(NA, 4) +
  theme(axis.text.y = element_text(size = 8))

ggsave(prior_post_dm, file = "plots/prior_post_dm.jpg", dpi = 500, width = 6, height = 7.5)

# Figures S2 and S3 -------------------------------------------------------


# figure S3
ggsave(prior_post_b,
       file = "plots/SI3_prior_post_b.jpg",
       dpi = 500,
       width = 7,
       height = 4)

#figure S4
ggsave(prior_post_dm,
       file = "plots/SI4_prior_post_dm.jpg",
       dpi = 500,
       width = 7,
       height = 4)



# Figure S2 ---------------------------------------------------------------

# prepare data for panel a and b
draws_b_wide <- draws_b_long %>%
  pivot_wider(names_from = "parameter",
              values_from = "value")

mean_site_sds <- draws_b_wide %>%
  group_by(model) %>%
  summarize(mean_site_sds = mean(sd_site_id))
mean_year_sds <- draws_b_wide %>%
  group_by(model) %>%
  summarize(mean_year_sds = mean(sd_year))

#add predictor variable(s)
draws_b_withx <- draws_b_wide %>%
  group_by(model) %>% 
  sample_n(1000) %>%  #randomly choose 2000 iterations from each group to save space
  expand_grid(mat.c = unique(b.mod.data$mat.c)) %>% #add temperature
  left_join(b.mod.data %>%
              distinct(mat.c, siteID, year)) #add info for site and year


#calculate random offsets for intercepts
site_offsets <- draws_b_withx %>%
  distinct(model, iter) %>% 
  left_join(mean_site_sds) %>% 
  mutate(site_offset = rnorm(nrow(.), 0, mean_site_sds))

year_offsets <- draws_b_withx %>%
  distinct(model, iter) %>% 
  left_join(mean_year_sds) %>% 
  mutate(year_offset = rnorm(nrow(.), 0, mean_year_sds))

#solve model equation at each iteration
post_prior_predict <- draws_b_withx %>% 
  left_join(site_offsets) %>% 
  left_join(year_offsets) %>% 
  mutate(fitted = intercept + b*mat.c, #solve equation at each iteration. log-link is indicated by the exp().
         pred = case_when(model == "prior" ~ intercept + site_offset + year_offset + b*mat.c,
                       TRUE ~ fitted)) #same but with random offsets added.


# Figure S2 panels a and b
plot_slope_postprior <- post_prior_predict %>% 
  #filter(iter < 1000) %>%
  mutate(model = fct_relevel(model, "prior")) %>% 
  mutate(model = case_when(model == "prior"~ "a) Prior", TRUE ~ "b) Posterior")) %>% 
  ggplot(aes(x = mat.c, y = pred)) + 
  geom_line(aes(group = interaction(model,b)), alpha = 0.1) +
  facet_wrap(~model) +
  guides(color = F) +
  labs(y = expression(lambda),
       x = "Std(Mean Annual Temperature") +
  geom_point(
    data = b.mod.data %>%
      mutate(model = "b) Posterior",
             model = fct_relevel(model,"b) Posterior",
                                 after = 2)),
    aes(y = b),
    size = 0.5) +
  theme_classic() +
  NULL

# prepare data for panel c and d
post_prior_predict <- posterior_samples(biomass.mod) %>% select(!starts_with("r_")) %>% 
  mutate(iter = 1:nrow(.)) %>% 
  sample_n(1000) %>% 
  expand_grid(mat.c = biomass.mod.data$mat.c) %>% 
  mutate(prior_y = exp(prior_Intercept + prior_b*mat.c),
         posterior_y = exp(b_Intercept + b_mat.c*mat.c)) %>% 
  pivot_longer(cols = c(prior_y, posterior_y)) %>% 
  separate(name, c("model", "y"))


# fig S1 panel c and d
plot_biomass_postprior <- post_prior_predict %>%
  #filter(iter<= 1000) %>%
  mutate(model = fct_relevel(model, "prior")) %>% 
  mutate(model = 
           case_when(model == "prior" ~ "c) Prior",
                     TRUE ~ "d) Posterior")) %>% 
  ggplot(aes(x = mat.c,
             y = value)) + 
  geom_line(aes(group = iter),
            alpha = 0.1) +
  facet_wrap(~model) +
  guides(color = F) +
  labs(y = expression(paste("Macroinvertebrate dry mass (g/",m^2,")")),
       x = "Std(Mean Annual Temperature") +
  geom_point(
    data = biomass.mod.data %>%
      mutate(model = "d) Posterior",
             model = fct_relevel(model,"d) Posterior",
                                 after = 2)),
    aes(y = biomass_g),
    size = 0.5) +
  theme_classic() +
  scale_y_log10() +
  NULL

# put fig S1 panels together
prior_post_preds <- plot_grid(plot_slope_postprior,
                              plot_biomass_postprior,
                              ncol = 1,
                              align = "v")

# save figure S1
ggsave(prior_post_preds,
       file = "plots/SI2_prior_post_preds.jpg",
       dpi = 500,
       width = 7,
       height = 7)


# Posterior Predictive Checks # Figs S4 and S5 ---------------------------------------------

pp_bmod <- pp_check(b.mod, type = "boxplot")
ggsave(pp_bmod, file = "plots/SI5_b_pp.jpg")

pp_clim <- pp_check(biomass.clim, type = "boxplot") +
  labs(title = "a) Climate Model")
pp_clim
pp_resource <- pp_check(biomass.resource, type = "boxplot") +
  labs(title = "b) Resources Model")

pp_climresource <- plot_grid(pp_clim, pp_resource, nrow = 2, align = "h")

ggsave(pp_climresource, file = "plots/SI6_biomass_pp.jpg", width = 5, height = 5)




# Prior Sensitivity -------------------------------------------------------

# slope model
# update the slope model by halving the SD prior values

b.mod_sd0.5 <- readRDS("results/b.mod_sd0.5.rds")
b.mod_sd2 <- readRDS("results/b.mod_sd2.rds")
biomass.clim_sd0.5 <- readRDS("results/biomass.clim_sd0.5.rds")
biomass.clim_sd2 <- readRDS("results/biomass.clim_sd2.rds")
biomass.res_sd0.5 <- readRDS("results/biomass.res_sd0.5.rds")
biomass.res_sd2 <- readRDS("results/biomass.res_sd2.rds")

# b.mod_sd0.5 <- brm(data = b.mod.data,
#                    b ~ mat.c +
#                      (1 |siteID) + (1|year), 
#                    family = gaussian(),
#                    prior =
#                      c(prior(normal(0,0.12),
#                              class = "b"),
#                        prior(normal(-1.5, .5), 
#                              class = "Intercept"),
#                        prior(exponential(2),
#                              class="sd"),
#                        prior(exponential(2),
#                              class="sigma")),
#                    chains = 4, 
#                    sample_prior = FALSE,
#                    iter = 1000,
#                    cores = 4)
# 
# saveRDS(b.mod_sd0.5, file = "results/b.mod_sd0.5.rds")
# 
# # update the slope model by doubling the SD prior values
# b.mod_sd2 <- brm(data = b.mod.data,
#                  b ~ mat.c +
#                    (1 |siteID) + (1|year), 
#                  family = gaussian(),
#                  prior =
#                    c(prior(normal(0,0.5),
#                            class = "b"),
#                      prior(normal(-1.5, 2), 
#                            class = "Intercept"),
#                      prior(exponential(2),
#                            class="sd"),
#                      prior(exponential(2),
#                            class="sigma")),
#                  chains = 4, 
#                  sample_prior = FALSE,
#                  iter = 1000,
#                  cores = 4)
# 
# saveRDS(b.mod_sd2, file = "results/b.mod_sd2.rds")
# 
# 
# # biomass model
# # update the biomass model by halving the SD prior values
# biomass.clim_sd0.5 <- brm(data = biomass.mod.data,
#                           biomass_g ~ mat.c + map.mm +
#                             (1 |siteID) + (1|year), 
#                           family = Gamma(link = "log"),
#                           prior =
#                             c(prior(normal(0, 0.5),
#                                     class = "b"),
#                               prior(normal(0,1), 
#                                     class = "Intercept"),
#                               prior(exponential(1),
#                                     class="sd")),
#                           chains = 4, 
#                           sample_prior = FALSE,
#                           iter = 1000,
#                           cores = 4)
# 
# saveRDS(biomass.clim_sd0.5, file = "results/biomass.clim_sd0.5.rds")
# 
# # update the biomass model by doubling the SD prior values
# biomass.clim_sd2 <- brm(data = biomass.mod.data,
#                           biomass_g ~ mat.c + map.mm +
#                             (1 |siteID) + (1|year), 
#                           family = Gamma(link = "log"),
#                           prior =
#                             c(prior(normal(0, 2),
#                                     class = "b"),
#                               prior(normal(0,4), 
#                                     class = "Intercept"),
#                               prior(exponential(1),
#                                     class="sd")),
#                           chains = 4, 
#                           sample_prior = FALSE,
#                           iter = 1000,
#                           cores = 4)
# 
# saveRDS(biomass.clim_sd2, file = "results/biomass.clim_sd2.rds")
# 
# 
# # update the biomass model by halving the SD prior values
# biomass.res_sd0.5 <- brm(data = biomass.mod.data,
#                           biomass_g ~ mat.c + map.mm + tdn + tdp + canopy +
#                             (1 |siteID) + (1|year), 
#                           family = Gamma(link = "log"),
#                           prior =
#                             c(prior(normal(0, 0.5),
#                                     class = "b"),
#                               prior(normal(0,1), 
#                                     class = "Intercept"),
#                               prior(exponential(1),
#                                     class="sd")),
#                           chains = 4, 
#                           sample_prior = FALSE,
#                           iter = 1000,
#                           cores = 4)
# 
# saveRDS(biomass.res_sd0.5, file = "results/biomass.res_sd0.5.rds")
# 
# # update the biomass model by doubling the SD prior values
# biomass.res_sd2 <- brm(data = biomass.mod.data,
#                         biomass_g ~ mat.c + map.mm + tdn + tdp + canopy +
#                           (1 |siteID) + (1|year), 
#                         family = Gamma(link = "log"),
#                         prior =
#                           c(prior(normal(0, 2),
#                                   class = "b"),
#                             prior(normal(0,4), 
#                                   class = "Intercept"),
#                             prior(exponential(1),
#                                   class="sd")),
#                         chains = 4, 
#                         sample_prior = FALSE,
#                         iter = 1000,
#                         cores = 4)
# 
# saveRDS(biomass.res_sd2, file = "results/biomass.res_sd2.rds")
# 
# 
# 



# compare posterior samples
# slope model
posts_b_model <- posterior_samples(b.mod) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "b exponent model",
         version = "Temperature")

posts_b_sd0.5 <- posterior_samples(b.mod_sd0.5) %>%
  clean_names() %>%
  mutate(model = "sdx0.5",
         response = "b exponent model",
         version = "Temperature")

posts_b_sd2 <- posterior_samples(b.mod_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "b exponent model",
         version = "Temperature")

# biomass models
posts_dmmodel_res <- posterior_samples(biomass.clim) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "Community Biomass",
         version = "Climate")

posts_dmsd0.5_res <- posterior_samples(biomass.clim_sd0.5) %>%
  clean_names() %>% 
  mutate(model = "sdx0.5",
         response = "Community Biomass",
         version = "Climate")

posts_dmsd2_res <- posterior_samples(biomass.clim_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "Community Biomass",
         version = "Climate")


posts_dmmodel_clim <- posterior_samples(biomass.resource) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "Community Biomass",
         version = "Resources")

posts_dmsd0.5_clim <- posterior_samples(biomass.res_sd0.5) %>%
  clean_names() %>% 
  mutate(model = "sdx0.5",
         response = "Community Biomass",
         version = "Resources")

posts_dmsd2_clim <- posterior_samples(biomass.res_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "Community Biomass",
         version = "Resources")


all_b_sens <- bind_rows(posts_b_model,
                        posts_b_sd0.5,
                        posts_b_sd2,
                        posts_dmmodel_clim,
                        posts_dmsd0.5_clim,
                        posts_dmsd2_clim,
                        posts_dmmodel_res,
                        posts_dmsd0.5_res,
                        posts_dmsd2_res) %>% 
  select(!contains(c("site_id", "lp", "year"))) %>% 
  pivot_longer(cols = c(-model, -response, -version)) %>% 
  mutate(wrap = paste(response, "-", version))

(prior_sens_plot <- all_b_sens %>% 
  mutate(name = case_when(name == "b_mat_c" ~ "\u03b2_mat", 
                          name == "b_intercept" ~ "\u03b1",
                          name == "b_canopy" ~ "\u03b2_canopy",
                          name == "b_tdn" ~ "\u03b2_tdn",
                          name == "b_tdp" ~ "\u03b2_tdp",
                          name == "b_map_mm" ~ "\u03b2_map_mm")) %>% 
  ggplot(aes(x = value,
             color = model,
             y = ..scaled..)) +
  geom_density() + 
  facet_grid(name~wrap,
             scales = "free") +
  scale_color_brewer(type = "qual",
                     palette = 7) +
  theme_classic() +
  labs(y = "Scaled Density",
       x = "Parameter Value"))


ggsave(prior_sens_plot,
       file = "plots/SI7_prior_sens.jpg",
       dpi = 500,
       width = 8,
       height = 7)

# Fig S8 ####
# effect of canopy cover on log10 biomass in streams across the NEON sites

s8_plot <- readRDS("plots/log_mg_can.RDS")

s8_plot <- s8_plot +
  theme_bw() +
  labs(y = expression(
    "Log10 Dry Mass mg/"~m^2),
    x = "Standardized Mean Canopy Cover %")

ggsave(s8_plot,
       file = "plots/SI8_biomass_canopy.jpg",
       dpi = 500,
       width = 7,
       height = 5)

# fig S9 ####
# mean body mass

#library(ggridges)


# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS")

# read in site info, includes mean annual temp (mat.c)
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,6,10)]
names(site.info) <- c("siteID","latitude", "mat.c")
ID_key <- readRDS("data/sample_ID_key.RDS")

ID_key <- left_join(ID_key, site.info)


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


# # add latitude
# biomass <- left_join(biomass,
#                      site.info)
biomass <- left_join(biomass,
                     ID_key)

biomass_mean <- biomass %>%
  group_by(ID) %>%
  summarize(u_biomass = mean(sample_biomass),
            sd_biomass = sd(sample_biomass))

biomass_mean <- left_join(biomass_mean,
                          ID_key)

dw <- dw %>%
  select(ID, dw, no_m2)

body <- dw %>% group_by(ID) %>%
  summarize(mean_body = weighted.mean(dw, no_m2))

d <- left_join(biomass_mean, body)

summary(lm(log10(mean_body) ~ scale(mat.c), data = d))
summary(lm((mean_body) ~ scale(mat.c), data = d))

# SI figure ####
# mean body size across temperature
s9_plot <- ggplot(d,
       aes(x = scale(mat.c),
           y = mean_body,
           color = scale(mat.c))) +
  geom_point(size = 3, alpha = 0.7)+
  scale_y_log10() +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Mean Body Size",
       x = "Standardized Mean Annual Temperature") +
  NULL

ggsave(s9_plot,
       file = "plots/SI9_body_biomass.jpg",
       dpi = 500,
       width = 7,
       height = 5)


# figure S10 --------------------------------------------------------------


# fig S10####
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
MLEbins <- readRDS("data/MLEbins.RDS")
MLEbins <- MLEbins[,c("siteID", "ID", "year", "sampleEvent", "b")]

biomass <- readRDS("data/mean_biomass.RDS")
biomass <- biomass[,c("siteID", "ID", "year",
                      "sampleEvent", "u_biomass")]

# join data sets
d <- left_join(MLEbins, biomass)
d <- left_join(d, abiotic_s)
d <- left_join(d, ID_key)

# log 10 mean dry weight estimate
d$log_mg <- log10(d$u_biomass)

s_10 <- ggplot(d,
       aes(y = log_mg,
           x = b,
           color = mat.c)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  labs(y = "Log10 Community biomass", 
       x = "ISD exponent") +
  #stat_smooth(aes(x = b, y = log_mg),method = "lm", inherit.aes = FALSE) +
  NULL

ggsave(s_10,
       "plots/potential_SI_biomass_by_ISD.png",
       width = 10,
       height = 8,
       dpi = 300)

summary(lm(log_mg ~ b, data = d))
summary(lm(b ~ log_mg, data = d))