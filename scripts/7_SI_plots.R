# SI figures

# load libraries
library(brms)
library(tidyverse)
library(janitor)
library(cowplot)


# load data and models 
b.mod <- readRDS("results/b_mat_c_brms.RDS")
b.mod.data <- readRDS("data/MLEbins.RDS")
b.mod.data$mat.c <- scale(b.mod.data$mat.c)

ID_mat <- b.mod.data %>%
  select(ID, mat.c)


biomass.mod <- readRDS("results/log_mg_mat_can_brms.RDS")
biomass.mod.data <- readRDS("data/mean_biomass.RDS")

abiotic <- read.csv("data/abiotic.csv")[,c(1,5)]
names(abiotic) <- c("siteID", "canopy")
abiotic$canopy <- scale(abiotic$canopy)

biomass.mod.data <- left_join(biomass.mod.data, abiotic)
biomass.mod.data <- biomass.mod.data %>%
  left_join(ID_mat)


fixef(b.mod)
fixef(biomass.mod)



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



draws_dm_long <- posterior_samples(biomass.mod) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  select(!contains(c("r_site", "r_year", "lp"))) %>% 
  pivot_longer(cols = -iter) %>% 
  mutate(parameter = case_when(
    name == "b_intercept" ~ "intercept",
    name == "b_mat_c" ~ "b_mat_c",
    name == "b_canopy" ~ "b_canopy",
    name == "sd_site_id_intercept" ~ "sd_site_id",
    name == "sd_year_intercept" ~ "sd_year",
    grepl("shape", name) ~ "shape",
    name == "prior_intercept" ~ "intercept",
    name == "prior_b" ~ "b_mat_c",
    name == "prior_sd_site_id" ~ "sd_site_id",
    name == "prior_sd_year" ~ "sd_year"),
    model = case_when(grepl("prior", name) ~ "prior", TRUE ~ "posterior")) %>% 
  select(-name)

prior_b_2 <- draws_dm_long %>%
  filter(parameter == "b_mat_c",
         model == "prior") %>%
  mutate(parameter = case_when(
    parameter == "b_mat_c" ~ "b_canopy"
  ))

draws_dm_long <- bind_rows(draws_dm_long, prior_b_2)


prior_post_dm <- draws_dm_long %>% 
  mutate(parameter = case_when(
    parameter == "b_mat_c" ~ "\u03b2_mat_c",
    parameter == "b_canopy" ~ "\u03b2_canopy",
    parameter == "intercept" ~ "\u03b1",
    parameter == "sd_site_id" ~ "\u03c3_s",
    parameter == "sd_year" ~ "\u03c3_y",
    parameter == "shape" ~ "shape")) %>% 
  ggplot(aes(x = value, y = ..scaled.., fill = model)) + 
  geom_density() +
  facet_wrap(~parameter, scales = "free") +
  scale_fill_manual(values = c("black", "NA")) +
  labs(y = "Scaled Density", 
       x = "Parameter Value") +
  theme_classic()


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
draws_dm_wide <- draws_dm_long %>%
  pivot_wider(names_from = "parameter", values_from = "value")

mean_site_sds <- draws_dm_wide %>%
  group_by(model) %>%
  summarize(mean_site_sds = mean(sd_site_id))
mean_year_sds <- draws_dm_wide %>%
  group_by(model) %>%
  summarize(mean_year_sds = mean(sd_year))

#add predictor variable(s)
draws_dm_withx <- draws_dm_wide %>%
  group_by(model) %>% 
  sample_n(1000) %>%  #randomly choose 2000 iterations from each group to save space
  expand_grid(mat.c = unique(b.mod.data$mat.c),
              canopy = 0) %>% #add degree days
  left_join(b.mod.data %>% distinct(mat.c, siteID, year)) #add info for site and year


#calculate random offsets for intercepts
site_offsets <- draws_dm_withx %>%
  distinct(model, iter) %>% 
  left_join(mean_site_sds) %>% 
  mutate(site_offset = rnorm(nrow(.), 0, mean_site_sds))

year_offsets <- draws_dm_withx %>%
  distinct(model, iter) %>% 
  left_join(mean_year_sds) %>% 
  mutate(year_offset = rnorm(nrow(.), 0, mean_year_sds))

#solve model equation at each iteration
post_prior_predict <- draws_dm_withx %>% 
  left_join(site_offsets) %>% 
  left_join(year_offsets) %>% 
  mutate(fitted = exp(intercept + b_mat_c*mat.c + b_canopy*canopy), #solve equation at each iteration. log-link is indicated by the exp().
         pred = case_when(
           model == "prior" ~
             exp(intercept + site_offset + year_offset + 
                   b_mat_c*mat.c + b_canopy*canopy),
                          TRUE ~ fitted)) #same but with random offsets added.

# fig S1 panel c and d
plot_biomass_postprior <- post_prior_predict %>%
  #filter(iter<= 1000) %>%
  mutate(model = fct_relevel(model, "prior")) %>% 
  mutate(model = 
           case_when(model == "prior" ~ "c) Prior",
                     TRUE ~ "d) Posterior")) %>% 
  ggplot(aes(x = mat.c,
             y = pred)) + 
  geom_line(aes(group = interaction(model,b_mat_c)),
            alpha = 0.1) +
  facet_wrap(~model) +
  guides(color = F) +
  labs(y = expression(paste("Community Biomass (mgDM/",m^2,")")),
       x = "Std(Mean Annual Temperature") +
  geom_point(
    data = biomass.mod.data %>%
      mutate(model = "d) Posterior",
             model = fct_relevel(model,"d) Posterior",
                                 after = 2)),
    aes(y = log10(u_biomass)),
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

pp_check(b.mod, type = "boxplot")
ggsave("plots/SI5_b_pp.jpg")

pp_check(biomass.mod, type = "boxplot")
ggsave("plots/SI6_biomass_pp.jpg")




# Prior Sensitivity -------------------------------------------------------

# slope model
# update the slope model by halving the SD prior values
b.mod_sd0.5 <- update(b.mod,
                      b.mod,
                      prior = c(prior(normal(0,0.125),
                                      class = "b"),
                                prior(normal(-1.5,0.5),
                                      class = "Intercept"),
                                prior(exponential(1),
                                      class = "sigma"),
                                prior(exponential(1),
                                      class = "sd")),
                      iter = 1000,
                      cores = 2,
                      chains = 2)

# update the slope model by doubling the SD prior values
b.mod_sd2 <- update(b.mod,
                    b.mod,
                    prior = c(prior(normal(-0,0.5),
                                    class = "b"),
                              prior(normal(-1.5,2),
                                    class = "Intercept"),
                              prior(exponential(4),
                                    class = "sigma"),
                              prior(exponential(4),
                                    class = "sd")),
                    iter = 1000,
                    cores = 2,
                    chains = 2)



# biomass model
# update the biomass model by halving the SD prior values
biomass.mod_sd0.5 <- brm(log10(u_biomass) ~ mat.c + canopy +
                           (1|siteID) +
                           (1|year),
                       data = biomass.mod.data,
                       family = Gamma(link = "log"),
                       prior = c(prior(normal(0, 0.25),
                                       class = "b"),
                                 prior(normal(1.1, 0.125),
                                       class = "Intercept"),
                                 prior(gamma(0.01, 0.01),
                                       class = "shape"),
                                 prior(exponential(1),
                                       class = "sd")),
                       chains = 2,
                       cores = 2,
                       iter = 1000,
                       sample_prior = TRUE)

# update the biomass model by doubling the SD prior values
biomass.mod_sd2 <- brm(log10(u_biomass) ~ mat.c + canopy +
                         (1|siteID) +
                         (1|year),
                   data = biomass.mod.data,
                   family = Gamma(link = "log"),
                   prior = c(prior(normal(0, 1),
                                   class = "b"),
                             prior(normal(1.1, 0.5),
                                   class = "Intercept"),
                             prior(gamma(0.01, 0.01),
                                   class = "shape"),
                             prior(exponential(4),
                                   class = "sd")),
                   chains = 2,
                   cores = 2,
                   iter = 1000,
                   sample_prior = TRUE)


# compare posterior samples
# slope model
posts_b_model <- posterior_samples(b.mod) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "b exponent model")

posts_b_sd0.5 <- posterior_samples(b.mod_sd0.5) %>%
  clean_names() %>%
  mutate(model = "sdx0.5",
         response = "b exponent model")

posts_b_sd2 <- posterior_samples(b.mod_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "b exponent model")

# biomass model
posts_dmmodel <- posterior_samples(biomass.mod) %>%
  clean_names() %>%
  mutate(model = "original model",
         response = "Community Biomass")

posts_dmsd0.5 <- posterior_samples(biomass.mod_sd0.5) %>%
  clean_names() %>% 
  mutate(model = "sdx0.5",
         response = "Community Biomass")

posts_dmsd2 <- posterior_samples(biomass.mod_sd2) %>%
  clean_names() %>%
  mutate(model = "sdx2",
         response = "Community Biomass")


all_b_sens <- bind_rows(posts_b_model,
                        posts_b_sd0.5,
                        posts_b_sd2,
                        posts_dmmodel,
                        posts_dmsd0.5,
                        posts_dmsd2) %>% 
  pivot_longer(cols = c("b_intercept", "b_mat_c", "b_canopy"))

(prior_sens_plot <- all_b_sens %>% 
  mutate(name = case_when(name == "b_mat_c" ~ "\u03b2_mat", 
                          name == "b_intercept" ~ "\u03b1",
                          name == "b_canopy" ~ "\u03b2_canopy")) %>% 
  ggplot(aes(x = value,
             color = model,
             y = ..scaled..)) +
  geom_density() + 
  facet_wrap(response~name,
             scales = "free") +
  scale_color_brewer(type = "qual",
                     palette = 7) +
  theme_classic() +
  labs(y = "Scaled Density",
       x = "Parameter Value"))


ggsave(prior_sens_plot,
       file = "plots/SI7_prior_sens.jpg",
       dpi = 500,
       width = 7,
       height = 5)

