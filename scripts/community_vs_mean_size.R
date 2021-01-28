# mean body mass

library(tidyverse)
library(brms)
library(janitor)
library(ggridges)


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


ggplot(d,
       aes(x = mean_body,
           y = u_biomass,
           color = mat.c)) +
  geom_point(size = 3, alpha = 0.7)+
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  stat_smooth(aes(x = mean_body, y = u_biomass),
              method = "lm", inherit.aes = FALSE) +
  geom_abline(slope = 0.25, intercept = 3.7,
              color = "red", size = 2) +
  NULL

summary(lm(log10(u_biomass) ~ log10(mean_body), data = d))
summary(lm((u_biomass) ~ (mean_body), data = d))

ggplot(d,
       aes(x = mean_body,
           y = u_biomass,
           color = mat.c,
           group = mat.c)) +
  geom_point(size = 3, alpha = 0.7)+
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(option = "plasma") +
  theme_bw() +
  stat_smooth(method = "lm", se = FALSE) +
  NULL
