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

biomass_mean <- biomass %>%
  group_by(ID) %>%
  summarize(u_biomass = mean(sample_biomass),
            sd_biomass = sd(sample_biomass))

biomass_mean <- left_join(biomass_mean, ID_key)
biomass_mean <- left_join(biomass_mean, site.info)
saveRDS(biomass_mean, "data/mean_biomass.RDS")
