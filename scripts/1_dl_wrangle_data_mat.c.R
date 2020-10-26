# 1 - Download and wrangle Neon Data

library(neonUtilities)
library(tidyverse)

source("scripts/custom_functions.R")

# load Length Weight coefficient table (used in part C below)
coeff <- read.csv("data/LW_coeffs.csv")


# A) Download data ---------------------------------------------------------

# these are large files, and may take a while to fully download dependent on web traffic and connection speed. 

# download all macroinvertebrate colllection data from January 2017 to december 2019
macro <- loadByProduct(dpID = "DP1.20120.001",
                       startdate = "2017-01",
                       enddate = "2019-12",
                       check.size = FALSE,
                       nCores = 3)


# C) Wrangle macroinvertebrate data ------------------------------------

# add length weight coefficients by taxon
MN.lw <- LW_coef(x = macro$inv_taxonomyProcessed,
                 lw_coef = coeff,
                 percent = TRUE)
# 96% of observations have length weight equations

# # estimate dry weights based on length weight coefficients
MN <- est_dw(MN.lw, fieldData = macro$inv_fieldData)

# questionable measurements ####
# filter out individuals that were "damaged" and measurement was affected
MN.no.damage <- MN.lw %>%
  filter(!str_detect(sampleCondition, "measurement")) %>%
  est_dw( fieldData = macro$inv_fieldData)

# save object for analysis in script 2 and 3
saveRDS(MN.no.damage, "data/macro_dw.RDS")


# Data summary ------------------------------------------------------------


# what percent of data is damaged?
message(paste(
  round(nrow(MN.no.damage) / nrow(MN.lw), 4)*100,
  "percent of individuals were not damaged"))

# 90% not damaged


# total number of measured individuals
length(na.omit(MN.no.damage$sizeClass))
# > 82 k individual measurements

# clean up MN.no.damage object
MN.no.damage <- MN.no.damage %>%
  # calculate total count / per m2 for each size class
  group_by(siteID, collectDate, dw) %>%
  summarise(count = sum(no_m2)) %>%
  ungroup() %>%
  # remove counts that are NA
  na.omit() %>%
  # duplicate number of rows based on count
  #i.e. if count = 10, duplicate that row 10 times
  uncount(count)

nrow(MN.no.damage)
# 18 million rows


