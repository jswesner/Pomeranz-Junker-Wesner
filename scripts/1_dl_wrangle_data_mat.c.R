# 1 - Download and wrangle Neon Data

# this is the script that was used to download the data presented in the manuscript. 
# !!!!!! HOWEVER !!!!!!
# the NEON dataportal is live, and data may be updated or ammended in the future. I include this script here for future use, and to aid in the access to NEON data. 
# If you would like to recreate the results presented in the manuscript, I recommend you DO NOT RUN THIS SCRIPT!
# instead, you can work directly with the data used for the analysis by downloading it from [<insert data dryad URL here once archived>] and skipping to script 2_mle_bin_full_data.R

library(neonUtilities)
library(tidyverse)
library(lubridate)

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

MN.no.damage <- readRDS("data/macro_dw.RDS")

# remove lake and large river sites
lakes <- c("BARC", "CRAM", "LIRO", "PRLA", "PRPO",
           "SUGG", "TOOK", "BLWA", "FLNT", "TOMB")
MN.no.damage <- MN.no.damage %>%
  filter(!siteID %in% lakes)

# filter out NA values in dw
MN.no.damage <- MN.no.damage %>%
  filter(!is.na(dw), !is.na(no_m2)) 

########################################################
########################################################
# double check this!!!! ---------------------------------------------------
########################################################
########################################################


# when the data was downloaded in October 2020, some of the collections were conducted across 2 days. Here, I combine them so that there is one unique siteID:collectDate combination for each collection. This may change as NEON data is updated

# fix some collectDate 
# i.e. HOPB has 2 samples from same date (2017-4-12) but different times
# MAYF was collected on 2 dates in 2017-07, combining to one sample
# run the following and view the output in the console
MN.no.damage %>%
  filter(siteID == "HOPB" | siteID == "MAYF") %>%
  select(siteID, collectDate) %>%
  unique() %>%
  arrange(siteID, collectDate)

# note row 1 and 2 are the same date, but a different time, and line 12 and 13 were collected 2 days apart

# the following code chunks "fix" these so they are the same collectDate

# change HOPB 2017-04-12 13:23:00 --> 2017-04-12 16:00
MN.no.damage[
  MN.no.damage$siteID == "HOPB" &
    MN.no.damage$collectDate == 
    as.POSIXct("2017-04-12 13:23:00", tz = "GMT"),
  "collectDate"] <- as.POSIXct("2017-04-12 16:00:00", tz = "GMT")

# change MAYF 2017-07-18 14:51:00 to 2017-07-20 14:37:00
MN.no.damage[
  MN.no.damage$siteID == "MAYF" &
    MN.no.damage$collectDate == 
    as.POSIXct("2017-07-18 14:51:00", tz = "GMT"),
  "collectDate"] <- as.POSIXct("2017-07-20 14:37:00", tz = "GMT")

# make unique group index for siteID:collectDate combo
MN.no.damage <- MN.no.damage %>%  
  mutate(ID = group_indices(
    ., siteID, collectDate)) %>%
  arrange(ID)


# make a key for which ID goes with whoch siteID:collectDate combo
ID_key <- MN.no.damage %>%
  select(ID, siteID, collectDate) %>%
  distinct()

# add sample year, month, and event (1, 2, 3) to ID_key
# most sites have first sample event before June, and last in August or later
ID_key <- ID_key %>%
  mutate(year = year(collectDate),
         month = month(collectDate),
         sampleEvent = case_when(month < 6 ~ "first",
                                 month < 8 ~ "second",
                                 month >=8 ~ "third"))

# need to adjust some sample events.
# i.e. June at OKSR (Alaska) is first sample, not second
# COMO (Colorado) have second and third sample in early and late august, respectively
# MCRA 2017-8-01 is second sample (usually collected in July, but this year on august 1st)
ID_key[c(40, 42, 43, 45, 124, 132, 135), "sampleEvent"] <- 
  c("second", "first", "second", "second", "second", "first", "first")



# save objects for analysis in script 2 and 3
saveRDS(MN.no.damage, "data/macro_dw.RDS")
saveRDS(ID_key, "data/sample_ID_key.RDS")


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


