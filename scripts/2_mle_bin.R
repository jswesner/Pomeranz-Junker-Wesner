# mle_bin to full data

# packages
library(sizeSpectra)
library(tidyverse)
library(brms)
library(janitor)
library(ggridges)
library(viridis)
library(lubridate)

source("scripts/custom_functions.R")

# process downloaded NEON data to include wmin and wmax for the MLEbin method presented in Edwards et al. 2020


# load and wrangle data--------------------------------------------------


# latitude and mean annual temperature provided in NEON site descriptions
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,6,10)]
names(site.info) <- c("siteID","latitude", "mat.c")


# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS") 
ID_key <- readRDS("data/sample_ID_key.RDS")


# estimate the range of estimated dw's based on width of sizeclass
dw_range <- est_dw_range(dw)

# select and rename columns for MLEbins method
dw_range <- dw_range %>%
  select(siteID, collectDate, ID,
         no_m2, dw_min, dw_max) %>%
  rename(Number = no_m2,
         wmin = dw_min, 
         wmax = dw_max)


# MLE ---------------------------------------------------------------------



filt_min_size <- 0.0026
collections <- sort(unique(dw_range$ID))

for( iii in 1:length(collections)){
  dataBinForLike = 
    filter(dw_range,
           ID == collections[iii])
  
  dataBinForLike = select(dataBinForLike,
                          wmin,
                          wmax,
                          Number)
  dataBinForLike = dataBinForLike %>%
    filter(wmin >= filt_min_size)
  
  n = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  MLEbins.new  = 
  try({
    calcLike(negLL.fn = negLL.PLB.bins.species,
             p = -1.9,
             suppress.warnings = TRUE,
             dataBinForLike = dataBinForLike,
             n = n,
             xmin = xmin,
             xmax = xmax)
    }
    )
  
  if(class(MLEbins.new) == "try-error"){
    MLEbins.new = list(b = NA, conf = c(NA, NA))
  }
  
  if(iii == 1)
  {
    MLEbins = 
      data.frame(ID = collections[iii],
                 xmin=xmin,
                 xmax=xmax,
                 n=n,
                 b=MLEbins.new$MLE,
                 confMin=MLEbins.new$conf[1],
                 confMax=MLEbins.new$conf[2])
  } else {
    MLEbins =
      rbind(MLEbins,
            c(ID = collections[iii],
              xmin,
              xmax,
              n,
              MLEbins.new$MLE,
              MLEbins.new$conf[1],
              MLEbins.new$conf[2]))
  }
}


MLEbins <- left_join(MLEbins, ID_key, by = "ID")
MLEbins <- left_join(MLEbins, site.info,
                     by = "siteID")
saveRDS(MLEbins, "data/MLEbins.RDS")


