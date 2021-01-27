# MLE bins example from edwards
library(sizeSpectra)
library(tibble)
library(tidyverse)

est_dw_range <- function(x){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001" with LW coefficients added using the LW_coeff() function, and the pertinent data columns from the NEON "fieldData" table (data product DP1.20120.001) already added
  
  # function uses tidyverse verbs and functions for ease of programming
  require(tidyverse)
  
  # caluclate bin windths
  # sizeClass = length of individual in mm
  # NEON data are measured to the nearest 1 mm. Therefore, the minimum possible measurement is sizeclass -0.5 and the maximum is size class + 0.49
  # e.g., if sizeClass == 5 mm, the minimum would be 5 - 0.5 = 4.5 and the maximum would be 5 + 0.49 = 5.49, and min/max for sizeclass == 6 is 5.5 and 6.49. 
  # note, for sizeclass == 1 we assume that the bin width is 0.5 to 1.49 (i.e. anything smaller than 0.5 is unlikely to be seen and measured). This assumes that all bin widths are the same size. 
  x <- x %>% 
    mutate(
      LngtMin = sizeClass - 0.5,
      LngtMax = sizeClass + 0.49,
      # calculate dw based on different formula types
      dw_bin = # this is the dry weight estimate using the measured sizeClass
        case_when(
          formula_type == 1 ~ a * sizeClass^b,
          formula_type == 2 ~ exp(a + b * log(sizeClass))),
      wmin = # this is the dry weight estimate using the minimum edge of the measured sizeClass bin
        case_when(
          formula_type == 1 ~ a * LngtMin^b,
          formula_type == 2 ~ exp(a + b * log(LngtMin))),
      wmax = # this is the dry weight estimate using the maximum edge of the measured sizeClass bin
        case_when(
          formula_type == 1 ~ a * LngtMax^b,
          formula_type == 2 ~ exp(a + b * log(LngtMax))))
  return(x)
}

MLE_tidy <- function(df, rsp_var){
  # define variables
  x <- df[[rsp_var]]
  xmin = min(x)
  xmax = max(x)
  log.x = log(x)
  sum.log.x = sum(log.x)
  
  # initial starting point for parameter estimate
  PL.bMLE = 1/(log(min(x)) - sum.log.x/length(x)) - 1
  
  # non-linear minimization  
  PLB.minLL = nlm(negLL.PLB, 
                  p = PL.bMLE, x = x, n = length(x), 
                  xmin = xmin, xmax = xmax,
                  sumlogx = sum.log.x)
  
  # estimate for b
  PLB.bMLE = PLB.minLL$estimate
  # minimum estimate of b
  PLB.minNegLL = PLB.minLL$minimum
  
  ## 95% CI calculation
  bvec = seq(PLB.bMLE - 1, PLB.bMLE + 1, 1e-05) # original =-0.5
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB(bvec[i],
                              x = x,
                              n = length(x), 
                              xmin = xmin,
                              xmax = xmax,
                              sumlogx = sum.log.x)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  # confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if (PLB.MLE.bConf[1] == min(bvec) | 
      PLB.MLE.bConf[2] == max(bvec)) {
    dev.new()
    plot(bvec, PLB.LLvals)
    abline(h = critVal, col = "red")
    stop("Need to make bvec larger - see R window")
  }
  # return b estimate and min/max 95% CI
  return(data.frame(b = PLB.bMLE,
                    minCI = min(bIn95),
                    maxCI = max(bIn95)))
}


# mean annual temperature provided in NEON site descriptions
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,10)]
names(site.info) <- c("siteID", "mat.c")


# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS") 
# remove lake and large river sites
lakes <- c("BARC", "CRAM", "LIRO", "PRLA", "PRPO",
           "SUGG", "TOOK", "BLWA", "FLNT", "TOMB")
dw <- dw %>%
  filter(!siteID %in% lakes)

dw$collectDate <- as.character(dw$collectDate)

# add 

# filter out one BLUE sample to troubleshoot
blue <- dw %>%
  filter(siteID == "BLUE",
         collectDate == "2018-03-12 18:05:00") 
blue <- blue %>%
  select(acceptedTaxonID, sizeClass, no_m2, formula_type, a, b) %>%
  rename(Number = no_m2) %>%
  mutate(SpecCode = as.factor(acceptedTaxonID))

blue <- blue %>% est_dw_range()

dataBin <-  blue %>%
  select(SpecCode, LngtMin, LngtMax, wmin, wmax, Number) %>%
  mutate(SpecCode = as.numeric(SpecCode))

#species_bins_plots(specCodeHighlight = NULL)

filt_min_size <- 0.1
dataBinForLike <- dataBin %>%
  filter(wmin >= filt_min_size)
n = sum(dataBinForLike$Number)
xmin = min(dataBinForLike$wmin)
xmax = max(dataBinForLike$wmax)

MLEbins = calcLike(
  negLL.fn = negLL.PLB.bins.species,
  p = -1.9,
  suppress.warnings = TRUE,
  dataBinForLike = dataBinForLike,
  n = n,
  xmin = xmin,
  xmax = xmax)

data.frame(#Year = fullYears[iii],
           xmin=xmin,
           xmax=xmax,
           n=n,
           b=MLEbins$MLE,
           confMin=MLEbins$conf[1],
           confMax=MLEbins$conf[2],
           min_size = filt_min_size)

blue %>%
  filter(dw_bin >= filt_min_size) %>%
  MLE_tidy("dw_bin")


# multiple samples per site -----------------------------------------------

arik <- dw %>% filter(siteID == "ARIK")

arik <-  arik %>%
  select(sizeClass, no_m2, formula_type,
         a, b, collectDate) %>%
  rename(Number = no_m2) 

arik <- arik %>% est_dw_range()

arik_samples <- sort(unique(arik$collectDate))

filt_min_size <- 0.0026
for( iii in 1:length(arik_samples)){
  dataBinForLike = 
    filter(arik, collectDate == arik_samples[iii])
    
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
    calcLike(negLL.fn = negLL.PLB.bins.species,
             p = -1.9,
             suppress.warnings = TRUE,
             dataBinForLike = dataBinForLike,
             n = n,
             xmin = xmin,
             xmax = xmax)
  
  if(iii == 1)
  {
    MLEbins = 
      data.frame(sample = arik_samples[iii],
                 xmin=xmin,
                 xmax=xmax,
                 n=n,
                 b=MLEbins.new$MLE,
                 confMin=MLEbins.new$conf[1],
                 confMax=MLEbins.new$conf[2])
  } else {
    MLEbins =
      rbind(MLEbins,
            c(arik_samples[iii],
              xmin,
              xmax,
              n,
              MLEbins.new$MLE,
              MLEbins.new$conf[1],
              MLEbins.new$conf[2]))
  }
}

range(MLEbins$b)
arik_mle_bins <- MLEbins

