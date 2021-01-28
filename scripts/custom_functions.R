# custom functions for analysis 
# written by Justin Pomeranz
# jfpomeranz@gmail.com
# started June 2020

# function to add published length weight regression coefficients to macroinvertebrate taxonomic data from NEON sites. Only returns observations which has a corresponding LW coefficient. i.e. output may be smaller than input
LW_coef <- function(x, lw_coef, percent = FALSE){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001"
  # lw_coef = table of taxon specific length_weight coefficients including higher taxonomic classification. Note that this table must have a column named "taxon" which is used to match up with the different taxonomic levels in x
  # percent = Report percent of rows in x which have a taxon-specific LW equation? default is FALSE
  
  
  # record number of rows in x to calculate percent of observations with a LW equation
  x.nrow <- nrow(x)
  # make indexes for each taxonomic level and extract those rows
  # remove rows from x that are accounted for in taxononmic tables
  
  # genus
  genus_index <- which(x$genus %in% lw_coef$taxon)
  # make table of observations which have a genus-specific LW equation
  genus_table <- x[genus_index,]
  # remove these from global taxonomic data
  x <- x[-genus_index,]
  
  # repeat process for other taxonomic levels
  
  # family
  family_index <- which(x$family %in% lw_coef$taxon)
  family_table <- x[family_index,]
  x <- x[-family_index,]
  
  # order
  order_index <- which(x$order %in% lw_coef$taxon)
  order_table <- x[order_index,]
  x <- x[-order_index,]
  
  # subclass
  subclass_index <- which(x$subclass %in% lw_coef$taxon)
  subclass_table <- x[subclass_index,]
  x <- x[-subclass_index,]
  
  # class
  class_index <- which(x$class %in% lw_coef$taxon)
  class_table <- x[class_index,]
  x <- x[-class_index,]
  
  # combine coefficients to taxonomic tables
  lw_cols <- c("taxon",
               "a",
               "b",
               "L_units",
               "dw_units",
               "formula_type",
               "formula")
  genus_table <- merge(genus_table, lw_coef[,lw_cols], 
                       by.x = "genus",
                       by.y = "taxon", all.x = TRUE)
  family_table <- merge(family_table, lw_coef[,lw_cols], 
                        by.x = "family",
                        by.y = "taxon", all.x = TRUE)
  order_table <- merge(order_table, lw_coef[,lw_cols], 
                       by.x = "order",
                       by.y = "taxon", all.x = TRUE)
  subclass_table <- merge(subclass_table, lw_coef[,lw_cols], 
                          by.x = "subclass",
                          by.y = "taxon", all.x = TRUE)
  class_table <- merge(class_table, lw_coef[,lw_cols], 
                       by.x = "class",
                       by.y = "taxon", all.x = TRUE)
  datout <- rbind(genus_table, family_table, order_table, subclass_table, class_table)
  
  if (percent == TRUE){
    message(paste(nrow(datout)/x.nrow * 100,
                "percent of input data had Length-Weight equations"))
  }
  return(datout)
}

# estimate dry weight (dw) values from length-weight regression coefficients
est_dw <- function(x, fieldData){
  # x = inv_taxonomyProcessed table from NEON data product "DP1.20120.001" with LW coefficients added using the LW_coeff() function
  # fieldData = inv_fieldData table from NEON data product "DP1.20120.001"
  
  # functions uses tidyverse functions
  require(tidyverse)
  
  # simplify fieldData to three columns
  field = fieldData %>%
    select(siteID, sampleID, benthicArea)
  # add benthicArea column from fieldData to x. This is necessary to calculate number per m2 below
  
  # join by siteID and sampleID
  x <- left_join(x, field, by = c("siteID", "sampleID"))
  
  # correct counts to per meter squared. Round to nearest integer
  x <- x %>%
    mutate(no_m2 = round(estimatedTotalCount / benthicArea))
  
  # calculate dw based on different formula types
  # sizeClass = length of individual in mm
  x <- x %>% 
    mutate(dw = case_when(formula_type == 1 ~ a * sizeClass^b,
                          formula_type == 2 ~ exp(a + b * log(sizeClass))))
  return(x)
}

# function modified from Pomeranz, Warburton, and Harding. Anthropogenic mining alters macroinvertebrate size spectra in streams. Freshwater Biology. 2019. 
bin_and_center <- function(data, var, breaks, ...){
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist(). Can be a vector giving breakpoints of bins [i.e. Log10 breaks = 10^seq(min, max)], a single number designating number of bins, etc. See ?hist for details. If not supplied by the user a default of log2 width bins are calculated 
  
  # are breaks supplied by the user?
  if(exists("breaks") == FALSE){
    
    # calculate log2 width bins which are inclusive of the range of data supplied
    breaks = 2^seq(floor(range(log2(data[[var]]))[1]),
                   ceiling(range(log2(data[[var]]))[2]))
    message("breaks not supplied, using log2 width bins")
  }
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks,
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  bin_width = breaks_offset - breaks_orig
  count = binned_hist$counts
  log_mids = log10(binned_hist$mids)
  biomass = count * 10**log_mids
  nbiomass = log10(biomass / bin_width)
  dataout = data.frame(
    count = count,
    log_count = log10(count),
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / bin_width),
    # original midpoint of bin log10 transformed
    log_mids = log_mids,
    bin_width = bin_width,
    biomass = biomass,
    nbiomass = nbiomass)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # # recenter data at x=0
  # mid_row = ceiling(nrow(dataout)/2)
  # # subtract value of mid row from all mids
  # dataout$log_mids_center = 
  #   dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}

deg_days <- function(x){
  # x = TSW_30min object from NEON surface water temperature data product dpID = "DP1.20053.001"
  
  require(tidyverse)
  # separate startDateTime into components
  x <- x %>%
    as_tibble() %>%
    # first, need to filter out missing data
    filter(!is.na(surfWaterTempMean)) %>%
    # select only columns that we will need below
    select(siteID, startDateTime, surfWaterTempMean) %>%
    separate(startDateTime,
             c("year", "month", "day", "hour", "min", "sec")) %>%
    # group by individual days
    group_by(siteID, year, month, day) %>%
    # caluculate daily mean temperature
    summarize(daily = mean(surfWaterTempMean)) %>%
    ungroup() %>% # ungroup for next step
    # only want the positive degree days
    # i.e. a negative temperature day doesn't "remove" degree days
    filter(daily >=0) %>%
      # group by site ID and year
    group_by(siteID, year) %>%
      # calculate cumulutaive sum of mean daily temperatures
    mutate(degree_days = cumsum(daily))
    return(x)
}

# plot prior and posterior distributions
# this function is currently only for Intercept and b prioirs
plot_prior_vs_posterior <- function(model,
                                    prior_intercept_mu,
                                    prior_intercept_sigma,
                                    prior_b_mu,
                                    prior_b_sigma){
  require(dplyr)
  message("this function is currently programed for 3 beta coefficients")
  # int_coef <- fixef(model)[1,1:2]
  # slope_coef <- fixef(model)[2,1:2]
  
  # get distributions of coefficients
  fixed_effects <- row.names(fixef(model))
  prior.list = NULL
  for (i in seq_along(fixed_effects)){
    prior.list[[i]] =
      data.frame(value =
                   rnorm(1000, fixef(model)[i,1], fixef(model)[i,2]),
                 sample = fixed_effects[i])
  }
  # sample prior distributions
  post.dat <- data.frame(
       value = c(rnorm(1000, prior_intercept_mu, prior_intercept_sigma),
                 rnorm(1000, prior_b_mu, prior_b_sigma)),
       sample = rep(c("Intercept prior", "b priors"), each = 1000))
  # combine prior and posts for plot
  plot.dat <- bind_rows(prior.list, post.dat)
  # plot
  ggplot(plot.dat, aes(x = value, fill = sample))+
    geom_density(alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = c("mediumorchid4", "mediumorchid1",
                                 "turquoise4", "steelblue1", "cadetblue1",
                                 "skyblue1"),
                      breaks = c("Intercept prior", "Intercept",
                                 "b priors", "log_mids_center",
                                 "log10degree_days",
                                 "log_mids_center:log10degree_days"))
}


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
      sizeClassMin = sizeClass - 0.5,
      sizeClassMax = sizeClass + 0.49,
      # calculate dw based on different formula types
      dw_bin = # this is the dry weight estimate using the measured sizeClass
        case_when(
          formula_type == 1 ~ a * sizeClass^b,
          formula_type == 2 ~ exp(a + b * log(sizeClass))),
      dw_min = # this is the dry weight estimate using the minimum edge of the measured sizeClass bin
        case_when(
          formula_type == 1 ~ a * sizeClassMin^b,
          formula_type == 2 ~ exp(a + b * log(sizeClassMin))),
      dw_max = # this is the dry weight estimate using the minimum edge of the measured sizeClass bin
        case_when(
          formula_type == 1 ~ a * sizeClassMax^b,
          formula_type == 2 ~ exp(a + b * log(sizeClassMax))))
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

MLE_bin_tidy <- function(df){
  test1 <- all( c("wmin", "wmax", "Number") %in% names(df))
  if(test1 == FALSE){
    stop("DF needs to have 'wmin', 'wmax' and 'Number' columns")
  }
  n = sum(df$Number)
  xmin = min(df$wmin)
  xmax = max(df$wmax)
  
  like  <- calcLike(negLL.fn = negLL.PLB.bins.species,
                    p = -1.9,
                    suppress.warnings = TRUE,
                    dataBinForLike = df,
                    n = n,
                    xmin = xmin,
                    xmax = xmax)
  
  return(data.frame(b = like$MLE,
                    bL95 = like$conf[1],
                    bU95 = like$conf[2]))
}

# modified function to plot MLE estimate of size spectra (ISD)
isd_plot <- function (x, b, confVals = NULL,
                      panel = "b", log.xy = "xy",
                      mgpVals = c(1.6, 0.5, 0),
                      inset = c(0, -0.04),
                      xlim_global = NA,
                      ylim_global = NA, ...) 
{
  if (is.na(xlim_global[1])) {
    xlim_global = c(min(x), max(x))
  }
  if (is.na(ylim_global[1])) {
    ylim_global = c(1, length(x))
  }
  plot(sort(x, decreasing = TRUE), 1:length(x), log = log.xy, 
       xlab = expression(paste("Values, ", italic(x))), 
       ylab = expression(
         paste("Number of ", values >= x), sep = ""),
       mgp = mgpVals, xlim = xlim_global, 
       ylim = ylim_global, axes = FALSE)#, ...)
  xLim = 10^par("usr")[1:2]
  yLim = 10^par("usr")[3:4]
  if (log.xy == "xy") {
    logTicks(xLim, yLim, xLabelSmall = c(5, 50, 500))
  }
  if (log.xy == "x") {
    mgpVal = c(2, 0.5, 0)
    logTicks(xLim, yLim = NULL, xLabelSmall = c(5, 50, 500), 
             mgpVal = mgpVal)
    yBig = c(0, 500, 1000)
    axis(2, at = yBig, labels = yBig, mgp = mgpVal)
    axis(2, seq(yBig[1], yBig[length(yBig)], by = 100),
         labels = rep("", 11), tcl = -0.2, mgp = mgpVal)
  }
  x.PLB = seq(min(x), max(x), length = 1000)
  y.PLB = (1 - pPLB(x = x.PLB,
                    b = b,
                    xmin = min(x.PLB),
                    xmax = max(x.PLB))) * 
    length(x)
  lines(x.PLB, y.PLB, col = "red")
  if (panel == "b" & !is.null(confVals)) {
    for (i in c(1, length(confVals))) {
      lines(x.PLB,
            (1 - pPLB(x = x.PLB, b = confVals[i],
                      xmin = min(x.PLB),
                      xmax = max(x.PLB))) * length(x), 
            col = "red", lty = 2)
    }
    #legend("topright", "(b)", bty = "n", inset = inset)
  }
  if (panel == "h") {
    legJust(c("(h) MLE",
              paste("b=", signif(b, 3), sep = "")),
            inset = inset, logxy = TRUE)
  }
}

# helper function to plot across lists of results
plot_b_est <- function(dat, b, ...){
  isd_plot(dat$dw,
           b = b$b,
           confVals = c(b$confMin, b$confMax))
  # add labels
  
  mtext(paste(dat$siteID[1]#, 
              #str_sub(as.character(
              #  dat$collectDate[1], 1, 7))
  ),
  side = 3, line = -6, adj = 0.01)
  mtext(paste0("b = ", round(b$b, digits = 2)),
        side = 3, line = -7, adj = 0.05)
  # some of the plots have a text error in b-estimate, not sure why
}