# isd plots

library(tidyverse)
library(sizeSpectra)
library(lubridate)

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
  
  mtext(paste(dat$siteID[1],
              b$year,
              "-",
              as.character(b$month)
              #, 
              #str_sub(as.character(
              #  dat$collectDate[1], 1, 7))
  ),
  side = 3, line = -6, adj = 0.01)
  mtext(paste0("b = ", round(b$b, digits = 2)),
        side = 3, line = -7, adj = 0.05)
  # some of the plots have a text error in b-estimate, not sure why
}

png_isd_plot <- function(dat, b, ...){
  # 1. Open png file
  png(paste0(
    "results/isd_plots_SI/isd_ID00", # added 00, haven't tested it yet
    dat$ID[1],
    ".png"),
    width = 600,
    height = 600)
  # 2. Create the plot
  plot_b_est(dat = dat, b = b)
  # 3. Close the file
  dev.off()
}


# ####
# read in data ####
MLEbins <- readRDS("data/MLEbins.RDS")
#MLEbins$month <- month(MLEbins$collectDate)

# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS") 
# remove lake and large river sites
lakes <- c("BARC", "CRAM", "LIRO", "PRLA", "PRPO",
           "SUGG", "TOOK", "BLWA", "FLNT", "TOMB")
dw <- dw %>%
  filter(!siteID %in% lakes)

dw <- dw %>%
  filter(!is.na(dw), !is.na(no_m2)) 

# # i.e. HOPB has 2 samples from same date (2017-4-12) but different times
# # MAYF was collected on 2 dates in 2017-07, combining to one sample
# dw %>%
#   filter(siteID == "HOPB" | siteID == "MAYF") %>%
#   select(siteID, collectDate) %>%
#   unique() %>%
#   arrange(siteID, collectDate)
# 
# # change HOPB 2017-04-12 13:23:00 --> 2017-04-12 16:00
# dw[
#   dw$siteID == "HOPB" &
#     dw$collectDate == 
#     as.POSIXct("2017-04-12 13:23:00", tz = "GMT"),
#   "collectDate"] <- as.POSIXct("2017-04-12 16:00:00", tz = "GMT")
# 
# # change MAYF 2017-07-18 14:51:00 to 2017-07-20 14:37:00
# dw[
#   dw$siteID == "MAYF" &
#     dw$collectDate == 
#     as.POSIXct("2017-07-18 14:51:00", tz = "GMT"),
#   "collectDate"] <- as.POSIXct("2017-07-20 14:37:00", tz = "GMT")

# modify data structure
dw <- dw %>%
  filter(!is.na(dw), dw>=0.0026) %>%
  # calculate total count / per m2 for each body size 
  group_by(siteID, collectDate, dw) %>%
  summarise(count = sum(no_m2, na.rm = TRUE)) %>%
  ungroup() %>%
  # duplicate number of rows based on count
  #i.e. if count = 10, duplicate that row 10 times
  uncount(count)

dw_plot_data <- dw %>% 
  select(siteID, collectDate, dw) %>%
  mutate(ID = group_indices(
    ., siteID, collectDate)) %>%
    #,
    # date_id = str_sub(
    #   as.character(collectDate)),
    # 1, 7) %>%
  arrange(ID) %>%
  select(-collectDate)


# plots --------------------------------------------------

# split objects into equal sized lists
plot_data_list <- split(dw_plot_data, dw_plot_data$ID)
MLE_b_list <- split(MLEbins, MLEbins$ID)

# save all isd plots
walk2(plot_data_list, MLE_b_list, png_isd_plot)
