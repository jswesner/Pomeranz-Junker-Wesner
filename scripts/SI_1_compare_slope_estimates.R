#compare NAS,LBNBiom, Woodward

library(tidyverse)
library(sizeSpectra)
source("scripts/custom_functions.R")
library(lubridate)
library(cowplot)
# helper function for MLE method
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
       ylim = ylim_global, axes = FALSE, ...)
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
  if (panel == "b") {
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


# function to compare slopes####
compare_slopes <- function(data, dw_range, rsp_var, ...){
  # define breaks widths and mid bins for log2 and 6 equal log bins
  # log2
  breaks2 <- 2^(floor(log2(min(dw_range))):
                  ceiling(log2(max(dw_range))) )
  mid_bin_2 <- log10(breaks2[floor(length(breaks2)/2)]) 
  # 6 equal log bins
  breaks_log_6 <- exp(seq(floor(log(min(dw_range))),
                          ceiling(log(max(dw_range))),
                          length.out = 7))
  mid_bin_log_6 <- log10(breaks_log_6[floor(length(breaks_log_6)/2)])
  
  # Normalized Abundance Spectra
  NAS <- bin_and_center(data, var = rsp_var, breaks = breaks2)
  NAS$log_mids_center <- NAS$log_mids - mid_bin_2
  NAS_coefs <- coef(lm(log_count_corrected~log_mids_center,
                       data = NAS))
  # non-normalized abundance
  AS_coefs<- coef(lm(log_count~log_mids_center,
                     data = NAS))
  # equal logarithmic binning method
  
  # breaks_log_6 <- exp(seq(floor(log(min(dw_range))),
  #                         ceiling(log(max(dw_range))),
  #                         length.out = 7))
  # mid_bin_log_6 <- log10(breaks_log_6[floor(length(breaks_log_6)/2)])
  
  perkins <- bin_and_center(data, var = rsp_var, breaks = breaks_log_6)
  perkins$log_mids_center <- perkins$log_mids - mid_bin_log_6
  p_coef <- coef(lm(log_count~log_mids_center, data = perkins))
  pn_coef <- coef(lm(log_count_corrected~log_mids_center, data = perkins))
  
  # MLE method from Edwards et al. 2018
  mle_b <- MLE_tidy(data, rsp_var)$b
  
  slopes = data.frame(NAS = NAS_coefs[2],
                      AS = AS_coefs[2],
                      Perkins = p_coef[2],
                      PN = pn_coef[2],
                      mle = mle_b)
  slopes
}

# compare_slopes(data = d, dw_range = range(d$dw), rsp_var = "dw")



# # latitude and mean annual temperature provided in NEON site descriptions
site.info <- read.csv("data/aquatic-field-sites.csv")[,c(2,10)]
names(site.info) <- c("siteID","mat.c")

# read in estimated dry weight data
dw <- readRDS("data/macro_dw.RDS") 

# change data structure
dw <- dw %>%
  # filter small dw and NA values
  filter(!is.na(dw), dw >= 0.0026) %>%
  # calculate total count / per m2 for each body size 
  group_by(siteID, collectDate, dw) %>%
  summarise(count = sum(no_m2, na.rm = TRUE)) %>%
  ungroup() %>%
  # duplicate number of rows based on count
  #i.e. if count = 10, duplicate that row 10 times
  uncount(count)


# # fix sample dates in HOPB and MAYF ####
# 
# # i.e. HOPB has 2 samples from same date (2017-4-12) but different times
# # MAYF was collected on 2 dates in 2017-07, combining to one sample
# dw %>%
#   filter(siteID == "HOPB" | siteID == "MAYF") %>%
#   select(siteID, collectDate) %>%
#   unique() %>%
#   arrange(siteID, collectDate)
# 
# # make unique group index for siteID:collectDate combo
# dw <- dw %>%  
#   mutate(ID = group_indices(
#     ., siteID, collectDate))
# 
# # the following should result in 160 rows
# dw %>% select(ID) %>% unique() %>% nrow() 



# compare methods ####
dw_range = range(dw$dw, na.rm = TRUE)

test_method <- dw %>%
    # mutate(date = as.Date(collectDate)) %>%
    group_by(siteID, collectDate) %>%
    #filter(siteID == "CARI") %>%
    # create list-column
    nest() %>% 
    mutate(method_compare = 
             map(data,
                 compare_slopes,
                 rsp_var = "dw",
                 dw_range = dw_range)) %>%
    ungroup() %>%
    select(-data) %>%
    unnest(cols = method_compare)

#saveRDS(test_method, "results/slope_estimates_methods_SI.RDS")
#test_method <- readRDS("results/slope_estimates_methods_SI.RDS")
comp <- readRDS("results/perkins-compare.RDS")
test_method <- left_join(test_method, comp)
test_method <- left_join(test_method, site.info)

# plot all slope estimates and add regression line by method
test_method %>%
  mutate(MLBExponent = MLBExponent - 1) %>%
  pivot_longer(cols = 3:8,
               names_to = "method",
               values_to = "value") %>%
  left_join(site.info) %>%
  ggplot(aes(x = mat.c, 
             y = value, 
             color = method)) +
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE) +
  theme_bw()




MLEbins <- readRDS("data/MLEbins.RDS")

MLEbins <- MLEbins %>%
  select(siteID, b) %>%
  rename(MLEbins = b) %>%
  pivot_longer(2,
               names_to = "method",
               values_to = "value") %>%
  left_join(site.info)


slope_range <- test_method %>%
  pivot_longer(cols = 3:8,
               names_to = "method",
               values_to = "value") %>%
  bind_rows(MLEbins)

slope_q <- slope_range %>%
  group_by(method) %>%
  summarize(
    q_val = list(quantile(
      value,
      probs = c(0.025, 0.5, 0.975)))) %>%
  unnest(cols = q_val)

slope_q$q <- rep(c("lower", "median", "upper"), 7)

slope_q %>% pivot_wider(names_from = q, values_from = q_val) %>%
  arrange(median) %>%
  mutate(width = abs(lower - upper)) %>%
  write.csv("results/slope_est_quants.csv",
            row.names = FALSE)


(panelA <- test_method %>%
  rename(ELB = Perkins,
         MLE = mle,
         MLE_tail = MLBExponent) %>%
  pivot_longer(cols = 3:8,
               names_to = "method",
               values_to = "value") %>%
  bind_rows(MLEbins) %>%
  #left_join(site.info) %>%
  ggplot(aes(x = mat.c, 
             y = value, 
             color = method,
             linetype = method)) +
  geom_point(position = position_jitter(width = 0.3),
             size = 2, alpha = 0.4) +
  geom_line(stat = "smooth",
            method = "lm",
            size = 1.5, 
            alpha = 0.75) +
  theme_bw()+
  ylim(c(-2.8, 0.5)) +
  theme(legend.position = "none") +
  ylab("Estimated coefficient")
)

panelB <- test_method %>%
    mutate(MLBExponent = MLBExponent - 1) %>%
    rename(ELB = Perkins,
           NELB = PN,
           MLE = mle,
           MLE_tail = MLBExponent) %>%
    mutate(AS = AS - 1,
           ELB = ELB - 1) %>%
    pivot_longer(cols = 3:8,
                 names_to = "method",
                 values_to = "value") %>%
    bind_rows(MLEbins) %>%
    #left_join(site.info) %>%
    ggplot(aes(x = mat.c, 
               y = value, 
               color = method,
               linetype = method)) +
    geom_point(position = position_jitter(width = 0.3),
               size = 2, alpha = 0.4) +
    geom_line(stat = "smooth",
              method = "lm",
              size = 1.5, 
              alpha = 0.75) +
    # stat_smooth(method = "lm",
    #             se = FALSE) +
    theme_bw() +
    ylim(c(-2.8, 0.5)) +
    ylab("Corrected coefficient")

plot_grid(panelA, panelB, rel_widths = c(1.5, 2))

ggsave("plots/SI_compare_beta_mat.png")
