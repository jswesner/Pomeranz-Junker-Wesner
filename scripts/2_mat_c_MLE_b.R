# new analysis 10/8/20
# using mean annual temperature in Celsius as predictor variable

# load packages
# packages
library(sizeSpectra)
library(tidyverse)
library(brms)
library(janitor)
library(ggridges)

### Custom Functions ###
# the following function is modified from the eightMethodsMEE() function in the sizeSpectra package
# it is modified to only return the MLE estimate and 95 %CI for the bounded power law exponent, b, and is modified to work with list-columns in the tidyverse.  

# df is the data frame with every observation in it's own row
# rsp_var is the column to perform MLE on
# i.e. in this analysis the column is "dw" for estimated dry weight
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

# helper function to plot across lists of results
plot_b_est <- function(dat, b, grp_var, ...){
  isd_plot(sample(dat$dw, 10000),
           b = b$b,
           confVals = c(b$minCI, b$maxCI))
  # add labels
  mtext(paste0('Site = ',
               dat$siteID[1],
               as.character(dat[[grp_var]][1])),
        side = 3, line = -6, adj = 0.01)
  mtext(paste0("b = ", round(b, digits = 2)),
        side = 3, line = -7, adj = 0.05)
  # some of the plots have a text error in b-estimate, not sure why
}

# load and wrangle data--------------------------------------------------


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

# modify data structure
dw <- dw %>%
  filter(!is.na(dw)) %>%
  # calculate total count / per m2 for each body size 
  group_by(siteID, collectDate, dw) %>%
  summarise(count = sum(no_m2, na.rm = TRUE)) %>%
  ungroup() %>%
  # duplicate number of rows based on count
  #i.e. if count = 10, duplicate that row 10 times
  uncount(count)

# dw %>% group_by(siteID, collectDate) %>%
#   summarize(min_dw = min(dw),
#             max_dw = max(dw)) %>% View


# MLE ---------------------------------------------------------------------


# estimate MLE for each site and collection 
# note, this is running across 18+ million rows of data, it can take a while
mle.dw <- dw %>%
  #filter(!is.na(dw)) %>%
  mutate(date = as.Date(collectDate)) %>%
  group_by(siteID, collectDate) %>%
  # create list-column
  nest() %>% 
  # estimate b and 95% CI
  mutate(mle = map(data, MLE_tidy, rsp_var = "dw")) %>%
  unnest(cols = mle)

median(mle.dw$b)

b_filt <- dw %>%
  filter(dw >=0.0064) %>%
  mutate(date = as.Date(collectDate)) %>%
  group_by(siteID, collectDate) %>%
  # create list-column
  nest() %>% 
  # estimate b and 95% CI
  mutate(mle = map(data, MLE_tidy, rsp_var = "dw")) %>%
  unnest(cols = mle)
median(b_filt$b)

mle.dw <- b_filt
# separate year out of collectDate column
mle.dw <- mle.dw %>%
  separate(collectDate,
           c("year", "month", "day",
             "hour", "min", "sec")) %>%
  select(-hour, -min, -sec)

# join b estimates with site.info
# this should join by = c("siteID", "year")
mle.dw <- inner_join(mle.dw, site.info)

saveRDS(mle.dw, "data/mle_b_data.RDS")

# Bayesian Model ----------------------------------------------------------

b.mod <- brm(data = mle.dw,
               b ~ mat.c + (1 |siteID) + (1|year), 
               family = gaussian(),
               prior = c(prior(normal(0,0.1),
                               class = "b"),
                         prior(normal(-1.5,1),
                               class = "Intercept"),
                         prior(exponential(2),
                               class = "sigma"),
                         prior(exponential(2),
                               class = "sd")),
               control = list(adapt_delta = 0.999),
               chains = 4, 
               sample_prior = TRUE,
               iter = 6000,
               cores = 4)

# save model
saveRDS(b.mod, "data/mle_b_mod.RDS")

# model outputs and summaries
b.mod
fixef(b.mod)
plot(conditional_effects(b.mod), points = TRUE)

# posterior predictive check
# Supplemental figure S4 (saved in script 5)
pp_check(b.mod, type = "boxplot")


# slope positive or negative?
post.mod <- posterior_samples(b.mod)
# probability that slope is less than 0
sum(post.mod$b_mat.c < 0)/ nrow(post.mod)
# probability that slope is positive 
sum(post.mod$b_mat.c > 0)/ nrow(post.mod)

# main figure -------------------------------------------------------------

# main figure ####
# data to condition on
newdat <- tibble(
  mat.c = seq(min(mle.dw$mat.c),
                    max(mle.dw$mat.c),
                    length.out = 100)) 

# posterior samples to estimate from
posts <- posterior_samples(b.mod) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  sample_n(1000) %>% 
  expand_grid(newdat) %>% 
  mutate(b = b_intercept + b_mat_c*mat.c) %>% 
  group_by(mat.c) %>% 
  summarize(median = median(b),
            lower = quantile(b, probs = 0.025),
            upper = quantile(b, probs = 0.975))

# plot
# remove legend
(plot_b <- posts %>% 
    ggplot(aes(x = mat.c, y = median)) +
    geom_point(data = mle.dw,
               aes(y = b),
               size = 1.7) +
    geom_ribbon(alpha = 0.2,
                aes(ymax = upper,
                    ymin = lower)) + 
    geom_line() +
    labs(y = "b",
         x = "Mean annual temperature C") +
    # scale_color_brewer(type = "qual") +
    theme_classic() +
    # ylim(-5,5) +
    geom_text(aes(x = -Inf, y = Inf, label = "a)"),
              hjust = -.2,
              vjust = 1.2)  +
    NULL)

# save panel for use in script 4 to make the plots
saveRDS(plot_b, file = "plots/mat_c_b.rds")

# posterior distributions of slope estimate -------------------------------

site_b_est <- posterior_samples(b.mod) %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  #select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"), names_to = "site_ranef", values_to = "site_offset") %>%
  pivot_longer(contains("r_year"), names_to = "year_ranef", values_to = "year_offset") %>%
  mutate(siteID = str_sub(site_ranef, 11,14),
         siteID = toupper(siteID),
         year = str_sub(year_ranef, 8, 11)) %>%
  inner_join(site.info) %>%
  mutate(b_est = b_intercept + site_offset + year_offset +
           b_mat_c * mat.c) %>%
  group_by(siteID) %>%
  mutate(median_b = median(b_est, na.rm = TRUE)) %>%
  ungroup()

# global median
site_b_est %>% 
  summarize(
    med_b = median(b_est, na.rm = TRUE),
    l95 = quantile(b_est, probs = 0.025, na.rm = TRUE),
    u95 = quantile(b_est, probs = 0.975, na.rm = TRUE)) %>%
  arrange(med_b)

# smallest and largest slope
site_b_est %>% 
  group_by(siteID, mat.c) %>%
  summarize(
    med_b = median(b_est, na.rm = TRUE),
    l95 = quantile(b_est, probs = 0.025, na.rm = TRUE),
    u95 = quantile(b_est, probs = 0.975, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(med_b) %>%
  slice(c(1, n()))


(slope_post_A_panel <- site_b_est %>%
    filter(iter <= 2000) %>%
    ggplot() +
    geom_density_ridges_gradient(aes(
      x = b_est,
      y = reorder(siteID, mat.c),
      group = interaction(siteID, mat.c),
      fill = mat.c
      ),
      scale = 2,
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2) +
    scale_fill_viridis_c(alpha = 0.5, option = "plasma") +
    theme_bw() +
    #xlim(c(-1.5, -1)) +
    labs(y = "Site",
         x = "Slope estimate") +
    geom_text(aes(x = -Inf, y = Inf, label = "a)"),
              hjust = -.2,
              vjust = 1.2)  +
    NULL)

# save panel for use in script 4 to make the plots
saveRDS(slope_post_A_panel,
        file = "plots/slope_post_panel_A.RDS")

