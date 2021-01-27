# mle_bin to full data

# packages
library(sizeSpectra)
library(tidyverse)
library(brms)
library(janitor)
library(ggridges)
library(viridis)
library(lubridate)

# process downloaded NEON data to include wmin and wmax for the MLEbin method presented in Edwards et al. 2020

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

# range(MLEbins$b)
# median(MLEbins$b)
# MLEbins %>% select(ID, b)

MLEbins <- left_join(MLEbins, ID_key, by = "ID")
MLEbins <- left_join(MLEbins, site.info,
                     by = "siteID")

# ggplot(MLEbins,
#        aes(x = latitude,
#            y = b,
#            ymin = confMin,
#            ymax = confMax,
#            color = reorder(siteID, latitude))) +
#   geom_pointrange() +
#   stat_smooth(method = "lm", inherit.aes = FALSE,
#               aes(latitude, b)) +
#   scale_color_viridis_d(option = "plasma", direction = -1) +
#   theme_bw()

saveRDS(MLEbins, "data/MLEbins.RDS")
#MLEbins <- readRDS("data/MLEbins.RDS")


# Bayesian regression -----------------------------------------------------

### incorporate variation, confMin/Max into b estimate?
b.mod <- brm(data = MLEbins,
             b ~ latitude + (1 |siteID) + (1|year), 
             family = gaussian(),
             prior = c(prior(normal(0,0.1),
                             class = "b"),
                       prior(normal(-1.5,2),
                             class = "Intercept"),
                       prior(exponential(2),
                             class = "sigma"),
                       prior(exponential(2),
                             class = "sd")),
             control = list(
               adapt_delta = 0.999,
               max_treedepth = 11),
             chains = 4, 
             sample_prior = TRUE,
             iter = 6000,
             cores = 4)

# save model
saveRDS(b.mod, "results/MLEbins_latitude_brm.RDS")
# b.mod <- readRDS("results/MLEbins_latitude_brm.RDS")
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
sum(post.mod$b_latitude < 0)/ nrow(post.mod)
# probability that slope is positive 
sum(post.mod$b_latitude > 0)/ nrow(post.mod)


# main figure -------------------------------------------------------------

# main figure ####
# data to condition on
newdat <- tibble(
  latitude = seq(min(MLEbins$latitude),
              max(MLEbins$latitude),
              length.out = 100)) 

# posterior samples to estimate from
posts <- posterior_samples(b.mod) %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  sample_n(500) %>% 
  expand_grid(newdat) %>% 
  mutate(b = b_intercept + 
           b_latitude * latitude) %>% 
  group_by(latitude) %>% 
  summarize(median = median(b),
            lower = quantile(b, probs = 0.025),
            upper = quantile(b, probs = 0.975))
# plot
# remove legend
(plot_b <- posts %>% 
    ggplot(aes(x = latitude,
               y = median)) +
    geom_point(data = MLEbins,
               aes(y = b),
               size = 1.7) +
    geom_ribbon(alpha = 0.2,
                aes(ymax = upper,
                    ymin = lower)) + 
    geom_line() +
    labs(y = "b",
         x = "Latitude") +
    # scale_color_brewer(type = "qual") +
    theme_classic() +
    # ylim(-5,5) +
    geom_text(aes(x = -Inf, y = Inf, label = "a)"),
              hjust = -.2,
              vjust = 1.2)  +
    NULL)


# save panel for use in script 4 to make the plots
saveRDS(plot_b, file = "plots/latitude_b.rds")

# posterior distributions of slope estimate -------------------------------

site_b_est <- posterior_samples(b.mod) %>%
  clean_names() %>%
  mutate(iter = 1:nrow(.)) %>% 
  #select(!contains("cor_site")) %>% 
  pivot_longer(contains("r_site"),
               names_to = "site_ranef",
               values_to = "site_offset") %>%
  pivot_longer(contains("r_year"),
               names_to = "year_ranef",
               values_to = "year_offset") %>%
  mutate(siteID = str_sub(site_ranef, 11,14),
         siteID = toupper(siteID),
         year = str_sub(year_ranef, 8, 11)) %>%
  inner_join(site.info) %>%
  mutate(
    b_est = 
      b_intercept + site_offset + year_offset +
           b_latitude * latitude) %>%
  group_by(siteID) %>%
  mutate(median_b = median(b_est, na.rm = TRUE)) %>%
  ungroup()

# global median
site_b_est %>% 
  summarize(
    med_b = median(b_est, na.rm = TRUE),
    l95 = quantile(b_est, probs = 0.025,
                   na.rm = TRUE),
    u95 = quantile(b_est, probs = 0.975,
                   na.rm = TRUE)) %>%
  arrange(med_b)

# smallest and largest slope
site_b_est %>% 
  group_by(siteID, latitude) %>%
  summarize(
    med_b = median(b_est, na.rm = TRUE),
    l95 = quantile(b_est, probs = 0.025,
                   na.rm = TRUE),
    u95 = quantile(b_est, probs = 0.975,
                   na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(med_b) %>%
  slice(c(1, n()))

# panel A fig 3(?)
(slope_post_A_panel <- site_b_est %>%
    filter(iter <= 2000) %>%
    ggplot() +
    geom_density_ridges_gradient(aes(
      x = b_est,
      y = reorder(siteID, latitude),
      group = interaction(siteID, latitude),
      fill = latitude
    ),
    scale = 2,
    rel_min_height = 0.01,
    quantile_lines = TRUE,
    quantiles = 2) +
    scale_fill_viridis_c(alpha = 0.5,
                         option = "plasma",
                         direction = -1) +
    theme_bw() +
    #xlim(c(-1.5, -1)) +
    labs(y = "Site",
         x = "Slope estimate") +
    geom_text(aes(x = -Inf, y = Inf, label = "a)"),
              hjust = -.2,
              vjust = 1.2)  +
    NULL)

saveRDS(slope_post_A_panel, "plots/slope_post_panel_A.RDS")

# update models with month ------------------------------------------------

# # make month variable an ordered factor
# MLEbins <- MLEbins %>%
#   mutate(sampleEvent = factor(sampleEvent, ordered = TRUE,
#                               levels = c("first", "second", "third"))) 


mod2 <- brm(data = MLEbins,
                     b ~ latitude + sampleEvent + 
              (1 |siteID) + (1|year), 
                     family = gaussian(),
                     prior = c(prior(normal(0,0.25),
                                     class = "b"),
                               prior(normal(-1.5,2),
                                     class = "Intercept"),
                               prior(exponential(1),
                                     class = "sigma"),
                               prior(exponential(1),
                                     class = "sd")),
                     control = list(adapt_delta = 0.999),
                     chains = 4, 
                     sample_prior = TRUE,
                     iter = 6000,
                     cores = 4)
saveRDS(mod2, "results/MLEbins_latitude_season_brm.RDS")

mod3 <- brm(data = MLEbins,
            b ~ latitude * sampleEvent + (1 |siteID) +
              (1|year), 
            family = gaussian(),
            prior = c(prior(normal(0,0.25),
                            class = "b"),
                      prior(normal(-1.5,2),
                            class = "Intercept"),
                      prior(exponential(1),
                            class = "sigma"),
                      prior(exponential(1),
                            class = "sd")),
            control = list(adapt_delta = 0.999),
            chains = 4, 
            sample_prior = TRUE,
            iter = 6000,
            cores = 4)
saveRDS(mod3, "results/MLEbins_interaction_brm.RDS")

fixef(b.mod)
fixef(mod2)
fixef(mod3)


b.mod <- add_criterion(b.mod, "waic")
mod2 <- add_criterion(mod2, "waic")
mod3 <- add_criterion(mod3, "waic")
loo_compare(b.mod, mod2, mod3, criterion = "waic")

loomod <- loo(b.mod, reloo = TRUE)
loomod2 <- loo(mod2, reloo = TRUE)
loomod3 <- loo(mod3, reloo = TRUE)

# fixef(b.mod)
# fixef(mod2)
# fixef(mod3)
loo_model_weights(list(loomod, loomod2, loomod3))


plot(conditional_effects(b.mod), points = TRUE)
plot(conditional_effects(mod2), points = TRUE)
plot(conditional_effects(mod3), points = TRUE)


b.mod <- readRDS("results/MLEbins_latitude_brm.RDS")
b.mod2 <- readRDS("results/MLEbins_latitude_season_brm.RDS")
b.mod3 <- readRDS("results/MLEbins_interaction_brm.RDS")


library(future)
plan(multiprocess)
#plan(sequential)
b.mod1_k <- kfold(b.mod, k = 10, chains = 1)
b.mod2_k <- kfold(b.mod2, k = 10, chains = 1)
#plan(multiprocess)
b.mod3_k <- kfold(b.mod3, k = 10, chains = 1)

mod_kfold_point <- cbind(
  b.mod1_k$pointwise[,"elpd_kfold"],
  b.mod2_k$pointwise[,"elpd_kfold"],
  b.mod3_k$pointwise[,"elpd_kfold"]
)

stacking_weights(mod_kfold_point)
