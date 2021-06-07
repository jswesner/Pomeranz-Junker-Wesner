# plots for manuscript

# all of the necessary data objects are saved in this R project. If making them from scratch, you need to first run scripts 2 and 3 in order to save the individual plot panels. This script just puts them together and modifies them for publication. 

# libraries
library(tidyverse)
library(cowplot)
library(maps)
library(ggthemes)
library(brms)
library(viridis)
library(ggridges)
library(here)

# create plots directory, if it doesn't already exist.
# this should have already been run in script 4, but including here in case people run scripts out of order. 
if(!dir.exists("plots")){
        dir.create("plots")
}

# Figure 1 main text ------------------------------------------------------

# plot map of sites
# get coordinates for polygons - function in ggplot
world <- map_data("world")
states <- map_data("state")

# read in site info with lat long
site.info <- read.csv("data/aquatic-field-sites.csv")

#make the map
(map <- ggplot() +
        geom_polygon(data = world,
                     aes(x = long, y = lat, group = group),
                             color = "white",
                     fill = "grey") +
        geom_polygon(data = states,
                     aes(x = long,
                         y = lat,
                         group =group),
                     color = "white",
                     fill = "grey") +
        coord_quickmap(ylim = c(18,70),
                       xlim = c(-160,-50))+
        geom_point(data = site.info,
                   aes(x = Longitude,
                       y = Latitude,
                       fill = mat.c,
                       color = mat.c),
                   size = 4,
                   alpha = 0.9,
                   shape = 21) +
        scale_fill_viridis_c(option = "plasma") +
        scale_color_viridis_c(option = "plasma") +
                labs(color = expression("Mean Annual\nTemperature " ( degree*C))) +
                guides(fill = F) +
        theme_void()
)

ggsave(map,
       file = "plots/map.png",
       width = 10,
       height = 8,
       dpi = 600)


# Figure 2 main text ------------------------------------------------------


#plot main figure

c2 <- readRDS("plots/b_mat_c.RDS")
d2 <- readRDS("plots/log_mg_mat.RDS")


main_plot <- plot_grid(
        c2 ,
        d2 +guides(color = F) +
                labs(y = expression(
                        "Macroinvertebrate Dry Mass g/"~m^2)),
        ncol = 1,
        align = "vh",
        labels = "auto")


main_plot

ggsave(main_plot,
       file = "plots/main_plot.jpg",
       dpi = 600,
       width = 5,
       height = 7)


# Figure 3 main text ------------------------------------------------------

# plot posterior distributions at each site

#extract posts with custom function

e3 <- readRDS("plots/b_mat_c_post_dist.RDS")
f3 <- readRDS("plots/log_mg_mat_canopy_post_dist.RDS")


post_plot <- plot_grid(
        e3 + 
                theme_bw() +
                theme(legend.position="none") +
                labs(x ="ISD exponent"),
        f3 + 
                theme_bw() + 
        labs(x = expression(
                "Macroinvertebrate Dry Mass g/"~m^2)),
        ncol = 2,
        align = "h",
        axis = "b",
        rel_widths = c(1, 1.55),
        labels = "auto")
post_plot

ggsave(post_plot,
       file = "plots/post_plot.jpg",
       dpi = 600,
       width = 7,
       height = 5)






# Literature Figure Comparisons -------------------------------------------

# literature comparison

lit <- read_csv("data/temp_summaries_table.csv") %>% 
        filter(Include == "Y")

mod_best <- readRDS("results/b_mat_c_brms.RDS")

range_bmat <- fitted(mod_best, newdata = data.frame(mat.c = c(min(mod_best$data$mat.c),
                                                                  max(mod_best$data$mat.c))),
                                                        re_formula = NA, summary = F) %>% 
        as.data.frame() %>% rename(min = V1, max = V2) %>% 
        mutate(abs_diff = abs(min - max))


range_bmat_summary <- range_bmat %>% summarize(mean = median(abs_diff),
                                               lower = quantile(abs_diff, probs = 0.025),
                                               upper = quantile(abs_diff, probs = 0.975),
                                               mean_low = mean(min),
                                               mean_high = mean(max))



(lit_plot <- lit %>% 
                mutate(Driver = fct_relevel(Driver, "Temperature", "Land Use")) %>% 
                ggplot(aes(x = reorder(Author, -b_diff), y = b_diff)) +
                coord_flip() +
                # geom_segment(aes(y= 0, yend = b_diff, xend = reorder(Author, -b_diff))) +
                geom_pointrange(aes(ymin =b_diff - error, ymax = b_diff + error, shape = Driver, 
                                    fill = Driver), size = 1) +
                geom_hline(yintercept = range_bmat_summary$mean) +
                annotate("text", x = 12, y = 0.7, label = "This study (median and 95% CrI)") +
                geom_rect(aes(xmin = 0, xmax = 13, ymin = range_bmat_summary$lower, 
                              ymax = range_bmat_summary$upper), color = NA, alpha = 0.01) +
                # scale_fill_brewer(type = "qual") +
                # facet_grid(Driver ~ .) +
                scale_shape_manual(values = c(21, 22, 23, 24)) +
                scale_fill_manual(values = c("black", "white", "white", "white")) +
                labs(y = "Absolute change in ISD exponent (or slope)") +
                theme_classic() +
                theme(axis.title.y = element_blank(),
                      text = element_text(size = 15),
                      axis.text.y = element_text(size = 10)) +
                annotate("segment", x = 11.6, y = 0.43, xend = 11.6, yend = 0.24, arrow=arrow(type = "closed")) +
                ylim(0,1))


ggsave(lit_plot, file = "plots/lit_plot.jpg", width = 9, height = 3.5, units = "in", dpi = 500)




