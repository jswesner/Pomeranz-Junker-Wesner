# plots for manuscript

# all of the necessary data objects are saved in this R project. If making them from scratch, you need to first run scripts 2 and 3 in order to save the individual plot panels. This script just puts them together and modifies them for publication. 

# libraries
library(tidyverse)
library(cowplot)
library(maps)
library(ggthemes)
library(brms)
library(viridis)


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
                             color = "white") +
        geom_polygon(data = states,
                     aes(x = long,
                         y = lat,
                         group =group),
                     color = "white") +
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
        theme_map()
)
ggsave(map,
       file = "plots/map.png",
       width = 10,
       height = 8,
       dpi = 600)


# Figure 2 main text ------------------------------------------------------

#plot main figure
c2 <- readRDS("plots/b_mat_c.RDS")
d2 <- readRDS("plots/log_mg_mat_c.RDS")


main_plot <- plot_grid(
        c2 +
        theme_bw() +
        theme(legend.position = "none") +
        labs(y = "ISD exponent",
             x = "Standardized Mean Annual Temperature"),
        d2 +
        theme_bw() +
        labs(y = expression(
                "Log10 Dry Mass mg/"~m^2),
             x = "Standardized Mean Annual Temperature"),
        ncol = 2,
        align = "v",
        rel_widths = c(1, 1.3))

main_plot

ggsave(main_plot,
       file = "plots/main_plot.jpg",
       dpi = 600,
       width = 7,
       height = 5)


# Figure 3 main text ------------------------------------------------------

# plot posterior distributions
e3 <- readRDS("plots/b_mat_c_post_dist.RDS")
f3 <- readRDS("plots/log_mg_mat_c_post_dist.RDS")

post_plot <- plot_grid(
        e3 + 
                theme_bw() +
                theme(legend.position="none") +
                labs(x ="ISD exponent"),
        f3 + 
                theme_bw() + 
                guides(fill = guide_legend(
                        title = "Std(Temp.)")) +
        labs(x = expression(
                "Macroinvertebrate Dry Mass mg/"~m^2)),
        ncol = 2,
        align = "h",
        axis = "b",
        rel_widths = c(1, 1.3))
post_plot

ggsave(post_plot,
       file = "plots/post_plot.jpg",
       dpi = 600,
       width = 7,
       height = 5)



