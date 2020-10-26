# plots for manuscript

# all of the necessary data objects are saved in this R project. If making them from scratch, you need to first run scripts 2 and 3 in order to save the individual plot panels. This script just puts them together and modifies them for publication. 

# libraries
library(tidyverse)
library(cowplot)
library(maps)
library(ggthemes)
library(brms)


# Figure 1 main text ------------------------------------------------------

# plot map of sites
# get coordinates for polygons - function in ggplot
world <- map_data("world")
states <- map_data("state")

# read in site info with lat long
site.info <- read.csv("data/aquatic-field-sites.csv")

#make the map
map <- ggplot() + 
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
                       y = Latitude),
                   size = 3,
                   shape = 21,
                   fill = "yellow") +
        theme_map()

ggsave(map,
       file = "plots/map.png",
       width = 10,
       height = 8,
       dpi = 600)


# Figure 2 main text ------------------------------------------------------

#plot main figure
c <- readRDS("plots/mat_c_b.rds")
d <- readRDS("plots/mat_c_biomass.rds")


main_plot <- plot_grid(
        c +
        theme_bw() +
        labs(y = "Size spectrum exponent",
             x = expression("Mean annual temperature "( degree*C))),
        d +
        theme_bw() +
        labs(y = expression("Macroinvertebrate Dry Mass mg/"~m^2),
             x = expression("Mean annual temperature "( degree*C))),
        ncol = 2,
        align = "v")

main_plot

ggsave(main_plot,
       file = "plots/main_plot.jpg",
       dpi = 600,
       width = 7,
       height = 3)


# Figure 3 main text ------------------------------------------------------

# plot posterior distributions
e <- readRDS("plots/slope_post_panel_A.RDS")
f <- readRDS("plots/biomass_post_panel_B.rds")

post_plot <- plot_grid(
                e + 
                theme_bw() +
                theme(legend.position="none") +
                labs(x ="Size spectrum exponent"),
                f + 
                theme_bw() + 
                guides(fill = guide_legend(title = 
                        expression(atop("Temperature ",
                        paste(( degree*C)))))) +
                labs(x = expression("Macroinvertebrate Dry Mass mg/"~m^2)),
        ncol = 2,
        align = "h",
        axis = "b",
        rel_widths = c(1, 1.3))

ggsave(post_plot,
       file = "plots/post_plot.jpg",
       dpi = 600,
       width = 7,
       height = 3)



