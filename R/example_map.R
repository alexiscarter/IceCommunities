#
##
### Example for mapping IceCommunities sites
## Alexis Carteron
# 25/05/2021

##
## Load packages packages and data ####
##
library(ggplot2)
library(ggedit)
library(ggrepel)
library(ggspatial)
library(tidyverse)
library(dplyr)
library(rnaturalearth)

## Sites
sites <- read.csv(file = "data/IC_sampled_points_may2021.csv") # Matrix from Alessia 17/05/2021

##
## Selection and manipulation ####
##
str(sites)
sites$dating <- as.numeric(sites$dating)
sites$spot <- as.factor(sites$spot)
sites$date <- as.Date(sites$date, "%d/%m/%Y") # Entered as d/m/y

## Remove columns with technical information, keep 'notes'
sites = select(sites, -time, -logger_n:-X)

## Select all sites to map
sites_map <- sites %>%
  dplyr::filter(!is.na(lon)) %>% # Only sites that have available coordinates
  group_by(glacier) %>% 
  summarise_all(mean)

## Select only sites in the Alps for further mapping
sites_alps_map <- sites %>%
  select(glacier, lon, lat) %>%
  dplyr::filter(lon > 0 & lon < 20 & lat > 44 & lat < 50) %>% 
  dplyr::filter(!is.na(lon)) %>% 
  group_by(glacier) %>% 
  summarise_all(mean)

##
## Make maps ####
##

## Select world frontiers
world <- ne_countries(scale = "medium", returnclass = "sf")

sites_map_simple <- ggplot(data = world) +
  geom_sf(color = "darkgrey", fill = "white") +
  geom_point(data = sites_map, mapping = aes(x = lon, y = lat), shape = 21, size = 2, alpha = .7, color = "black", fill = "red") +
  #geom_label_repel(data = glaciers_map, mapping = aes(x = lon, y = lat, label = glacier), seed = 1, force_pull = 0, max.overlaps = Inf) +
  coord_sf(expand = FALSE) + # Zoom in or out
  labs(x = 'Longitude (°)', y = 'Latitude (°)') + # Label x and y axis
  theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.background = element_rect(fill = 'aliceblue'), text = element_text(family = 'Times', size = 14))
sites_map_simple
#ggsave(file = "figures/sites_map_simple.png", width = 10, height = 8)

## Alps
alps_simple <- ggplot(data = world) +
  geom_sf(color = "darkgrey", fill = "white") +
  geom_point(data = sites_alps_map, mapping = aes(x = lon, y = lat), shape = 21, size = 2, alpha = .7, color = "black", fill = "red") +
  coord_sf(xlim = c(5, 15), ylim = c(44, 48), expand = FALSE, clip = "on") +
  annotate("text", x = 5.8, y = 47.7, label = "Alps") +
  theme(axis.title=element_blank(),axis.ticks=element_blank(),  axis.text=element_blank(), panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.background = element_rect(fill = 'aliceblue'), text = element_text(family = 'Times', size = 14))
alps_simple

## World + Alps inlet
sites_map_simple +
  annotation_custom(ggplotGrob(alps_simple), xmin = 100, xmax = 180, 
                    ymin = 30, ymax = 90)
#ggsave(file = "figures/sites_map_world_alps.png", width = 10, height = 8)
