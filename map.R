# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            map.R                                  ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-08                                       ║
# ║ Last Modified  : 2024-07-08                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(sf)
library(mapSpain)
library(magrittr)
library(tidyverse)
library(ggrepel)


## Load map resources

esp <- esp_get_prov(moveCAN = FALSE)
can_box <- esp_get_can_box()

locations <- data.frame(
  lng = c(-4.21, -5.38, -4.70),
  lat = c(37.15, 41.25, 41.7),
  ID = c("Riofrío", "Fuentelapeña", "Zamadueñas"))

locations %<>% st_as_sf(coords = c("lng", "lat"), remove = FALSE, 
                      crs = 4326, agr = "constant")


## Plot map

ggplot(esp) +
  geom_sf(fill = "white") +
  geom_sf(data = locations) +
  geom_label_repel(data = locations, aes(x = lng, y = lat, label = ID), 
                  fontface = "bold", nudge_x = c(1, -1.5, 2),
                  nudge_y = c(0.25, -0.25, 0.5)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(-11, 6), ylim = c(35.5, 44.5), expand = FALSE)

