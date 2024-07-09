# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            map.R                                  ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-08                                       ║
# ║ Last Modified  : 2024-07-09                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(sf)
library(mapSpain)
library(magrittr)
library(tidyverse)
library(ggrepel)

## Setup

outdir <- "/home/sergio/scratch/diversity-cereal/maps"
source("/home/sergio/projects/diversity-cereal/colors.R")

## Load map resources

esp <- esp_get_prov(moveCAN = FALSE)
esp %<>%
  mutate(n_samples = case_when(
    iso2.prov.code == "ES-GR" ~ "1",
    iso2.prov.code == "ES-VA" ~ "2",
    iso2.prov.code == "ES-ZA" ~ "3",
    .default = "0"
  ))

locations <- data.frame(
  lng = c(-4.21, -5.38, -4.70),
  lat = c(37.15, 41.25, 41.7),
  ID = c("Riofrío (Granada)\n· S14 (ROT)",
         "Fuentelapeña (Zamora)\n· S15 (ECO)\n· S16 (ECO)\n· S17 (CON)",
         "Zamadueñas (Valladolid)\n· S18 (CON)\n· S19 (CON)"))

locations %<>% st_as_sf(coords = c("lng", "lat"), remove = FALSE, 
                      crs = 4326, agr = "constant")


## Plot map

pdf(file.path(outdir, "sample_map.pdf"),
    width = 9)

ggplot(esp) +
  geom_sf(aes(fill = n_samples)) +
  scale_fill_discrete(type = sample_map_colors,
                      guide = "none") +
  geom_sf(data = locations) +
  geom_label_repel(data = locations,
                   aes(x = lng, y = lat,
                       label = ID),
                   parse = FALSE,
                   hjust = 0,
                   nudge_x = c(-2.5, 2, 3),
                   nudge_y = c(1.25, -1, 1)) +
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.2),
        panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(-11, 6), ylim = c(35.5, 44.5), expand = TRUE)

dev.off()
