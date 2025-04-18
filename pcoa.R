# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            pcoa.R                                 ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-02                                       ║
# ║ Last Modified  : 2025-04-03                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(rlang)
library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(patchwork)


## Import QIIME 2 files

project_name <- "micofood_24" # micofood_24 or cereal_16S
local_metadata <- "diversity-cereal" # diversity-cereal or diversity-cereal-16S
out <- "paper_ready_diversity"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir <- file.path(cluster_path,
                         "scratch/salias/projects",
                         project_name)
outdir <- file.path("/home/sergio/scratch",
                 out,
                 "beta")

jaccard_file_path <- file.path(project_dir,
                             "qiime2/diversity/jaccard_pcoa_results.qza")
bray_curtis_file_path <- file.path(project_dir,
                                   "qiime2/diversity/bray_curtis_pcoa_results.qza")
aitchison_file_path <- file.path(project_dir,
                                 "qiime2/diversity/aitchison_pcoa_results.qza")

jaccard <- read_qza(jaccard_file_path)
bray_curtis <- read_qza(bray_curtis_file_path)
aitchison <- read_qza(aitchison_file_path)

metadata <- read.csv(file.path("/home/sergio/scratch",
                               local_metadata,
                               "metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

metadata %<>%
  mutate(Treatment = case_when(
    Fertilization == "MFI" ~ "CON",
    Fertilization == "ORG" ~ "ORG",
    Fertilization == "ROT" ~ "PLO",
    Fertilization == "CON" ~ "CON",
    Fertilization == "ECO" ~ "ORG"
  ))


## Get variance explained by each PCo 

jaccard_pco1 <- round(jaccard[["data"]]$ProportionExplained$PC1 * 100, 2)
jaccard_pco2 <- round(jaccard[["data"]]$ProportionExplained$PC2 * 100, 2)

bray_curtis_pco1 <- round(bray_curtis[["data"]]$ProportionExplained$PC1 * 100, 2)
bray_curtis_pco2 <- round(bray_curtis[["data"]]$ProportionExplained$PC2 * 100, 2)

aitchison_pco1 <- round(aitchison[["data"]]$ProportionExplained$PC1 * 100, 2)
aitchison_pco2 <- round(aitchison[["data"]]$ProportionExplained$PC2 * 100, 2)


## Colors and shapes

source("/home/sergio/projects/diversity-cereal/colors.R")


## PCoA plots

### Jaccard

pdf(file.path(outdir, "pcoa_jaccard_location.pdf"),
    height = 5,
    width = 7)

p_j1 <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Location`, shape = `Location`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = location_colors, name = "Location") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  xlab(paste0("PCo-1 | ", jaccard_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", jaccard_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Location`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_j1

dev.off()


pdf(file.path(outdir, "pcoa_jaccard_treatment.pdf"),
    height = 5,
    width = 7)

p_j2 <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Treatment`, shape = `Treatment`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  xlab(paste0("PCo-1 | ", jaccard_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", jaccard_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Treatment`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_j2

dev.off()

pdf(file.path(outdir, "pcoa_jaccard_movement.pdf"),
    height = 5,
    width = 7)

p_j3 <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Movement`, shape = `Movement`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = movement_colors, name = "Movement") +
  scale_shape_manual(values = movement_shapes, name = "Movement") +
  xlab(paste0("PCo-1 | ", jaccard_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", jaccard_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Movement`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_j3

dev.off()

### Bray-Curtis

pdf(file.path(outdir, "pcoa_bray_curtis_location.pdf"),
    height = 5,
    width = 7)

p_b1 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Location`, shape = `Location`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = location_colors, name = "Location") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  xlab(paste0("PCo-1 | ", bray_curtis_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", bray_curtis_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Location`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_b1

dev.off()


pdf(file.path(outdir, "pcoa_bray_curtis_treatment.pdf"),
    height = 5,
    width = 7)

p_b2 <- bray_curtis$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Treatment`, shape = `Treatment`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  xlab(paste0("PCo-1 | ", bray_curtis_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", bray_curtis_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Treatment`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_b2

dev.off()


### Aitchison

pdf(file.path(outdir, "pcoa_aitchison_location.pdf"),
    height = 5,
    width = 7)

p_a1 <- aitchison$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Location`, shape = `Location`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = location_colors, name = "Location") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  xlab(paste0("PCo-1 | ", aitchison_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", aitchison_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Location`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_a1

dev.off()


pdf(file.path(outdir, "pcoa_aitchison_treatment.pdf"),
    height = 5,
    width = 7)

p_a2 <- aitchison$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = PC1, y = PC2, color = `Treatment`, shape = `Treatment`)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  xlab(paste0("PCo-1 | ", aitchison_pco1, "% of variance explained")) +
  ylab(paste0("PCo-2 | ", aitchison_pco2, "% of variance explained")) +
  stat_ellipse(aes(group = `Treatment`)) +
  guides(color = guide_legend(override.aes = list(linetype = 0, alpha=1)))

p_a2

dev.off()


### Grouped plots

pdf(file.path(outdir, "patched_pcoa_jaccard.pdf"),
    width = 8,
    height = 3.5)

(p_j2 + p_j1 &
    theme(plot.tag.position = "topright")) +
  plot_layout(axis_titles = "collect",
              guides = "collect") # +
  # plot_annotation(tag_levels = 'A')

dev.off()

pdf(file.path(outdir, "patched_pcoa_bray_curtis.pdf"),
    width = 8,
    height = 3.5)

(p_b2 + p_b1 &
  theme(plot.tag.position = "topright")) +
  plot_layout(axis_titles = "collect",
              guides = "collect") # +
  # plot_annotation(tag_levels = 'A')
  
dev.off()

pdf(file.path(outdir, "patched_pcoa_aitchison_vertical.pdf"),
    width = 5,
    height = 8.5)

(p_a2 / p_a1 &
    theme(plot.tag.position = "topright")) +
  plot_layout(axis_titles = "collect") # ,
              # guides = "collect")  # +
  # plot_annotation(tag_levels = 'A')

dev.off()


### Supplementary figures plots

j_both <- (p_j2 + p_j1) + 
  plot_layout(axis_titles = "collect")

b_both <- (p_b2 + p_b1) +
  plot_layout(axis_titles = "collect")


pdf(file.path(outdir, "patched_supp_fungi.pdf"),
    width = 8,
    height = 7)

(j_both / b_both & theme(plot.tag.position = "topleft")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "", "B", "")))

dev.off()

# PNG equivalent

png(file.path(outdir, "patched_supp_fungi.png"),
    width = 8,
    height = 7,
    units = "in",
    res = 300)

(j_both / b_both & theme(plot.tag.position = "topleft")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "", "B", "")))

dev.off()
