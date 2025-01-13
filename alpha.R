# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            alpha.R                                ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-02                                       ║
# ║ Last Modified  : 2025-01-13                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(ggpubr)
library(patchwork)


## Import QIIME 2 files

project_name <- "micofood_24" # micofood_24 or cereal_16S
amplicon <- "ITS" # ITS or 16S
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
                    "alpha")

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

shannon_file_path <- file.path(project_dir,
                               "qiime2/diversity/shannon_vector.qza")
simpson_file_path <- file.path(project_dir,
                               "qiime2/diversity/simpson_vector.qza")
chao1_file_path <- file.path(project_dir,
                               "qiime2/diversity/chao1_vector.qza")

shannon <- read_qza(shannon_file_path)
shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(shannon)

simpson <- read_qza(simpson_file_path)
simpson <- simpson$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(simpson)

chao1 <- read_qza(chao1_file_path)
chao1 <- chao1$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(chao1)

## Colors and shapes

source("/home/sergio/projects/diversity-cereal/colors.R")


## Comparisons

comparisons_treatment <- combn(unique(metadata$Treatment), 2, simplify = FALSE)
comparisons_location <- combn(unique(metadata$Location), 2, simplify = FALSE)


## Alpha boxplots

if (amplicon == "ITS"){
  shannon_pos_stat <- 7.8
  shannon_pos_stat_paired <- c(7.3, 7.55, 7.4)
  simpson_pos_stat <- 0.993
  simpson_pos_stat_paired <- c(0.982, 0.988, 0.9845)
  chao1_pos_stat <- 820
  chao1_pos_stat_paired <- c(730, 775, 750)
}
if (amplicon == "16S"){
  shannon_pos_stat <- 11.2
  shannon_pos_stat_paired <- c(10.55, 10.7, 10.625)
  simpson_pos_stat <- 0.99885
  simpson_pos_stat_paired <- c(0.99852, 0.99866, 0.99859)
  chao1_pos_stat <- 5500
  chao1_pos_stat_paired <- c(3100, 3500, 3300)
}

### Shannon

pdf(file.path(outdir, "shannon_treatment.pdf"))

shannon_t <- metadata %>%
  ggboxplot("Treatment", "shannon_entropy",
            color = "Treatment",
            fill = "Treatment",
            alpha = 0.1,
            palette = treatment_colors,
            add = "jitter",
            shape = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  ylab("Shannon") +
  stat_compare_means(label.y = shannon_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_treatment,
                     label.y = shannon_pos_stat_paired)

shannon_t

dev.off()


pdf(file.path(outdir, "shannon_location.pdf"))

shannon_l <- metadata %>%
  ggboxplot("Location", "shannon_entropy",
            color = "Location",
            fill = "Location",
            alpha = 0.1,
            palette = location_colors,
            add = "jitter",
            shape = "Location") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  ylab("Shannon") +
  stat_compare_means(label.y = shannon_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_location,
                     label.y = shannon_pos_stat_paired)

shannon_l

dev.off()

### Simpson

pdf(file.path(outdir, "simpson_treatment.pdf"))

simpson_t <- metadata %>%
  ggboxplot("Treatment", "simpson",
            color = "Treatment",
            fill = "Treatment",
            alpha = 0.1,
            palette = treatment_colors,
            add = "jitter",
            shape = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  ylab("Inverse Simpson") +
  stat_compare_means(label.y = simpson_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_treatment,
                     label.y = simpson_pos_stat_paired)

simpson_t

dev.off()


pdf(file.path(outdir, "simpson_location.pdf"))

simpson_l <- metadata %>%
  ggboxplot("Location", "simpson",
            color = "Location",
            fill = "Location",
            alpha = 0.1,
            palette = location_colors,
            add = "jitter",
            shape = "Location") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  ylab("Inverse Simpson") +
  stat_compare_means(label.y = simpson_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_location,
                     label.y = simpson_pos_stat_paired)

simpson_l

dev.off()

### Chao1

pdf(file.path(outdir, "chao1_treatment.pdf"))

chao1_t <- metadata %>%
  ggboxplot("Treatment", "chao1",
            color = "Treatment",
            fill = "Treatment",
            alpha = 0.1,
            palette = treatment_colors,
            add = "jitter",
            shape = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  ylab("Chao1") +
  stat_compare_means(label.y = chao1_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_treatment,
                     label.y = chao1_pos_stat_paired)

chao1_t

dev.off()


pdf(file.path(outdir, "chao1_location.pdf"))

chao1_l <- metadata %>%
  ggboxplot("Location", "chao1",
            color = "Location",
            fill = "Location",
            alpha = 0.1,
            palette = location_colors,
            add = "jitter",
            shape = "Location") +
  scale_shape_manual(values = location_shapes, name = "Location") +
  ylab("Chao1") +
  stat_compare_means(label.y = chao1_pos_stat) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_location,
                     label.y = chao1_pos_stat_paired)

chao1_l

dev.off()

### Grouped plots

pdf(file.path(outdir, "patched_treatment.pdf"),
    width = 8.7,
    height = 4.5)

(chao1_t + theme(legend.position="none") +
 shannon_t + theme(legend.position="none") +
 simpson_t + theme(legend.position="none") &
    theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = 'A')

dev.off()


pdf(file.path(outdir, "patched_location.pdf"),
    width = 8.7,
    height = 4.5)

(chao1_l + theme(legend.position="none") +
    shannon_l + theme(legend.position="none") +
    simpson_l + theme(legend.position="none") &
    theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = 'A')

dev.off()

