# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            alpha.R                                ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-02                                       ║
# ║ Last Modified  : 2024-07-05                                       ║
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

project_name <- "micofood_24"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir <- file.path(cluster_path,
                         "scratch/salias/projects",
                         project_name)
outdir <- "/home/sergio/scratch/diversity-cereal/alpha"

metadata <- read.csv(file.path(cluster_path,
                               "home/salias/projects/sporeflow/metadata.tsv"),
                     sep = "\t")

colnames(metadata)[1] <- "SampleID"

metadata %<>%
  mutate(Treatment = case_when(
    Fertilization == "MFI" ~ "CON",
    Fertilization == "ORG" ~ "ECO",
    Fertilization == "ROT" ~ "ROT"
  ))

shannon_file_path <- file.path(project_dir,
                               "qiime2/diversity/shannon_vector.qza")
simpson_file_path <- file.path(project_dir,
                               "qiime2/diversity/simpson_vector.qza")

shannon <- read_qza(shannon_file_path)
shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(shannon)

simpson <- read_qza(simpson_file_path)
simpson <- simpson$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(simpson)

## Colors and shapes

source("/home/sergio/projects/diversity-cereal/colors.R")


## Kruskal-Wallis

kw_shannon <- kruskal.test(shannon_entropy ~ Location, data = metadata)
pw_shannon <- pairwise.wilcox.test(x = metadata$shannon_entropy,
                                   g = metadata$Location,
                                   p.adjust.method = "holm")

## Alpha boxplot

### Shannon

comparisons_treatment <- combn(unique(metadata$Treatment), 2, simplify = FALSE)

pdf(file.path(outdir, "shannon_treatment.pdf"))

shannon_t <- metadata %>%
  ggboxplot("Treatment", "shannon_entropy",
            color = "Treatment", palette = treatment_colors,
          add = "jitter", shape = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  ylab("Shannon") +
  stat_compare_means(label.y = 7.7) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_treatment,
                     label.y = c(7.3, 7.5, 7.4))

shannon_t


### Simpson

dev.off()

pdf(file.path(outdir, "simpson_treatment.pdf"))

simpson_t <- metadata %>%
  ggboxplot("Treatment", "simpson",
            color = "Treatment", palette = treatment_colors,
            add = "jitter", shape = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  ylab("Simpson") +
  stat_compare_means(label.y = 0.992) +
  stat_compare_means(aes(label = after_stat(paste0('p = ', p.format, '\n', p.signif))),
                     comparisons = comparisons_treatment,
                     label.y = c(0.982, 0.987, 0.9845))

simpson_t

dev.off()


### Grouped plots

pdf(file.path(outdir, "patched_treatment.pdf"),
    width = 12)

(shannon_t + theme(legend.position="none") + simpson_t + theme(legend.position="none") &
    theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") +
  plot_annotation(tag_levels = 'A')

dev.off()

