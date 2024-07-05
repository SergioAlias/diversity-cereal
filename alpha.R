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

library(magrittr)
library(tidyverse)
library(qiime2R)
library(ggpubr)


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
outdir <- "/home/sergio/scratch/diversity-cereal"

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

shannon <- read_qza(shannon_file_path)
shannon <- shannon$data %>% rownames_to_column("SampleID")
metadata %<>% left_join(shannon)


## Colors and shapes

source("/home/sergio/projects/diversity-cereal/colors.R")


## Kruskal-Wallis

kw_shannon <- kruskal.test(shannon_entropy ~ Location, data = metadata)
pw_shannon <- pairwise.wilcox.test(x = metadata$shannon_entropy,
                                   g = metadata$Location,
                                   p.adjust.method = "holm")

## Alpha boxplot

my_comparisons <- combn(unique(metadata$Treatment), 2, simplify = FALSE)

pdf(file.path(outdir, "shannon.pdf"))

metadata %>%
  ggboxplot("Treatment", "shannon_entropy",
            color = "Treatment", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "jitter", shape = "Treatment") +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  ylab("Shannon") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(7.3, 7.5, 7.4),
                     label = "p.signif") +
  stat_compare_means(label.y = 7.7)

dev.off()

