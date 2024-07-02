# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            pcoa.R                                 ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-02                                       ║
# ║ Last Modified  : 2024-07-02                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(tidyverse)
library(qiime2R)

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

jaccard_file_path <- file.path(project_dir,
                             "qiime2/diversity/jaccard_pcoa_results.qza")

jaccard <- read_qza(jaccard_file_path)

metadata <- read.csv(file.path(cluster_path,
                                    "home/salias/projects/sporeflow/metadata.tsv"),
                     sep = "\t")

## Get variance explained by each PCo 

jaccard_pco1 <- round(jaccard[["data"]]$ProportionExplained$PC1 * 100, 2)
jaccard_pco2 <- round(jaccard[["data"]]$ProportionExplained$PC2 * 100, 2)


## Plot PCoA

pdf(file.path(outdir, "pcoa_jaccard.pdf"))

jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Location`, shape=`Fertilization`)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(16,1,2), name="Fertilization") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Location") +
  xlab(paste0("PC1 | ", jaccard_pco1, "% of variance explained")) +
  ylab(paste0("PC2 | ", jaccard_pco2, "% of variance explained"))

dev.off()
