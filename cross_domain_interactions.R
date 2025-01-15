# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   cross_domain_interactions.R                     ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-01-15                                       ║
# ║ Last Modified  : 2025-01-15                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
  library(devtools)
  install_github("zdk123/SpiecEasi")
}

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(SpiecEasi)
library(qiime2R)

## Import QIIME 2 files

fungi_project <- "micofood_24"
bacteria_project <- "cereal_16S"
base_path <- "scratch/salias/projects"
table_path <- "qiime2/feature_tables/level_6_table.qza"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")

fungi_table <- file.path(cluster_path,
                         base_path,
                         fungi_project,
                         table_path)

bacteria_table <- file.path(cluster_path,
                           base_path,
                           bacteria_project,
                           table_path)

fungi_table %<>% read_qza() %>%
  {.$data} %>%
  as.data.frame()
