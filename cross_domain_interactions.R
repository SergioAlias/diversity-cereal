# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   cross_domain_interactions.R                     ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-01-15                                       ║
# ║ Last Modified  : 2025-01-28                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
  }

if (!requireNamespace("CompoCor", quietly = TRUE)) {
  library(devtools)
  install_github("IbTJensen/CompoCor")
}

if (!requireNamespace("qiime2R", quietly = TRUE)) {
  library(devtools)
  install_github("jbisanz/qiime2R")
}

library(magrittr, include.only = "%<>%")
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(CompoCor)

## Functions

fromQ2toPhyloTable <- function(project) {
  table <- file.path(cluster_path,
                     base_path,
                     project,
                     table_path)
  table %<>% read_qza() %>%
    {.$data} %>%
    as.matrix()
  colnames(table) <- sapply(strsplit(colnames(table), "_"), function(x) {
    paste(head(x, -1), collapse = "_")
  })
  table %<>% otu_table(taxa_are_rows = TRUE)
  return(table)
}

## Import QIIME 2 files

set.seed(1234)

fungi_project <- "micofood_24"
bacteria_project <- "cereal_16S"
base_path <- "scratch/salias/projects"
table_path <- "qiime2/feature_tables/filtered_table.qza"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")

metadata <- sample_data(read.table(file.path(cluster_path,
                                   "home/salias/projects/cross_domain/metadata.tsv"),
                                   row.names = 1,
                                   sep = "\t",
                                   quote = "",
                                   header = TRUE)
                        )

fungi_table <- fromQ2toPhyloTable(fungi_project)
sample_names(fungi_table) <- sample_names(metadata)

bacteria_table <- fromQ2toPhyloTable(bacteria_project)
sample_names(bacteria_table) <- sample_names(metadata)

fungi_table %<>% as.data.frame() %>% t()
bacteria_table %<>% as.data.frame() %>% t()


fungi_table_org <- fungi_table[grepl("_ORG", rownames(fungi_table)), ]
bacteria_table_org <- bacteria_table[grepl("_ORG", rownames(bacteria_table)), ]

fungi_table_con <- fungi_table[grepl("_CON", rownames(fungi_table)), ]
bacteria_table_con <- bacteria_table[grepl("_CON", rownames(bacteria_table)), ]

sparxcc_org <- SparXCC(fungi_table_org,
                       bacteria_table_org)
