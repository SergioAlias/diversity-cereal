# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   cross_domain_interactions.R                     ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-01-15                                       ║
# ║ Last Modified  : 2025-02-19                                       ║
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
library(pheatmap)
library(gridExtra)

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

doPheatmap <- function(mat, gaps_row, gaps_col) {
  return(pheatmap(mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           color = colorRampPalette(c("orange", "white", "purple"))(100),
           legend = TRUE,
           scale = "none",
           show_rownames = TRUE,
           show_colnames = TRUE,
           angle_col = 45,
           fontsize = 12,
           fontsize_row = 10,
           fontsize_col = 10,
           legend_breaks = c(-0.95, -0.5, 0, 0.5, 0.95),
           legend_labels = c(-1, -0.5, 0, 0.5, 1))
  )
}


## Import QIIME 2 files

set.seed(1234)
n_cores <- 10

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

list2env(lapply(list(fungi_table_org = fungi_table_org, 
                     bacteria_table_org = bacteria_table_org, 
                     fungi_table_con = fungi_table_con, 
                     bacteria_table_con = bacteria_table_con),
                function(mat) mat[, colSums(mat) != 0]),
         envir = .GlobalEnv)

sparxcc_org <- SparXCC_base(fungi_table_org,
                            bacteria_table_org,
                            pseudo_count = 1,
                            cores = n_cores)

sparxcc_con <- SparXCC_base(fungi_table_con,
                            bacteria_table_con,
                            pseudo_count = 1,
                            cores = n_cores)


cor_org <- sparxcc_org[["cor"]]
cor_con <- sparxcc_con[["cor"]]

fungi_to_keep <- c("Aspergillus",
                   "Penicillium",
                   "Fusarium",
                   "Alternaria")

bacteria_to_keep <- c("Bacillus",
                      "Arthrobacter",
                      "Pseudoarthrobacter",
                      "Streptomyces")


taxa_fungi <- parse_taxonomy(read_qza(file.path(cluster_path,
                                       base_path,
                                       fungi_project,
                                       "qiime2/taxonomy/taxonomy.qza"))[["data"]])

taxa_bacteria <- parse_taxonomy(read_qza(file.path(cluster_path,
                                                   base_path,
                                                   bacteria_project,
                                                   "qiime2/taxonomy/taxonomy.qza"))[["data"]])

rownames(cor_con) <- gsub("_", " ",
                          paste0(taxa_fungi$Species[match(rownames(cor_con),
                                                          rownames(taxa_fungi))], 
                                 " (",
                                 substr(rownames(cor_con), 1, 7),
                                 ")"))

colnames(cor_con) <- gsub("_", " ",
                          paste0(taxa_bacteria$Genus[match(colnames(cor_con),
                                                           rownames(taxa_bacteria))],
                                 " - ",
                                 taxa_bacteria$Species[match(colnames(cor_con),
                                                           rownames(taxa_bacteria))], 
                                 " (",
                                 substr(colnames(cor_con), 1, 7),
                                 ")"))

rownames(cor_org) <- gsub("_", " ",
                          paste0(taxa_fungi$Species[match(rownames(cor_org),
                                                          rownames(taxa_fungi))], 
                                 " (",
                                 substr(rownames(cor_org), 1, 7),
                                 ")"))

colnames(cor_org) <- gsub("_", " ",
                          paste0(taxa_bacteria$Genus[match(colnames(cor_org),
                                                           rownames(taxa_bacteria))],
                                 " - ",
                                 taxa_bacteria$Species[match(colnames(cor_org),
                                                             rownames(taxa_bacteria))], 
                                 " (",
                                 substr(colnames(cor_org), 1, 7),
                                 ")"))

cor_con <- cor_con[grepl(paste(fungi_to_keep, collapse = "|"), rownames(cor_con)),
                   grepl(paste(bacteria_to_keep, collapse = "|"), colnames(cor_con))]

cor_org <- cor_org[grepl(paste(fungi_to_keep, collapse = "|"), rownames(cor_org)),
                   grepl(paste(bacteria_to_keep, collapse = "|"), colnames(cor_org))]

cor_con <- cor_con[order(rownames(cor_con)), order(colnames(cor_con))]
cor_org <- cor_org[order(rownames(cor_org)), order(colnames(cor_org))]

bacterialRenamer <- function(x) {
  if (grepl("- NA", x)) {
    return(gsub("- NA", "sp", x))
  } else {
    return(strsplit(x, " - ")[[1]][2])
  }
}

colnames(cor_con) <- sapply(colnames(cor_con), bacterialRenamer)
colnames(cor_org) <- sapply(colnames(cor_org), bacterialRenamer)

con_plot <- doPheatmap(mat = cor_con,
                       gaps_col = c(4, 21),
                       gaps_row = c(4, 10, 19))

org_plot <- doPheatmap(mat = cor_org,
                       gaps_col = c(4, 20),
                       gaps_row = c(5, 10, 17))


