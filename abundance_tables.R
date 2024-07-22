# ╔═══════════════════════════════════════════════════════════════════╗
# ║                       abundance_tables.R                          ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-19                                       ║
# ║ Last Modified  : 2024-07-19                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(tidyverse)
library(openxlsx)

## Functions

source("/home/sergio/projects/diversity-cereal/parse_ancombc.R")

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
outdir <- "/home/sergio/scratch/diversity-cereal/abundance"

treatment_rot_file_path <- file.path(project_dir,
                                     "qiime2/abundance/Fertilization_ROT/filtered_ancombc.qza")
treatment_con_file_path <- file.path(project_dir,
                                     "qiime2/abundance/Fertilization_MFI/filtered_ancombc.qza")

rot <- import_ancombc(treatment_rot_file_path)
con <- import_ancombc(treatment_con_file_path)

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)$data

taxonomy_parsed <- taxonomy %>% parse_taxonomy() %>% rownames_to_column("id")

colnames(taxonomy)[1] <- "id"

taxonomy <- taxonomy[, c("id", "Confidence")]

taxonomy_parsed %<>% left_join(taxonomy)
rot %<>% left_join(taxonomy_parsed)
con %<>% left_join(taxonomy_parsed)

## Table filtering

qval_thr <- 0.05
lfc_thr <- 2

ECOvsROT <- rot[,!grepl("MFI", colnames(rot))]
CONvsROT <- rot[,!grepl("ORG", colnames(rot))]
ECOvsCON <- con[,!grepl("ROT", colnames(con))]


ECOvsROT_up <- ECOvsROT %>%
  rename_with(~str_remove(.x, "FertilizationORG_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

ECOvsROT_down <- ECOvsROT %>%
  rename_with(~str_remove(.x, "FertilizationORG_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)


CONvsROT_up <- CONvsROT %>%
  rename_with(~str_remove(.x, "FertilizationMFI_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

CONvsROT_down <- CONvsROT %>%
  rename_with(~str_remove(.x, "FertilizationMFI_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)


ECOvsCON_up <- ECOvsCON %>%
  rename_with(~str_remove(.x, "FertilizationORG_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc >= lfc_thr) %>%
  arrange(desc(lfc))

ECOvsCON_down <- ECOvsCON %>%
  rename_with(~str_remove(.x, "FertilizationORG_")) %>%
  filter(q_val <= qval_thr) %>%
  filter(lfc <= -lfc_thr) %>%
  arrange(lfc)


write.xlsx(list("ECOvsROT_up" = ECOvsROT_up,
                "ECOvsROT_down" = ECOvsROT_down,
                "CONvsROT_up" = CONvsROT_up,
                "CONvsROT_down" = CONvsROT_down,
                "ECOvsCON_up" = ECOvsCON_up,
                "ECOvsCON_down" = ECOvsCON_down),
           file = file.path(outdir, "abundance.xlsx"))
