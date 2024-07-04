# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            abundance.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-04                                       ║
# ║ Last Modified  : 2024-07-04                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(readr)
library(tidyverse)
library(ggplot2)


## Functions

processAncomCsv <- function(file_path) {
  df <- read_csv(file_path) %>% 
    select(-matches("\\(Intercept\\)")) %>%
    rename_with(~ if_else(.x == "id", .x, paste0(.x,
                                                 "_",
                                                 sub("_slice\\.csv$",
                                                     "",
                                                     basename(file_path)))))
  return(df)
}

import_ancombc <- function(qza_path) {
  csv_name <- unzip(qza_path, list = TRUE, exdir = tempdir()) %>%
    filter(str_detect(.data$Name, ".csv")) %>% {.$Name}
  unzip(qza_path, files = csv_name, exdir = tempdir())
  csv_files <- setNames(lapply(file.path(tempdir(), csv_name), processAncomCsv),
                        sub("_slice\\.csv$", "", basename(csv_name)))
  merged_df <- reduce(csv_files, full_join, by = "id")
  return(merged_df)
}


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

location_file_path <- file.path(project_dir,
                               "qiime2/abundance/Location/filtered_ancombc.qza")

location_ancombc <- import_ancombc(location_file_path)

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)

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
