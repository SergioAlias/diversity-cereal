# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            abundance.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-04                                       ║
# ║ Last Modified  : 2024-07-08                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(readr)
library(tidyverse)
library(EnhancedVolcano)

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
outdir <- "/home/sergio/scratch/diversity-cereal/abundance"

treatment_file_path <- file.path(project_dir,
                               "qiime2/abundance/Fertilization/level_7_ancombc.qza")

treatment_ancombc <- import_ancombc(treatment_file_path)

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)$data

taxonomy %<>% parse_taxonomy() %>% rownames_to_column("id")

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

treatment_ancombc %<>%
  mutate(Species = sapply(str_split(id, ";"), function(x) tail(x, 1)))
# treatment_ancombc %<>% left_join(taxonomy)

micotoxin_producers <- c("Fusarium",
                         "Aspergillus",
                         "Penicillium",
                         "Alternaria",
                         "Claviceps")


contains_producer <- function(taxa, producers) {
  if (any(str_detect(taxa, producers))) {
    return(taxa %>% str_remove("^[a-z]__") %>% str_replace_all("_", " "))
  } else {
    return(NA)
  }
}

treatment_ancombc %<>%
  mutate(Species = sapply(Species, contains_producer, producers = micotoxin_producers))

## Prepare for volcano

logfc_thr <- 2
qval_thr <- 0.05

filtered_ancombc <- treatment_ancombc[(treatment_ancombc$FertilizationORG_lfc > logfc_thr |
                                         treatment_ancombc$FertilizationORG_lfc < -logfc_thr) &
                                        treatment_ancombc$FertilizationORG_q_val < qval_thr, ]

sorted_desc <- filtered_ancombc[order(-filtered_ancombc$FertilizationORG_lfc), ]
sorted_asc <- filtered_ancombc[order(filtered_ancombc$FertilizationORG_lfc), ]

top_da <- head(sorted_desc, 10)
bottom_da <- head(sorted_asc, 10)

top_features <- c(top_da$Species, bottom_da$Species)


## Volcano plots

pdf(file.path(outdir, "volcano_treatment_eco.pdf"))

v_treatment_eco <- treatment_ancombc %>%
  EnhancedVolcano(lab = treatment_ancombc$Species,
                  x = "FertilizationORG_lfc",
                  y = "FertilizationORG_q_val",
                  selectLab = top_features,
                  pCutoff = 0.05,
                  FCcutoff = 2,
                  title = NULL,
                  subtitle = NULL,
                  caption = NULL,
                  ylab = bquote(~Log[10]~ "Q-value"),
                  legendPosition = "top",
                  drawConnectors = TRUE,
                  widthConnectors = 0.75,
                  min.segment.length = 1)

v_treatment_eco

dev.off()
