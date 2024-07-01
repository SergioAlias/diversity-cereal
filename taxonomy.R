# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            taxonomy.R                             ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-01                                       ║
# ║ Last Modified  : 2024-07-01                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias@ucm.es                                    ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr)
library(file2meco)
library(microeco)
library(ggplot2)
library(ggnested)
library(ggh4x)


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
 
 
dada2_file_path <- file.path(project_dir,
                             "qiime2/feature_tables/filtered_table.qza")
metadata_file_path <- file.path(cluster_path,
                                "home/salias/projects/sporeflow/metadata.tsv")
taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

meco <- qiime2meco(dada2_file_path,
                   sample_table = metadata_file_path,
                   taxonomy_table = taxonomy_file_path)


## Relabel unclassified

# meco$tax_table$Phylum[grepl("__$", meco$tax_table$Phylum)] %<>% paste0(., "Unclassified")
# meco$tax_table$Family[grepl("__$", meco$tax_table$Family)] %<>% paste0(., "Unclassified")
# meco$tax_table$Genus[grepl("__$", meco$tax_table$Genus)] %<>% paste0(., "Unclassified")

## Create trans_abund object

t_family <- trans_abund$new(dataset = meco,
                            taxrank = "Family",
                            ntaxa = 20)

t_stacked_phylum <- trans_abund$new(dataset = meco,
                                    taxrank = "Class",
                                    ntaxa = 20,
                                    delete_taxonomy_prefix = TRUE,
                                    high_level = "Phylum",
                                    prefix = "c__")

## Nested barplot

pdf(file.path(outdir, "barplot_class.pdf"))

t_stacked_phylum$plot_bar(ggnested = TRUE,
                          high_level_add_other = TRUE,
                          xtext_keep = FALSE,
                          # xtext_angle = 90,
                          # xtext_size = 6,
                          facet = c("Type", "Sampling")) + 
  theme(ggh4x.facet.nestline = element_line(colour = "grey95"))

dev.off()

## Boxplot

pdf(file.path(outdir, "boxplot_family.pdf"))

t_family$plot_box(group = "Type",
                  xtext_size = 10,
                  xtext_angle = 30,
                  ytitle_size = 15) + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white", size = 0.5),
    legend.position = c(0.9, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm")
    )

dev.off()

