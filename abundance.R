# ╔═══════════════════════════════════════════════════════════════════╗
# ║                            abundance.R                            ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2024-07-04                                       ║
# ║ Last Modified  : 2024-10-29                                       ║
# ║ GitHub Repo    : https://github.com/SergioAlias/diversity-cereal  ║
# ║ Contact        : salias[at]ucm[dot]es                             ║
# ╚═══════════════════════════════════════════════════════════════════╝

## Libraries

library(magrittr, include.only = "%<>%")
library(qiime2R)
library(readr)
library(tidyverse)
library(EnhancedVolcano)
library(patchwork)

## Functions

source("/home/sergio/projects/diversity-cereal/parse_ancombc.R")

volcanoFromAncombc <- function(qza_path,
                               log2fc_col,
                               pval_col,
                               up_color,
                               down_color,
                               up_shape,
                               down_shape,
                               up_legend,
                               down_legend,
                               ...,
                               lab_col = NA,
                               log2fc_cutoff = 2,
                               pval_cutoff = 0.05,
                               ns_color = "black",
                               ns_shape = 4,
                               ns_legend = "NS",
                               taxonomy_df = taxonomy)
{
  ancombc <- import_ancombc(qza_path)
  ancombc %<>% left_join(taxonomy_df)
  
  keyvals_col <- ifelse(
    ancombc[[log2fc_col]] < -log2fc_cutoff & ancombc[[pval_col]] < pval_cutoff, down_color,
    ifelse(ancombc[[log2fc_col]] > log2fc_cutoff & ancombc[[pval_col]] < pval_cutoff, up_color,
           ns_color))
  keyvals_col[is.na(keyvals_col)] <- ns_color
  names(keyvals_col)[keyvals_col == up_color] <- up_legend
  names(keyvals_col)[keyvals_col == ns_color] <- ns_legend
  names(keyvals_col)[keyvals_col == down_color] <- down_legend
  
  keyvals_shape <- keyvals_col
  keyvals_shape[keyvals_shape == up_color] <- up_shape
  keyvals_shape[keyvals_shape == ns_color] <- ns_shape
  keyvals_shape[keyvals_shape == down_color] <- down_shape
  keyvals_shape %<>% as.integer()
  names(keyvals_shape) <- names(keyvals_col)
  
  v_plot <- ancombc %>%
    EnhancedVolcano(lab = {{lab_col}},
                    x = {{log2fc_col}},
                    y = {{pval_col}},
                    pCutoff = pval_cutoff,
                    FCcutoff = log2fc_cutoff,
                    colCustom = keyvals_col,
                    shapeCustom = keyvals_shape,
                    ...) +
    guides(color = guide_legend("Combined Legend",
                                override.aes = list(alpha=1)),
           shape = guide_legend("Combined Legend")) +
    theme_classic() +
    theme(legend.title=element_blank(),
          legend.position="top")
  
  return(v_plot)
}


## Import QIIME 2 files

project_name <- "cereal_16S"
out <- "diversity-cereal-16S"
eco_tag <- "ECO" # ORG (ITS)
con_tag <- "CON" # MFI (ITS)
rot_tag <- "ROT"

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")
project_dir <- file.path(cluster_path,
                         "scratch/salias/projects",
                         project_name)
outdir <- file.path("/home/sergio/scratch",
                    out,
                    "abundance")

treatment_rot_file_path <- file.path(project_dir,
                                     paste0("qiime2/abundance/Fertilization_",
                                            rot_tag,
                                            "/filtered_ancombc.qza"))
treatment_con_file_path <- file.path(project_dir,
                                     paste0("qiime2/abundance/Fertilization_",
                                            con_tag,
                                            "/filtered_ancombc.qza"))

taxonomy_file_path <- file.path(project_dir,
                                "qiime2/taxonomy/taxonomy.qza")

taxonomy <- read_qza(taxonomy_file_path)$data

taxonomy %<>% parse_taxonomy() %>% rownames_to_column("id")

## Colors and shapes

source("/home/sergio/projects/diversity-cereal/colors.R")


## Volcano plots

### ECO vs ROT

pdf(file.path(outdir, "volcano_treatment_eco_vs_rot.pdf"))

v_treatment_eco_vs_rot <- volcanoFromAncombc(qza_path = treatment_rot_file_path,
                             log2fc_col = paste0("Fertilization", eco_tag, "_lfc"),
                             pval_col = paste0("Fertilization", eco_tag, "_q_val"),
                             up_color = treatment_colors[["ECO"]],
                             down_color = treatment_colors[["ROT"]],
                             up_shape = treatment_shapes[["ECO"]],
                             down_shape = treatment_shapes[["ROT"]],
                             up_legend = "DA (ECO)",
                             down_legend = "DA (ROT)",
                             ylab = bquote(~Log[10]~ "Q-value"),
                             title = NULL,
                             subtitle = NULL,
                             caption = NULL,
                             xlim = c(-10, 10),
                             ylim = c(0, 300))

v_treatment_eco_vs_rot

dev.off()

### CON vs ROT

pdf(file.path(outdir, "volcano_treatment_con_vs_rot.pdf"))

v_treatment_con_vs_rot <- volcanoFromAncombc(qza_path = treatment_rot_file_path,
                                      log2fc_col = paste0("Fertilization", con_tag, "_lfc"),
                                      pval_col = paste0("Fertilization", con_tag, "_q_val"),
                                      up_color = treatment_colors[["CON"]],
                                      down_color = treatment_colors[["ROT"]],
                                      up_shape = treatment_shapes[["CON"]],
                                      down_shape = treatment_shapes[["ROT"]],
                                      up_legend = "DA (CON)",
                                      down_legend = "DA (ROT)",
                                      ylab = bquote(~Log[10]~ "Q-value"),
                                      title = NULL,
                                      subtitle = NULL,
                                      caption = NULL,
                                      xlim = c(-10, 10),
                                      ylim = c(0, 300))

v_treatment_con_vs_rot

dev.off()

### ECO vs CON

pdf(file.path(outdir, "volcano_treatment_eco_vs_con.pdf"))

v_treatment_eco_vs_con <- volcanoFromAncombc(qza_path = treatment_con_file_path,
                                             log2fc_col = paste0("Fertilization", eco_tag, "_lfc"),
                                             pval_col = paste0("Fertilization", eco_tag, "_q_val"),
                                             up_color = treatment_colors[["ECO"]],
                                             down_color = treatment_colors[["CON"]],
                                             up_shape = treatment_shapes[["ECO"]],
                                             down_shape = treatment_shapes[["CON"]],
                                             up_legend = "DA (ECO)",
                                             down_legend = "DA (CON)",
                                             ylab = bquote(~Log[10]~ "Q-value"),
                                             title = NULL,
                                             subtitle = NULL,
                                             caption = NULL,
                                             xlim = c(-10, 10),
                                             ylim = c(0, 300))

v_treatment_eco_vs_con

dev.off()

### Grouped plots

pdf(file.path(outdir, "patched_abundance_treatment.pdf"),
    width = 12)

(v_treatment_eco_vs_rot +
    v_treatment_con_vs_rot +
    v_treatment_eco_vs_con &
    theme(plot.tag.position = "topleft")) +
  plot_layout(axis_titles = "collect") + # ,
              # guides = "collect") +
  plot_annotation(tag_levels = 'A')

dev.off()


