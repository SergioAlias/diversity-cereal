# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   cross_domain_interactions.R                     ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-01-15                                       ║
# ║ Last Modified  : 2025-02-28                                       ║
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
library(patchwork)

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
           angle_col = 90,
           fontsize = 10,
           fontsize_row = 6,
           fontsize_col = 6,
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
                      "Pseudarthrobacter",
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
                       gaps_col = c(10, 62, 64),
                       gaps_row = c(7, 19, 36))

org_plot <- doPheatmap(mat = cor_org,
                       gaps_col = c(9, 65, 68),
                       gaps_row = c(6, 14, 26))

## Permutation analysis

nperm = 10000

permutation_test <- function(con, org, n_perm = nperm) {
  obs_diff <- mean(con) - mean(org)
  all_values <- c(con, org)
  group_labels <- c(rep(1, length(con)), rep(2, length(org)))
  
  perm_diffs <- replicate(n_perm, {
    shuffled_labels <- sample(group_labels)
    mean_1 <- mean(all_values[shuffled_labels == 1])
    mean_2 <- mean(all_values[shuffled_labels == 2])
    mean_1 - mean_2
  })
  
  p_value_two_tailed <- mean(abs(perm_diffs) >= abs(obs_diff))
  p_value_greater <- mean(perm_diffs >= obs_diff)
  p_value_less <- mean(perm_diffs <= obs_diff)
  
  return(list(obs_diff = obs_diff, p_value_two_tailed = p_value_two_tailed,
              p_value_greater = p_value_greater, p_value_less = p_value_less))
}

### Overall Comparison

cor_con_values <- as.vector(cor_con)
cor_org_values <- as.vector(cor_org)

overall_test <- permutation_test(cor_con_values, cor_org_values)

cat("===== Overall Comparison =====\n")
cat("Observed Difference in Means:", overall_test$obs_diff, "\n")
cat("Two-tailed P-Value:", overall_test$p_value_two_tailed, "\n")
cat("P-Value (cor_con > cor_org):", overall_test$p_value_greater, "\n")
cat("P-Value (cor_con < cor_org):", overall_test$p_value_less, "\n\n")
if (overall_test$p_value_greater < 0.05) {
  cat("Conclusion: cor_con correlations are significantly larger than cor_org.\n")
} else if (overall_test$p_value_less < 0.05) {
  cat("Conclusion: cor_con correlations are significantly smaller than cor_org.\n")
} else {
  cat("Conclusion: No significant difference in direction.\n")
}

overall_df_pair <- data.frame(
  Correlation = c(cor_con_values, cor_org_values),
  Condition = rep(c("CON", "ORG"), times = c(length(cor_con_values), length(cor_org_values)))
)

overall_p <- ggplot(overall_df_pair, aes(x = Correlation, fill = Condition, color = Condition)) +
  geom_density(alpha = 0.6, size = 1) +
  theme_classic(base_size = 14) +
  labs(x = "Correlation",
       y = "Density") +
  scale_fill_manual(values = c("#C4A484", "lightgreen")) +
  scale_color_manual(values = c("brown", "darkgreen")) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.size = unit(1.25, "lines"),
    axis.text = element_text(size = 12),
  )

### Per-Fungi-Bacteria Pair Comparison

results <- data.frame(Fungus = character(), Bacteria = character(),
                      Obs_Diff = numeric(), P_TwoTailed = numeric(),
                      P_Greater = numeric(), P_Less = numeric(), stringsAsFactors = FALSE)

plot_list <- list()

for (fungus in fungi_to_keep) {
  for (bacterium in bacteria_to_keep) {
    cor_con_values <- as.vector(cor_con[grepl(fungus, rownames(cor_con)), grepl(bacterium, colnames(cor_con))])
    cor_org_values <- as.vector(cor_org[grepl(fungus, rownames(cor_org)), grepl(bacterium, colnames(cor_org))])
    test_result <- permutation_test(cor_con_values, cor_org_values)
    results <- rbind(results, data.frame(Fungus = fungus, Bacteria = bacterium,
                                         Obs_Diff = test_result$obs_diff,
                                         P_TwoTailed = test_result$p_value_two_tailed,
                                         P_Greater = test_result$p_value_greater,
                                         P_Less = test_result$p_value_less))
    df_pair <- data.frame(
      Correlation = c(cor_con_values, cor_org_values),
      Condition = rep(c("CON", "ORG"), times = c(length(cor_con_values), length(cor_org_values)))
    )
    
    x_axis_label <- NULL
    y_axis_label <- NULL
    
    if (length(plot_list) %in% c(0, 4, 8, 12)) {
      y_axis_label <- fungus
    }
    
    if (length(plot_list) %in% c(12, 13, 14, 15)) {
      x_axis_label <- bacterium
    }
    
    p <- ggplot(df_pair, aes(x = Correlation, fill = Condition, color = Condition)) +
      geom_density(alpha = 0.6, size = 1) +
      theme_classic(base_size = 14) +
      labs(x = x_axis_label,
           y = y_axis_label) +
      scale_fill_manual(values = c("#C4A484", "lightgreen")) +
      scale_color_manual(values = c("brown", "darkgreen")) +
      theme(
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(1.25, "lines"),
        axis.title = element_text(face = "italic"),
        axis.text = element_text(size = 12),
      )
    
    
    plot_list[[paste(fungus, bacterium, sep = "_")]] <- p
  }
}

cat("\n===== Per Fungi-Bacteria Pair Results =====\n")
print(results)


wp <- guide_area() / 
    wrap_plots(plot_list) / 
    grid::textGrob("Correlation", gp = grid::gpar(fontsize = 15)) + 
  plot_layout(guides = "collect",
              heights = c(1, 30, 1))

wp <- wrap_elements(grid::textGrob("Density", rot = 90, gp = grid::gpar(fontsize = 15))) + wp +
  plot_layout(widths = c(1, 30))


### Para enseñar ###

con_plot

dev.off()

org_plot

dev.off()

overall_p

View(overall_test)

View(results)

wp

