# ╔═══════════════════════════════════════════════════════════════════╗
# ║                   cross_domain_interactions.R                     ║
# ╠═══════════════════════════════════════════════════════════════════╣
# ║ Project        : diversity-cereal                                 ║
# ║ Author         : Sergio Alías-Segura                              ║
# ║ Created        : 2025-01-15                                       ║
# ║ Last Modified  : 2025-01-16                                       ║
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
library(phyloseq)

## Functions

fromQ2toPhyloTable <- function(project, to_keep) {
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
  table <- table[, grep(paste(to_keep, collapse = "|"), colnames(table))] # SUBSET
  table <- table[1:20,] # PARA HACER PRUEBAS
  table %<>% otu_table(taxa_are_rows = TRUE)
  return(table)
}

fromQ2toPhyloTaxa <- function(project, to_keep) {
  taxa <- file.path(cluster_path,
                    base_path,
                    project,
                    taxonomy_path)
  taxa %<>% read_qza() %>%
    {.$data} %>%
    parse_taxonomy() %>%
    as.matrix() %>%
    tax_table()
  return(taxa)
}
  
## Import QIIME 2 files

set.seed(1234)

fungi_project <- "micofood_24"
bacteria_project <- "cereal_16S"
base_path <- "scratch/salias/projects"
table_path <- "qiime2/feature_tables/filtered_table.qza"
taxonomy_path <- "qiime2/taxonomy/taxonomy.qza"

org_samples <- c("15", "16")
con_samples <- c("17", "18", "19")

readRenviron("/home/sergio/Renvs/.RenvBrigit")
brigit_IP <- Sys.getenv("IP_ADDRESS")
cluster_path <- paste0("/run/user/1001/gvfs/sftp:host=",
                       brigit_IP,
                       ",user=salias/mnt/lustre")

fungi_table_org <- fromQ2toPhyloTable(fungi_project, org_samples)
fungi_table_con <- fromQ2toPhyloTable(fungi_project, con_samples)
fungi_taxa <- fromQ2toPhyloTaxa(fungi_project)
fungi_taxa <- fungi_taxa[taxa_names(fungi_table_org),] # PARA HACER PRUEBAS


bacteria_table_org <- fromQ2toPhyloTable(bacteria_project, org_samples)
bacteria_table_con <- fromQ2toPhyloTable(bacteria_project, con_samples)
bacteria_taxa <- fromQ2toPhyloTaxa(bacteria_project)
bacteria_taxa <- bacteria_taxa[taxa_names(bacteria_table_org),] # PARA HACER PRUEBAS

phylo_fungi_org <- phyloseq(fungi_table_org, fungi_taxa)
phylo_fungi_con <- phyloseq(fungi_table_con, fungi_taxa)

phylo_bacteria_org <- phyloseq(bacteria_table_org, bacteria_taxa)
phylo_bacteria_con <- phyloseq(bacteria_table_con, bacteria_taxa)


se_org <- spiec.easi(list(phylo_fungi_org, phylo_bacteria_org),
                      method="mb",
                      nlambda=40,
                      lambda.min.ratio=1e-2,
                      pulsar.params = list(thresh = 0.05))

se_con <- spiec.easi(list(phylo_fungi_con, phylo_bacteria_con),
                     method="mb",
                     nlambda=40,
                     lambda.min.ratio=1e-2,
                     pulsar.params = list(thresh = 0.05))


dtype_org <- c(rep(1,ntaxa(phylo_fungi_org)), rep(2,ntaxa(phylo_bacteria_org)))
plot(adj2igraph(getRefit(se_org)), vertex.color=dtype_org+1, vertex.size=9)

dtype_con <- c(rep(1,ntaxa(phylo_fungi_con)), rep(2,ntaxa(phylo_bacteria_con)))
plot(adj2igraph(getRefit(se_con)), vertex.color=dtype_con+1, vertex.size=9)

assoMat_org <- symBeta(SpiecEasi::getOptBeta(se_org),
                       mode = "ave") %>% as.matrix()

assoMat_con <- symBeta(SpiecEasi::getOptBeta(se_con),
                       mode = "ave") %>% as.matrix()
