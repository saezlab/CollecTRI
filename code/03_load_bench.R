library(tidyverse)
knockTF <- read.table("data/differential expression of genes in all datasets.txt",
                      sep = "\t", header = T)

path <- file.path('data','bench')
network_name <- "dorotheaA_new"
network_path <- "data/dorothea/dorothea_A_new.rds"


knockTF_exp <- knockTF %>%
  select(c(Sample_ID, Gene, Log2FC)) %>%
  pivot_wider(names_from = Sample_ID, values_from = Log2FC) %>%
  column_to_rownames("Gene") %>%
  replace(is.na(.), 0) %>%
  as.matrix() %>%
  replace(is.infinite(.), 0)

knockTF_meta <- unique(knockTF[,c('Sample_ID','TF')]) %>%
  mutate(sign = -1) %>%
  rename(c("id" = Sample_ID, "target" = TF))

saveRDS(knockTF_exp, file.path(path, "knockTF_exp.rds"))
saveRDS(knockTF_meta, file.path(path, "knockTF_meta.rds"))

# Filter by overlap of TFs
network <- readRDS(file.path(network_path))
tfs <- network$source
knockTF_meta <- dplyr::filter(knockTF_meta, target %in% tfs)
knockTF_exp <- knockTF_exp[, colnames(knockTF_exp) %in% knockTF_meta$id]

# Save
saveRDS(knockTF_exp, paste0(path, '/knockTF_expr_', network_name,'.rds'))
saveRDS(knockTF_meta, paste0(path, '/knockTF_meta_', network_name, '.rds'))

get_rna_data <- function(path){
  # Download
  expr_url <- 'https://zenodo.org/record/5645208/files/rna_expr.rds?download=1'
  meta_url <- 'https://zenodo.org/record/5645208/files/rna_meta.rds?download=1'
  download.file(expr_url, file.path(path, 'rna_expr.rds'))
  Sys.sleep(1)
  download.file(meta_url, file.path(path, 'rna_meta.rds'))
  Sys.sleep(1)

  # Filter by overlap of TFs
  rna_expr <- readRDS(file.path(path, 'rna_expr.rds'))
  rna_meta <- readRDS(file.path(path, 'rna_meta.rds'))
  network <- readRDS(file.path(network_path))
  tfs <- network$source
  rna_meta <- dplyr::filter(rna_meta, target %in% tfs)
  rna_expr <- rna_expr[, colnames(rna_expr) %in% rna_meta$id]

  # Save
  saveRDS(rna_expr, paste0(path, '/rna_expr_', network_name,'.rds'))
  saveRDS(rna_meta, paste0(path, '/rna_meta_', network_name, '.rds'))
}
get_rna_data(path)
