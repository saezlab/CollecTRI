# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will try MatrixRider to calculate
#' binding affinities

library(MatrixRider)
library(TFBSTools)
library(Biostrings)
library(tidyverse)
library('org.Hs.eg.db')
source("code/weights_functions.R")

## Load data ---------------------------
# define version and load signed network
file.version <- "040722"
collecTRI_homogenized <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))

GRN <- collecTRI_homogenized %>%
  filter(TF.category == "DbTF")


# construct or load promoter sequence of target genes
prom_length <- "10000"
dir.create(file.path("data", "promoter_sequence"), showWarnings = FALSE)

if (length(list.files(file.path("data", "promoter_sequence"),
                      pattern = paste0("prom_sequence_", prom_length, ".rds"), )) == 0){
  prom_sequence <- getPromSeq(genes = unique(GRN$target), n_upstream = as.numeric(prom_length))
  saveRDS(prom_sequence, file.path("data", "promoter_sequence", paste0("prom_sequence_", prom_length, ".rds")))
} else {
  prom_sequence <- readRDS(file.path("data", "promoter_sequence", paste0("prom_sequence_", prom_length, ".rds")))
}

pwm_tfs <- getPWM_motifDB(unique(GRN$source), class = "TFBSTools-PFMatrix")

# Link Gene symbol to gene ID (used for promoters)
symbols <- unique(GRN$target)
ids <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ALIAS')
ids <- ids[!is.na(ids)]

# Filter GRN for TFs and genes with available promoter sequence and binding motif
GRN <- left_join(GRN, ids %>%
            as.data.frame() %>%
            rownames_to_column("target")) %>%
  dplyr::rename("id" = ".")


GRN_filtered <- GRN %>%
  filter(source %in% names(pwm_tfs) & id %in% names(prom_sequence)) %>%
  filter(!is.na(id))

dim(GRN_filtered)

## calculate probabilities ---------------------------
GRN_prob <- map_dfr(1:nrow(GRN_filtered), function(j){
  print(j)
  x <- GRN_filtered[j,]

  pfm <- pwm_tfs[[x$source]]
  motif_length <- pfm@profileMatrix %>% ncol()
  sequ <- prom_sequence[[x$id]][[1]]
  bind_aff <- log(MatrixRider::getSeqOccupancy(sequ, pfm, 0), base = 10)

  data.frame(source = x$source,
             target = x$target,
             bind_aff = bind_aff,
             motif_length = motif_length)

})

## save results ---------------------------
saveRDS(GRN_prob, file.path("output", file.version, "03_weighting_strategies", paste0("matrixRider_res_", prom_length, "bp.rds")))
