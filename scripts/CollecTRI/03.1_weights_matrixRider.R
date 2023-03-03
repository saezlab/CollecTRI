# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will try MatrixRider to calculate
#' binding affinities

library(MatrixRider)
library(TFBSTools)
library(Biostrings)
library(tidyverse)
library('org.Hs.eg.db')
source("scripts/helper/weights_functions.R")

## Load data ---------------------------
# define version and load signed network
output.folder <- "output"
raw.file <- "output/CollecTRI/CollecTRI.csv"

GRN <- read_csv(raw.file)

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

# load position matrix of TF
motif_list <- map(unique(GRN$source), function(source_i){
  motif <- MotifDb::MotifDb %>%
      MotifDb::query(andStrings=c(source_i, "hsapiens"))

  if(length(motif) != 0) {
    motif %>% universalmotif::convert_motifs(class = "TFBSTools-PFMatrix") %>%
      .[[1]]
  }
  else {
    NULL
  }
})

names(motif_list) <- unique(GRN$source)
pwm_tfs <- motif_list[!map_lgl(motif_list, is.null)]

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
dir.create(file.path(output.folder, "weighted_networks"), showWarnings = FALSE)
dir.create(file.path(output.folder, "weighted_networks", "raw"), showWarnings = FALSE)
write_csv(GRN_prob, file.path(output.folder, "weighted_networks", "raw", paste0("matrixRider_", prom_length, "bp.csv")))
