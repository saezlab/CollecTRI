# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will run FIMO to calculate
#' binding affinities

library(tidyverse)
library('org.Hs.eg.db')
library(GenomicRanges)
library(magrittr)
library(universalmotif)
library(memes)
source("code/weights_functions.R")

options(meme_bin = "/opt/local/bin/")

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
  prom_sequence <- getPromSeq(genes = unique(GRN$target), n_upstream = prom_length)
  saveRDS(prom_sequence, file.path("data", "promoter_sequence", paste0("prom_sequence_", prom_length, ".rds")))
} else {
  prom_sequence <- readRDS(file.path("data", "promoter_sequence", paste0("prom_sequence_", prom_length, ".rds")))
}

# Link Gene symbol to gene ID (used for promoters)
symbols <- unique(GRN$target)
ids <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ALIAS')
ids <- ids[!is.na(ids)]


## run FIMO ---------------------------
fimo_res <- map_dfr(unique(GRN$source), function(source_i){
  print(source_i)
  print(paste0(which(unique(GRN$source) == source_i), "/", length(unique(GRN$source))))
  targets <- GRN %>% filter(source == source_i) %>% pull(target)
  target_ids <- ids[names(ids) %in% targets]

  if(length(MotifDb::MotifDb %>%
            MotifDb::query(andStrings=c(source_i, "hsapiens"))) != 0) {
    motif <- MotifDb::MotifDb %>%
      MotifDb::query(andStrings=c(source_i, "hsapiens")) %>%
      universalmotif::convert_motifs() %>%
      .[[1]]
    motif["name"] <- source_i

    tmp <- prom_sequence[target_ids]
    idx <- map_lgl(tmp, is.null)
    tmp <- tmp[!idx]
    names(tmp) <- NULL
    sequences <- do.call(c,tmp)
    names(sequences) <- names(prom_sequence[target_ids])[!idx]

    fimo_results <- runFimo(sequences, motif, thresh = 1)

    map_ids <- target_ids %>%
      as.data.frame() %>%
      rownames_to_column("gene_name") %>%
      dplyr::rename("seqnames" = ".")

    fimo_results %>%
      data.frame() %>%
      group_by(seqnames) %>%
      dplyr::filter(pvalue == min(pvalue)) %>%
      dplyr::select(seqnames, motif_id, score, pvalue) %>%
      dplyr::left_join(map_ids)
  } else {
    data.frame(seqnames = "no_motif", motif_id = source_i, score = NA, gene_name = NA)
  }

})

## save results ---------------------------

saveRDS(fimo_res, file.path("output", file.version, "03_weighted_strategies", paste0("FIMO_res_", prom_length, "bp.rds")))
