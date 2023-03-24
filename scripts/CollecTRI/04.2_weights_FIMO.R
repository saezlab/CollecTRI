# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will run FIMO to calculate
#' binding affinities

library(tidyverse)
library('org.Hs.eg.db')
library(GenomicRanges)
library(magrittr)
#library(GenomicFeatures)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(universalmotif)
library(memes)
source("scripts/helper/weights_functions.R")

options(meme_bin = "/opt/local/bin/")

## Load data ---------------------------
# define version and load signed network
output.folder <- "output"
raw.file <- "output/CollecTRI/CollecTRI_GRN.csv"

GRN <- read_csv(raw.file)

# construct or load promoter sequence of target genes
prom_lengths <- c("1000", "10000") #number of base pairs
dir.create(file.path("data", "promoter_sequence"), showWarnings = FALSE)

for (prom_length in prom_lengths){
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

  # load position matrix of TF
  motif_list <- map(unique(GRN$source), function(source_i){
    motif <- MotifDb::MotifDb %>%
      MotifDb::query(andStrings=c(source_i, "hsapiens"))

    if(length(motif) != 0) {
      motif <- motif %>% universalmotif::convert_motifs() %>%
        .[[1]]
      motif["name"] <- source_i
      motif
    }
    else {
      NULL
    }
  })

  names(motif_list) <- unique(GRN$source)
  pwm_tfs <- motif_list[!map_lgl(motif_list, is.null)]

  # Filter GRN for TFs and genes with available promoter sequence and binding motif
  GRN <- left_join(GRN, ids %>%
                     as.data.frame() %>%
                     rownames_to_column("target")) %>%
    dplyr::rename("id" = ".")


  GRN_filtered <- GRN %>%
    filter(source %in% names(pwm_tfs) & id %in% names(prom_sequence)) %>%
    filter(!is.na(id))

  dim(GRN_filtered)

  ## run FIMO ---------------------------
  fimo_res <- map_dfr(unique(GRN_filtered$source), function(source_i){
    print(paste0(which(unique(GRN_filtered$source) == source_i), "/", length(unique(GRN_filtered$source))))

    target_ids <- GRN_filtered %>% filter(source == source_i) %>% pull(id)

    motif <- motif_list[[source_i]]

    tmp <- prom_sequence[target_ids]
    names(tmp) <- NULL
    sequences <- do.call(c,tmp)
    names(sequences) <- names(prom_sequence[target_ids])

    fimo_results <- runFimo(sequences, motif, thresh = 1)

    map_ids <- GRN_filtered %>%
      dplyr::filter(source == source_i) %>%
      dplyr::select(target, id) %>%
      dplyr::rename("seqnames" = "id")

    fimo_results %>%
      data.frame() %>%
      group_by(seqnames) %>%
      dplyr::filter(pvalue == min(pvalue)) %>%
      dplyr::select(seqnames, motif_id, score, pvalue) %>%
      dplyr::left_join(map_ids) %>%
      dplyr::rename("source" = "motif_id")
  })

  fimo_res <- fimo_res %>%
    mutate(edge = paste(source, target, sep = "."))
  fimo_res <- fimo_res[!duplicated(fimo_res$edge),]

  ## save results ---------------------------
  dir.create(file.path(output.folder, "weighted_networks"), showWarnings = FALSE)
  dir.create(file.path(output.folder, "weighted_networks", "raw"), showWarnings = FALSE)
  write_csv(fimo_res, file.path(output.folder, "weighted_networks", "raw", paste0("FIMO_", prom_length, "bp.csv")))

}
