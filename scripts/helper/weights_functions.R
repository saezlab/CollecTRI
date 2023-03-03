getPromSeq <- function(genes, n_upstream = 10000, n_downstream = 100, convert_to_ids = T){
  # convert gene names to ids
  if (convert_to_ids){
    ids <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'ALIAS')
    genes <- ids[!is.na(ids)]
  }

  # Transcript coordinated by gene
  transcriptCoordsByGene.GRangesList <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")

  prom_sequence <- map(genes, function(id){
    getPromSeq <- function(id){
      tscripts_gr <- stack(transcriptCoordsByGene.GRangesList[id], "gene_id")
      res <- getPromoterSeq (tscripts_gr[which.max(tscripts_gr@ranges@start)],
                             Hsapiens, upstream=n_upstream, downstream=n_downstream)
      return(res)
    }

    getPromSeq_poss <- possibly(.f = getPromSeq, otherwise = NA)
    getPromSeq_poss(id)

  })

  names(prom_sequence) <- unique(ids)
  prom_sequence_raw <- prom_sequence
  prom_sequence[!is.na(prom_sequence)]

}

getPWM_motifDB <- function(TFs, class = "universalmotif-universalmotif"){
  motifs <- map(TFs, function(TF){
    if (length(MotifDb::MotifDb %>%
              MotifDb::query(andStrings=c(TF, "hsapiens"))) == 0){
      NULL
    } else {
      motif <- MotifDb::MotifDb %>%
        MotifDb::query(andStrings=c(TF, "hsapiens")) %>%
        universalmotif::convert_motifs(class = class) %>%
        .[[1]]
      motif
    }
  })
  names(motifs) <- TFs
  return(motifs[!map_lgl(motifs, is.null)])
}
