# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will run FIMO to calculate
#' binding affinities

library(tidyverse)
library(RcisTarget)

## Load data ---------------------------
# define version and load signed network
file.version <- "040722"
collecTRI_homogenized <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))

GRN <- collecTRI_homogenized %>%
  filter(TF.category == "DbTF")

# define gene list
geneLists <- map(unique(GRN$source), function(x){
  GRN %>% filter(source == x) %>% pull(target)})
names(geneLists) <- unique(GRN$source)

geneLists <- geneLists[map_lgl(geneLists, function(x){length(x) >= 5})]

# Load data bases
motifRankings <- importRankings(file.path("data", "hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"))
data(motifAnnotations_hgnc)

## run RcisTarget ---------------------------
motifEnrichmentTable_wGenes <- cisTarget(geneLists,
                                         motifRankings,
                                         motifAnnot=motifAnnotations)

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable_wGenes$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

numberTF_identified <- map_lgl(names(anotatedTfs), function(x){x %in% anotatedTfs[[x]]}) %>% sum()
# for only 127/601 TFs we were able to identify the motif of the TF to be enriched in the gene set
signifMotifNames <- motifEnrichmentTable_wGenes$motif[1:3]


## save results ---------------------------
saveRDS(motifEnrichmentTable_wGenes, file.path("output", file.version, "03_weighted_strategies", paste0("RcisTarget_res.rds")))
