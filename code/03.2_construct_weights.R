# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will try MatrixRider to calculate
#' binding affinities

library(MatrixRider)
library(JASPAR2020)
library(TFBSTools)
library(Biostrings)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library('org.Hs.eg.db')

## Load data ---------------------------
GRN <- read.csv("output/signed_CollecTRI_dbTF_040722.csv") %>%
  mutate(TF = map_chr(str_split(TF.TG, ":"), 1),
         TG = map_chr(str_split(TF.TG, ":"), 2))


## Get promoter sequence
symbols <- unique(GRN$TG)
ids <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'ALIAS')
ids <- ids[!is.na(ids)]
transcriptCoordsByGene.GRangesList <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")

prom_sequence <- map(unique(ids), function(id){
  getPromSeq <- function(id){
    res <- getPromoterSeq (transcriptCoordsByGene.GRangesList[id],
                           Hsapiens, upstream=1000, downstream=0)
    return(res)
  }

  getPromSeq_poss <- possibly(.f = getPromSeq, otherwise = NA)
  getPromSeq_poss(id)

})
names(prom_sequence) <- unique(ids)
prom_sequence <- prom_sequence[!is.na(prom_sequence)]
sum(map_lgl(prom_sequence, function(x){class(x) == "DNAStringSetList"}))
sum(map_lgl(prom_sequence, function(x){!class(x) == "DNAStringSetList"}))


pwm_tfs <- map(unique(GRN$TF), function(tf){
  PWM_poss <- possibly(.f = getMatrixByName, otherwise = NA)

  tf_m <- paste0(substring(tf, 1, 1),
                 tolower(substring(tf, 2, nchar(tf))))


  pwm_h <- PWM_poss(JASPAR2020, tf)
  pwm_m <- PWM_poss(JASPAR2020, tf_m)

  if (class(pwm_h) == "PFMatrix"){
    pwm_h
  } else if (class(pwm_m) == "PFMatrix"){
    pwm_m
  } else {
    NA
  }
})
names(pwm_tfs) <- unique(GRN$TF)
pwm_tfs <- pwm_tfs[!is.na(pwm_tfs)]
sum(map_lgl(pwm_tfs, is.na))
sum(!map_lgl(pwm_tfs, is.na))


## calculate probabilities ---------------------------

GRN_prob <- map_dfr(1:nrow(GRN), function(j){
  print(j)
  x <- GRN[j,]
  if (any(names(ids) == x$TG)){
    gene_id <- ids[names(ids) == x$TG]
  } else {
    gene_id <- c("no gene id")
  }

  if (x$TF %in% names(pwm_tfs) & gene_id %in% names(prom_sequence)){
    pfm <- pwm_tfs[[x$TF]]
    bind_aff <- map_dbl(1:length(prom_sequence[[gene_id]][[1]]), function(i){
      sequ <- prom_sequence[[gene_id]][[1]][i]
      log(getSeqOccupancy(sequ[[1]], pfm, 0))
    })

    data.frame(TF = x$TF,
               TG = x$TG,
               bind_aff.mean = mean(bind_aff),
               bind_aff.sd = sd(bind_aff))
  } else {
    data.frame(TF = x$TF,
               TG = x$TG,
               bind_aff.mean = NA,
               bind_aff.sd = NA)
  }
})

GRN_prob$bind_aff.mean[is.infinite(GRN_prob$bind_aff.mean)] <- NA

saveRDS(GRN_prob, "output/binding_affinities_matrixRider.rds")

GRN_prob <- readRDS("output/binding_affinities_matrixRider.rds")

TF_weights <- GRN_prob %>%
  group_by(TF) %>%
  summarise(mean_weight = mean(bind_aff.mean, na.rm = T)) %>%
  arrange(desc(mean_weight)) %>%
  mutate(BindingMotif = case_when(
    is.na(mean_weight) ~ "Missing motif",
    !is.na(mean_weight) ~ "Motif",
  ))

# Number of TFs without motif
motif_p <- ggplot(data=TF_weights, aes(x=BindingMotif)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)

# Number of remaining TG without sequence
sequence_p <- GRN_prob %>%
  filter(TF %in% TF_weights$TF[TF_weights$BindingMotif == "Motif"]) %>%
  group_by(TG) %>%
  summarise(mean_weight = mean(bind_aff.mean, na.rm = T)) %>%
  mutate(PromoterSequence = case_when(
    is.na(mean_weight) ~ "Missing sequence",
    !is.na(mean_weight) ~ "Promoter sequence",
  )) %>%
  ggplot(aes(x=PromoterSequence)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)

# Total edges we cannot calculate weights for
edges_p <- GRN_prob %>%
  mutate(Missing = case_when(
    is.na(bind_aff.mean) ~ "Missing weight",
    !is.na(bind_aff.mean) ~ "Weight",
  )) %>%
  ggplot(aes(x=Missing)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)

# Distribution of weights
dist_weight_p <- ggplot(GRN_prob,aes(x = bind_aff.mean)) +
  geom_density() +
  geom_vline(xintercept = 1.5,
             color = "red")
dist_weight_tf_p <- ggplot(TF_weights,aes(x = mean_weight)) + geom_density()

GRN_prob_filter <- GRN_prob %>%
  filter(!is.na(bind_aff.mean))

# Based on distribution remove noise (cut off 1.5)
cut.off <- 1.5

GRN_prob_filter_nr <- GRN_prob_filter %>%
  mutate(bind_aff.mean = case_when(
    bind_aff.mean >= cut.off ~ bind_aff.mean,
    bind_aff.mean < cut.off ~0
  ))

# Number of edges considered noise
noise_p <- GRN_prob_filter_nr %>%
  mutate(Noise = case_when(
    bind_aff.mean == 0 ~ "Noise",
    bind_aff.mean > 0 ~ "No noise",
  )) %>%
  ggplot(aes(x=Noise)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1)

#
tmp <- GRN_prob_filter_nr %>%
  mutate(Noise = case_when(
    bind_aff.mean == 0 ~ 1,
    bind_aff.mean > 0 ~ 0,
  )) %>%
  group_by(TF) %>%
  summarise(n_noise = sum(Noise), n = n()) %>%
  mutate(perc_noise = n_noise/n)
ggplot(tmp, aes(x=perc_noise)) +
  geom_histogram(binwidth=0.1)

## Scaling
GRN_prob_filter_nr$bind_aff.mean[GRN_prob_filter_nr$bind_aff.mean == 0] <- NA
GRN_prob_filter2 <- GRN_prob_filter_nr %>% filter(!is.na(bind_aff.mean))

GRN_scaled <- map_dfr(unique(GRN_prob_filter2$TF), function(tf){
  GRN_tf <- GRN_prob_filter2 %>% dplyr::filter(TF == tf)

  GRN_tf_tmp <- GRN_tf %>% mutate(scaled_prob = (GRN_tf$bind_aff.mean) / (max(GRN_tf$bind_aff.mean, na.rm = T)))
  GRN_tf_tmp %>% select(TF, TG, scaled_prob)
})

# Distribution of scaledweights
dist_weight_scaled_p <- ggplot(GRN_scaled,aes(x = scaled_prob)) +
  geom_density()

TF_weights_scaled <- GRN_scaled %>%
  group_by(TF) %>%
  summarise(mean_weight = mean(scaled_prob, na.rm = T)) %>%
  arrange(desc(mean_weight))


dist_weight_scaled_tf_p <- ggplot(TF_weights_scaled,aes(x = mean_weight)) + geom_density()

## Make plot to show how much we loose by each step + how much per TF

