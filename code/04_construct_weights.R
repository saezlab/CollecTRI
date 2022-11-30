# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will construct weighted networks from
#' 1) matrixRider
#' 2) FIMO
#' 3) RcisTarget


library(tidyverse)

## Load data ---------------------------
# define version and load signed network
file.version <- "040722"
collecTRI_homogenized <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))

GRN <- collecTRI_homogenized %>%
  filter(TF.category == "DbTF")

matrixRider.res <-  readRDS(file.path("output", file.version, "03_weighting_strategies", "matrixRider_res_1000bp.rds"))
FIMO.res.1 <- readRDS(file.path("output", file.version, "03_weighting_strategies", "FIMO_res_1000bp.rds"))
FIMO.res.10 <- readRDS(file.path("output", file.version, "03_weighting_strategies", "FIMO_res_10000bp.rds"))

dir.create(file.path("output", file.version, "04_weighted_networks"), showWarnings = FALSE)

set.seed(123)

## matrixRider ---------------------------
# use weights raw
GRN_matRid <- matrixRider.res %>%
  mutate(TF.TG = paste(source, target, sep = ".")) %>%
  filter(!duplicated(TF.TG)) %>%
  mutate(weight = bind_aff) %>%
  filter(!is.na(weight)) %>%
  filter(weight > 0)

# use weights normalised by gene
GRN_matRid_scaled_gene <- map_dfr(unique(GRN_matRid$target), function(tg){
  GRN_tf <- GRN_matRid %>% dplyr::filter(target == tg)

  GRN_tf_tmp <- GRN_tf %>% mutate(weight = (GRN_tf$weight) / (max(GRN_tf$weight, na.rm = T)))
  GRN_tf_tmp %>% dplyr::select(source, target, weight)
})


# use weights normalised by TF
GRN_matRid_scaled_TF <- map_dfr(unique(GRN_matRid$source), function(tf){
  GRN_tf <- GRN_matRid %>% dplyr::filter(source == tf)

  GRN_tf_tmp <- GRN_tf %>% mutate(weight = (GRN_tf$weight) / (max(GRN_tf$weight, na.rm = T)))
  GRN_tf_tmp %>% dplyr::select(source, target, weight)
})


# use weights for filtering
quan_10 <- quantile(GRN_matRid$weight, probs = c(0.10), na.rm = T)
quan_20 <- quantile(GRN_matRid$weight, probs = c(0.20), na.rm = T)
quan_30 <- quantile(GRN_matRid$weight, probs = c(0.30), na.rm = T)

quantiles <- c(quan_10, quan_20, quan_30)
matRid_filtered <- map(names(quantiles), function(q_n){
  q <- quantiles[q_n]
  GRN_quant <- GRN_matRid %>%
    filter(weight >= q) %>%
    mutate(weight = 1)
})

GRN_random_10 <- GRN_matRid[sample(nrow(GRN_matRid), nrow(matRid_filtered[[1]])), ] %>%
  mutate(weight = 1)
GRN_random_20 <- GRN_matRid[sample(nrow(GRN_matRid), nrow(matRid_filtered[[2]])), ] %>%
  mutate(weight = 1)
GRN_random_30 <- GRN_matRid[sample(nrow(GRN_matRid), nrow(matRid_filtered[[3]])), ] %>%
  mutate(weight = 1)

# Save results
dir.create(file.path("output", file.version, "04_weighted_networks", "matrixRider"), showWarnings = FALSE)
write_csv(GRN_matRid %>% mutate(weight = 1), file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_res.csv"))
write_csv(GRN_matRid, file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_matrixRider_raw.csv"))
write_csv(GRN_matRid_scaled_gene, file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_matrixRider_scaled_gene.csv"))
write_csv(GRN_matRid_scaled_TF, file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_matrixRider_scaled_tf.csv"))

write_csv(matRid_filtered[[1]], file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_matrixRider_filtered_10.csv"))
write_csv(matRid_filtered[[2]], file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_matrixRider_filtered_20.csv"))
write_csv(matRid_filtered[[3]], file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_matrixRider_filtered_30.csv"))
write_csv(GRN_random_10, file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_random_filtered_10.csv"))
write_csv(GRN_random_20, file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_random_filtered_20.csv"))
write_csv(GRN_random_30, file.path("output", file.version, "04_weighted_networks", "matrixRider", "collecTRI_random_filtered_30.csv"))


## FIMO ---------------------------
# use weights raw
FIMO.res.1 <- FIMO.res.1 %>%
  mutate(TF.TG = paste(motif_id, gene_name, sep = "."))

FIMO.res.10 <- FIMO.res.10 %>%
  mutate(TF.TG = paste(motif_id, gene_name, sep = "."))

FIMO.res <- full_join(FIMO.res.1, FIMO.res.10, by = "TF.TG")
FIMO.res <- FIMO.res %>%
  filter(!duplicated(TF.TG))


GRN_fimo <- GRN %>%
  mutate(TF.TG = paste(source, target, sep = ".")) %>%
  dplyr::left_join(FIMO.res %>% dplyr::select(TF.TG, score.y)) %>%
  mutate(weight = weight*score.y) %>%
  filter(!is.na(weight)) %>%
  filter(weight > 0)

# use weights normalised by gene
GRN_fimo_scaled_gene <- map_dfr(unique(GRN_fimo$target), function(tg){
  GRN_tf <- GRN_fimo %>% dplyr::filter(target == tg)

  GRN_tf_tmp <- GRN_tf %>% mutate(weight = (GRN_tf$weight) / (max(GRN_tf$weight, na.rm = T)))
  GRN_tf_tmp %>% dplyr::select(source, target, weight)
})


# use weights normalised by TF
GRN_fimo_scaled_TF <- map_dfr(unique(GRN_fimo$source), function(tf){
  GRN_tf <- GRN_fimo %>% dplyr::filter(source == tf)

  GRN_tf_tmp <- GRN_tf %>% mutate(weight = (GRN_tf$weight) / (max(GRN_tf$weight, na.rm = T)))
  GRN_tf_tmp %>% dplyr::select(source, target, weight)
})

# use weights for filtering (fdr 40%)
FP_FIMO.1 <- FIMO.res.1 %>%
  filter(pvalue > 0.0002) %>%
  mutate(TF.TG = paste(motif_id, gene_name, sep = ".")) %>%
  pull(TF.TG)

FP_FIMO.10 <- FIMO.res.10 %>%
  filter(pvalue > 0.00002) %>%
  mutate(TF.TG = paste(motif_id, gene_name, sep = ".")) %>%
  pull(TF.TG)

FP_FIMO <- intersect(FP_FIMO.1, FP_FIMO.10)

GRN_FIMO_filtered <- GRN  %>%
  mutate(TF.TG = paste(source, target, sep = ".")) %>%
  dplyr::filter(!(TF.TG %in% FP_FIMO))

GRN_random_filtered <- GRN[sample(nrow(GRN), nrow(GRN)-length(FP_FIMO)), ]

# Save results
dir.create(file.path("output", file.version, "04_weighted_networks", "FIMO"), showWarnings = FALSE)
write_csv(GRN_fimo %>% mutate(weight = 1), file.path("output", file.version, "04_weighted_networks", "FIMO", "collecTRI_res.csv"))
write_csv(GRN_fimo, file.path("output", file.version, "04_weighted_networks", "FIMO", "collecTRI_FIMO_raw.csv"))
write_csv(GRN_fimo_scaled_gene, file.path("output", file.version, "04_weighted_networks", "FIMO", "collecTRI_FIMO_scaled_gene.csv"))
write_csv(GRN_fimo_scaled_TF, file.path("output", file.version, "04_weighted_networks", "FIMO", "collecTRI_FIMO_scaled_tf.csv"))

write_csv(GRN_FIMO_filtered, file.path("output", file.version, "04_weighted_networks", "FIMO", "collecTRI_FIMO_filtered.csv"))
write_csv(GRN_random_filtered, file.path("output", file.version, "04_weighted_networks", "FIMO", "collecTRI_random_filtered.csv"))
