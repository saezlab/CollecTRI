library(tidyverse)
source("code/constructNetwork_function.R")

# Construct different networks from curated NTNU information
# load curated information
GRNcuration <- read.csv('data/homogenized_resource.csv') %>%
  column_to_rownames(var = "TF.TG")

# Construct GRNs
GRN_agnostic <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = FALSE, weight = "none") %>% as_tibble()
saveRDS(GRN_agnostic, "data/GRN_agnostic.rds")

GRN_weighted_evidence <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = FALSE, weight = "evidence") %>% as_tibble()
saveRDS(GRN_weighted_evidence, "data/GRN_weighted_evidence.rds")

GRN_agnostic_unrestricted <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = FALSE, weight = "none", rm_bidirectional_rows = FALSE) %>% as_tibble()
saveRDS(GRN_agnostic_unrestricted, "data/GRN_agnostic_unrestricted.rds")

GRN_weighted_evidence_unrestricted <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = FALSE, weight = "evidence", rm_bidirectional_rows = FALSE) %>% as_tibble()
saveRDS(GRN_weighted_evidence_unrestricted, "data/GRN_weighted_evidence_unrestricted.rds")

GRN_signed <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, weight = "none") %>% as_tibble()
saveRDS(GRN_signed, "data/GRN_signed.rds")

GRN_signed <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, mor_filter = 0.8, weight = "none") %>% as_tibble()
saveRDS(GRN_signed, "data/GRN_signed_0.8.rds")

GRN_signed_weighted_evidence <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, weight = "evidence") %>% as_tibble()
saveRDS(GRN_signed_weighted_evidence, "data/GRN_signed_weighted_evidence.rds")

GRN_signed_weighted_evidence <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, mor_filter = 0.8, weight = "evidence") %>% as_tibble()
saveRDS(GRN_signed_weighted_evidence, "data/GRN_signed_weighted_evidence_0.8.rds")

GRN_signed_weighted_evidence_sign <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, weight = "evidence_sign") %>% as_tibble()
saveRDS(GRN_signed_weighted_evidence_sign, "data/GRN_signed_weighted_evidence_sign.rds")

GRN_signed_weighted_evidence_sign <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, mor_filter = 0.8, weight = "evidence_sign") %>% as_tibble()
saveRDS(GRN_signed_weighted_evidence_sign, "data/GRN_signed_weighted_evidence_sign_0.8.rds")

GRN_signed_weighted_discrepancy_sign <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, weight = "discrepancy_sign") %>% as_tibble()
saveRDS(GRN_signed_weighted_discrepancy_sign, "data/GRN_signed_weighted_discrepancy_sign.rds")

GRN_signed_weighted_discrepancy_sign <- constructNetwork(GRNcuration = GRNcuration, curated_counts = 'scaled', mor = TRUE, mor_filter = 0.8, weight = "discrepancy_sign") %>% as_tibble()
saveRDS(GRN_signed_weighted_discrepancy_sign, "data/GRN_signed_weighted_discrepancy_sign_0.8.rds")


# Construct different networks from dorothea
# load curated information
dorothea <- readRDS("data/dorothea_filtered.rds")
dorothea <- dorothea %>% filter(confidence == "A")
NTNU <- readRDS("data/GRN_weighted_evidence.rds")

merged_network <- merge(
  x=dorothea,
  y=NTNU,
  by.x=c("source","target"),
  by.y=c("source","target"))

# Construct GRNs
dorothea_signed <- merged_network %>% select(c(source, target, confidence.x, mor.x, likelihood.x)) %>% rename(confidence = confidence.x, likelihood = likelihood.x, mor = mor.x)
saveRDS(dorothea_signed, "data/dorothea_signed.rds")

dorothea_agnostic <- merged_network %>% select(c(source, target, confidence.x, mor.y, likelihood.x)) %>% rename(confidence = confidence.x, likelihood = likelihood.x, mor = mor.y)
saveRDS(dorothea_agnostic, "data/dorothea_agnostic.rds")

dorothea_agnostic_weighted <- merged_network %>% select(c(source, target, confidence.x, mor.y, likelihood.y)) %>% rename(confidence = confidence.x, likelihood = likelihood.y, mor = mor.y)
saveRDS(dorothea_agnostic_weighted, "data/dorothea_agnostic_weighted.rds")

dorothea_signed_weighted <- merged_network %>% select(c(source, target, confidence.x, mor.x, likelihood.y)) %>% rename(confidence = confidence.x, likelihood = likelihood.y, mor = mor.x)
saveRDS(dorothea_signed_weighted, "data/dorothea_signed_weighted.rds")
