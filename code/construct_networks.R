library(tidyverse)
source("code/constructNetwork_function.R")

#### NTNU v1.0 ####
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

#### Dorothea ####
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

#### NTNU v2.0, ExTRI ####
# Construct different networks for ExTri data
# First we load the GRN curation tables. GRNcuration_ExTRI_tot just contains the information for the total number of evidence but not the information for each resource separately. Here the duplicated pubmed IDs between resources are just count once.
# The GRNcuration_ExTRI contains the information per resource in addition to the sum. If different ressources take their information from the same pubmed ID this is counted multiple times. The advantage is that we can scale the information for each resource.
GRNcuration_ExTRI_tot <- read.table('data/evidence_counts_ExTRI_pairs.tsv', sep = "\t", header=FALSE) %>%
  rename(TF.TG = V1, totPositive = V2, totNeg = V3, totUnknown = V4) %>%
  column_to_rownames(var = "TF.TG")

GRNcuration_ExTRI <- read.csv('data/homogenized_ressource_ExTRI.csv') %>%
  column_to_rownames(var = "TF.TG")

# To test the effect of ExTRI that information is removed from the table.
GRNcuration_new <- GRNcuration_ExTRI %>%
  mutate(totNeg = totNeg-ExTRI_Negative) %>%
  mutate(totPositive = totPositive-ExTRI_Positive) %>%
  mutate(totUnknown = totUnknown-ExTRI_Unknown) %>%
  select(-c(ExTRI_Negative, ExTRI_Positive, ExTRI_Unknown)) %>%
  filter(!totUnknown + totPositive + totNeg == 0)

# The networks to construct are then collected in a data.frame. I used different data.frames to generate different collections I would like to test.
networks_weights <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                          curated_counts = "scaled",
                          mor = TRUE,
                          mor_filter = NA,
                          weight = c("none", "evidence", "evidence_sign", "discrepancy_sign"),
                          rm_bidirectional_rows = TRUE)

networks_sign <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                           curated_counts = "scaled",
                           mor = TRUE,
                           mor_filter = 0:10/10,
                           weight = "evidence",
                           rm_bidirectional_rows = TRUE)

networks_tot <- tibble(GRNcuration = list(GRNcuration_ExTRI_tot),
                           curated_counts = "standard",
                           mor = c(FALSE, TRUE, FALSE, TRUE),
                           mor_filter = NA,
                           weight = c("none", "none", "evidence_tot", "evidence_tot"),
                           rm_bidirectional_rows = TRUE)

networks_comparison <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                             curated_counts = "scaled",
                             mor = c(FALSE, TRUE, FALSE, TRUE),
                             mor_filter = NA,
                             weight = c("none", "none", "evidence", "evidence"),
                             rm_bidirectional_rows = TRUE)

networks_comparison_unrestricted <- networks_comparison %>% mutate(rm_bidirectional_rows = FALSE)

network_collection <- list(ExTRI_weights = networks_weights,
                           ExTRI_sign = networks_sign,
                           ExTRI_tot_comp = networks_tot,
                           ExTRI_comp = networks_comparison,
                           ExTRI_comp_unrestricted = networks_comparison_unrestricted)

for (network in names(network_collection)){
    network_collection[[network]] <- network_collection[[network]] %>% add_column(path = paste0("data/", paste(network, network_collection[[network]]$curated_counts,
                                                                                                               network_collection[[network]]$mor, network_collection[[network]]$mor_filter,
                                                                                                               network_collection[[network]]$weight, network_collection[[network]]$rm_bidirectional_rows, sep = "_"), ".rds"))
}

saveRDS(network_collection, "data/network_collection.rds")

# construct networks from the data.frames. the network within the constructNetwork function needs to be changed manually at the moment.
network_size <- map(names(network_collection), function(network){
  collection <- network_collection[[network]]
  res <- pmap(collection, function(GRNcuration, curated_counts, mor, mor_filter, weight, rm_bidirectional_rows, path){
    GRN <- constructNetwork(GRNcuration = GRNcuration,
                            curated_counts = curated_counts,
                            mor = mor, mor_filter = mor_filter,
                            weight = weight,
                            rm_bidirectional_rows = rm_bidirectional_rows) %>% as_tibble()
    saveRDS(GRN, path)
    network_size <- c(network = path,
                      edges = length(GRN$source),
                      TF = length(unique(GRN$source)))
    return(network_size)
})
  res <- as.data.frame(do.call(rbind, res))
  return(res)
})

network_size <- as.data.frame(do.call(rbind, network_size))
network_size <- network_size %>% mutate(network = str_remove(str_remove(network_size$network, ".rds"), "data/")) %>%
  add_row(network = "Dorothea_A", edges = "5378", TF = "94")

write_csv(network_size, "data/network_size.csv")

