library(tidyverse)
source("code/constructNetwork_function.R")

input_path <- 'data/homogenized_ressource_ExTRI_DbTF_v2.csv'
output_folder <- "data/networks_dbTF_v2/"

#### NTNU v2.0, ExTRI ####
# Construct different networks for ExTri data
# First we load the GRN curation tables. GRNcuration_ExTRI_tot just contains the information for the total number of evidence but not the information for each resource separately. Here the duplicated pubmed IDs between resources are just count once.
# The GRNcuration_ExTRI contains the information per resource in addition to the sum. If different ressources take their information from the same pubmed ID this is counted multiple times. The advantage is that we can scale the information for each resource.

# GRNcuration_ExTRI_tot <- read.table('data/evidence_counts_ExTRI_pairs.tsv', sep = "\t", header=FALSE) %>%
#   rename(TF.TG = V1, totPositive = V2, totNeg = V3, totUnknown = V4) %>%
#   column_to_rownames(var = "TF.TG")

GRNcuration_ExTRI <- read.csv(input_path) %>%
  column_to_rownames(var = "TF.TG")

# The networks to construct are then collected in a data.frame. I used different data.frames to generate different collections I would like to test.

networks_sign <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                           curated_counts = "scaled",
                           mor = TRUE,
                           mor_filter = 5:10/10,
                           weight = "evidence",
                           rm_no_sign_rows = TRUE)

networks_comparison <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                             curated_counts = "scaled",
                             mor = c(FALSE, TRUE, FALSE, TRUE),
                             mor_filter = 0.9,
                             weight = c("none", "none", "evidence", "evidence"),
                             rm_no_sign_rows = TRUE)

networks_sign_unrestricted <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                        curated_counts = "scaled",
                        mor = TRUE,
                        mor_filter = 5:10/10,
                        weight = "evidence",
                        rm_no_sign_rows = FALSE)

networks_comparison_unrestricted <- tibble(GRNcuration = list(GRNcuration_ExTRI),
                              curated_counts = "scaled",
                              mor = c(FALSE, TRUE, FALSE, TRUE),
                              mor_filter = 0.9,
                              weight = c("none", "none", "evidence", "evidence"),
                              rm_no_sign_rows = FALSE)


network_collection <- list(ExTRI_sign = networks_sign,
                           ExTRI_comp = networks_comparison,
                           ExTRI_sign_unrestricted = networks_sign_unrestricted,
                           ExTRI_comp_unrestricted = networks_comparison_unrestricted)

for (network in names(network_collection)){
    network_collection[[network]] <- network_collection[[network]] %>%
      add_column(path = paste0(output_folder, paste(network, network_collection[[network]]$curated_counts,
                                                    network_collection[[network]]$mor, network_collection[[network]]$mor_filter,
                                                    network_collection[[network]]$weight, network_collection[[network]]$rm_no_sign_rows, sep = "_"), ".rds"))
}

saveRDS(network_collection, paste0(output_folder, "network_collection_v2_dbTF.rds"))

# construct networks from the data.frames. the network within the constructNetwork function needs to be changed manually at the moment.
network_size <- map(names(network_collection), function(network){
  collection <- network_collection[[network]]
  res <- pmap(collection, function(GRNcuration, curated_counts, mor, mor_filter, weight, rm_no_sign_rows, path){
    GRN <- constructNetwork(GRNcuration = GRNcuration,
                            curated_counts = curated_counts,
                            mor = mor, mor_filter = mor_filter,
                            weight = weight,
                            rm_no_sign_rows = rm_no_sign_rows) %>% as_tibble()
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
network_size <- network_size %>% mutate(network = str_remove(str_remove(network_size$network, ".rds"), output_folder)) %>%
  add_row(network = "Dorothea_A", edges = "5378", TF = "94")

write_csv(network_size, paste0(output_folder, "network_size.csv"))
