library(tidyverse)
library(decoupleR)

setwd("/net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR")

# final networks
input_path_dorothea <- c(file.path('data', "dorothea", "dorothea_A_new.rds"),
                         file.path('data', "dorothea", "dorothea_ABC_new.rds"))
input_path_v1 <- file.path('data', 'networks_v1', 'network_collection_v1.rds')
input_path_v2 <- file.path('data', 'networks_v2', 'network_collection_v2.rds')
input_path_v2_dbTF <- file.path('data', 'networks_dbTF_v2', 'network_collection_v2_dbTF.rds')

output_path <- file.path('output')

network_collection_v1 <- readRDS(input_path_v1)$ExTRI_comp
network_collection_v2 <- readRDS(input_path_v2)$ExTRI_comp
network_collection_v2_dbTF <- readRDS(input_path_v2_dbTF)$ExTRI_comp

networks <- rbind(network_collection_v1, network_collection_v2, network_collection_v2_dbTF)
path_networks <- c(input_path_dorothea, networks$path)

final_networks <- map(path_networks, readRDS)
names(final_networks) <- c("Dorothea A", "Dorothea ABC",
                           "NTNU.1", "NTNU.1 s", "NTNU.1 w", "NTNU.1 s+w",
                           "NTNU.2", "NTNU.2 s", "NTNU.2 w", "NTNU.2 s+w",
                           "NTNU.2 dbTF", "NTNU.2 dbTF s", "NTNU.2 dbTF w", "NTNU.2 dbTF s+w")
bmeta_knockTF <- readRDS(file.path('data',"bench", "knockTF_meta.rds"))
bexpr_knockTF <-  readRDS(file.path('data',"bench", "knockTF_exp.rds"))

decoupler_act <- map(final_networks, function(net){
  net <- decoupleR::intersect_regulons(bexpr_knockTF, net, .source = "source", .target = "target", minsize = 5)
  decouple_res <- decouple(bexpr_knockTF, net)
  decouple_res <- decouple_res %>% filter(statistic == "consensus")

  ntargets <- table(net$source) %>% as.data.frame() %>% rename("source" = Var1)

  full_join(decouple_res, ntargets)
})

saveRDS(output_path, "/decoupler_act.rds")
