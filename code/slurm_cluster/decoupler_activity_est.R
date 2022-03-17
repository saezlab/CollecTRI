library(tidyverse)
library(decoupleR)

setwd("/net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR")

# final networks
path_networks <- c(dorothea_A = "data/dorothea/dorothea_A_new.rds",
                   dorothea_ABC = "data/dorothea/dorothea_ABC_new.rds",
                   v1_weighted = "data/networks_v1/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v1_weighted_signed = "data/networks_v1/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v1 = "data/networks_v1/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                   v1_signed = "data/networks_v1/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                   v2_weighted = "data/networks_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v2_weighted_signed = "data/networks_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v2 = "data/networks_v2/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                   v2_signed = "data/networks_v2/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                   v2_dbTF_weighted = "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v2_dbTF_weighted_signed = "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v2_dbTF = "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                   v2_dbTF_signed = "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds")

final_networks <- map(path_networks, readRDS)
names(final_networks) <- names(decoupler_act) <- c("Dorothea A", "Dorothea ABC",
                                                   "NTNU.1 w", "NTNU.1 w+s", "NTNU.1", "NTNU.1 s",
                                                   "NTNU.2 w", "NTNU.2 w+s", "NTNU.2", "NTNU.2 s",
                                                   "NTNU.2 dbTF w", "NTNU.2 dbTF w+s", "NTNU.2 dbTF", "NTNU.2 dbTF s")
bmeta_knockTF <- readRDS(file.path('data',"bench", "knockTF_meta.rds"))
bexpr_knockTF <-  readRDS(file.path('data',"bench", "knockTF_exp.rds")) %>% as.matrix()

decoupler_act <- map(final_networks, function(net){
  net <- decoupleR::intersect_regulons(bexpr_knockTF, net, .source = "source", .target = "target", minsize = 5)
  decouple_res <- decouple(bexpr_knockTF[,1:2], net)
  decouple_res <- decouple_res %>% filter(statistic == "consensus")

  ntargets <- table(net$source) %>% as.data.frame() %>% rename("source" = Var1)

  full_join(decouple_res, ntargets)
})

saveRDS(decoupler_act, "decoupler_act.rds")
