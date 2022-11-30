# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

library(tidyverse)

## load data ---------------------------
# define version and load homogenized table
file.version <- "040722"
collecTRI_homogenized <- read_csv(file.path("output", file.version, "01_homogenized_table", "PMID_CollecTRI.csv"))

# load TF prior knowledge
tf.keywords <- read_csv("data/tf_annotation_activators_repressors.csv")

all.keywords <- readxl::read_excel("data/CollecTRI_TF-role_new_draft_211122.xlsx") %>%
  rename("TF" = "TFC2_Associated.Gene.Name",
         "strict" = "STRICT_agreement (GO/UniProt-StructureFunction)",
         "relaxed_KB" = "RELAXED - KB (STRICT_add aggregated GO/UniProt)",
         "relaxed_SF" =  "RELAXED - StructureFunction (STRICT_add KRAB/Soto)",
         "relaxed" = "RELAXED_KB OR StructureFunction")


## Create homogenized table ---------------------------
assign_sign <- function(homogenized.table,
                        tf.keywords = NULL, #needs to be provided if use.keywords is set to TRUE
                        use.keywords = T,
                        use.TFregulon = T,
                        unknown.to.pos = T,
                        keywords = "key_GO"){

  # select sign based on more number of PMIDs
  network.signed <- homogenized.table %>%
    mutate(sign = case_when(
      activation > repression ~ 1,
      repression > activation ~ -1,
      activation == repression ~ 0
    ))

  network.signed$TF <- unlist(lapply(str_split(network.signed$TF.TG, ":"), `[[`, 1))
  network.signed <- network.signed %>%
    mutate(decision = case_when(
    sign != 0 ~ "PMID"
  ))

  print("Interactions PMID:")
  print(table(network.signed$sign))

  # for unknown interactions use keywords to select sign
  if (use.keywords){
    if (keywords == "key_GO"){
      tf.keywords.collapsed <- data.frame(TF = tf.keywords$Transcription.Factor..Associated.Gene.Name.,
                                          activator_evidence = rowSums(tf.keywords[, str_detect(colnames(tf.keywords), "Activator")]),
                                          repressor_evidence = rowSums(tf.keywords[, str_detect(colnames(tf.keywords), "Repressor") |
                                                                                     str_detect(colnames(tf.keywords), "KRAB")])) %>%
        mutate(type = case_when(
          activator_evidence > repressor_evidence ~ "activator",
          repressor_evidence > activator_evidence ~ "repressor",
          activator_evidence == repressor_evidence ~ "unknown"
        ))

      network.signed <- left_join(network.signed, tf.keywords.collapsed %>% dplyr::select(TF, type), by = "TF")

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "activator" ~ 1,
          type == "repressor" ~ -1,
          type == "unknown" ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(decision = case_when(
          decision == "PMID" ~ "PMID",
          sign != 0 & is.na(decision) ~ "keywords"
        ))

      print("Interactions PMID + keywords:")
      print(table(network.signed$sign))
    } else if (keywords == "strict"){
      network.signed <- left_join(network.signed, tf.keywords %>% dplyr::select(TF, strict), by = "TF") %>%
        rename("type" = "strict")

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "Act" ~ 1,
          type == "Repr" ~ -1,
          is.na(type) ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(decision = case_when(
          decision == "PMID" ~ "PMID",
          sign != 0 & is.na(decision) ~ "keywords"
        ))

    } else if (keywords == "relaxed_KB"){
      network.signed <- left_join(network.signed, tf.keywords %>% dplyr::select(TF, relaxed_KB), by = "TF") %>%
        rename("type" = "relaxed_KB")

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "Act" ~ 1,
          type == "Repr" ~ -1,
          is.na(type) ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(decision = case_when(
          decision == "PMID" ~ "PMID",
          sign != 0 & is.na(decision) ~ "keywords"
        ))
    } else if (keywords == "relaxed_SF"){
      network.signed <- left_join(network.signed, tf.keywords %>% dplyr::select(TF, relaxed_SF), by = "TF") %>%
        rename("type" = "relaxed_SF")

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "Act" ~ 1,
          type == "Repr" ~ -1,
          is.na(type) ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(decision = case_when(
          decision == "PMID" ~ "PMID",
          sign != 0 & is.na(decision) ~ "keywords"
        ))
    } else if (keywords == "relaxed"){
      network.signed <- left_join(network.signed, tf.keywords %>% dplyr::select(TF, relaxed), by = "TF") %>%
        rename("type" = "relaxed")

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "Act" ~ 1,
          type == "Repr" ~ -1,
          is.na(type) ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(decision = case_when(
          decision == "PMID" ~ "PMID",
          sign != 0 & is.na(decision) ~ "keywords"
        ))
    }

  }

  # for other interactions use information from other interaction in regulon to select sign
  if (use.TFregulon){
    TF.regulon <- map_dfr(unique(network.signed$TF), function(TF_selected){
      TF.interactions <- network.signed %>% filter(TF == TF_selected)
      data.frame(TF = TF_selected,
                 perc_repressor = sum(TF.interactions$sign == -1)/nrow(TF.interactions),
                 perc_activator = sum(TF.interactions$sign == 1)/nrow(TF.interactions),
                 perc_unknown = sum(TF.interactions$sign == 0)/nrow(TF.interactions))

    })

    network.signed <- left_join(network.signed, TF.regulon, by = "TF")

    network.signed$sign[network.signed$perc_repressor > 0.5 & network.signed$sign == 0] <- -1
    network.signed$sign[network.signed$perc_activator > 0.5 & network.signed$sign == 0] <- 1

    network.signed <- network.signed %>%
      mutate(decision = case_when(
        decision == "PMID" ~ "PMID",
        decision == "keywords" ~ "keywords",
        sign != 0 & is.na(decision) ~ "regulon"
      ))

    print("Interactions PMID + keywords + TF regulon:")
    print(table(network.signed$sign))
  }

  if (unknown.to.pos){
    network.signed$sign[network.signed$sign == 0] <- 1
    network.signed$decision[is.na(network.signed$decision)] <- "unknown"
  }

  network.signed <- network.signed %>%
    filter(sign != 0)
  network.signed <- network.signed %>%
    dplyr::select(c(TF.TG, sign, PMID, TF.category, decision))

  print("Interactions PMID final:")
  print(table(network.signed$sign))

  return(network.signed)
}


## Construct and save networks ---------------------------
dir.create(file.path("output", file.version), showWarnings = FALSE)
dir.create(file.path("output", file.version, "02_signed_networks"), showWarnings = FALSE)

map(c("key_GO", "strict", "relaxed_KB", "relaxed_SF", "relaxed"), function(keyword){
  print(keyword)
  if (keyword == "key_GO") {
    signed_network <- assign_sign(homogenized.table = collecTRI_homogenized,
                                  tf.keywords = tf.keywords,
                                  keywords = keyword)
  }  else {
    signed_network <- assign_sign(homogenized.table = collecTRI_homogenized,
                                  tf.keywords = all.keywords,
                                  keywords = keyword)
    }
    signed_network <- signed_network %>%
      mutate(source = map_chr(str_split(TF.TG, ":"), 1),
             target = map_chr(str_split(TF.TG, ":"), 2))%>%
      rename("weight" = "sign") %>%
      dplyr::select(c(source, target, weight, PMID, TF.category, decision))

    write_csv(signed_network, file.path("output", file.version, "02_signed_networks", paste0(keyword, "_signed_CollecTRI.csv")))
})



