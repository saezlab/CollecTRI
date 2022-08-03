library(tidyverse)

collecTRI_homogenized <- read_csv("output/PMID_CollecTRI_040722.csv")
collecTRI_dbTF_homogenized <- read_csv("output/PMID_CollecTRI_dbTF_040722.csv")
tf.keywords <- read_csv("data/tf_annotation_activators_repressors.csv")

assign_sign <- function(homogenized.table,
                        tf.keywords = NULL, #needs to be provided if use.keywords is set to TRUE
                        use.keywords = T,
                        use.TFregulon = T,
                        unknown.to.pos = T){

  # select sign based on more number of PMIDs
  network.signed <- homogenized.table %>%
    mutate(sign = case_when(
      activation > repression ~ 1,
      repression > activation ~ -1,
      activation == repression ~ 0
    ))

  network.signed$TF <- unlist(lapply(str_split(network.signed$TF.TG, ":"), `[[`, 1))

  print("Interactions PMID:")
  print(table(network.signed$sign))

  # for unknown interactions use keywords to select sign
  if (use.keywords){
    tf.keywords.collapsed <- data.frame(TF = tf.keywords$Transcription.Factor..Associated.Gene.Name.,
                                        activator_evidence = rowSums(tf.keywords[, str_detect(colnames(tf.keywords), "Activator")]),
                                        repressor_evidence = rowSums(tf.keywords[, str_detect(colnames(tf.keywords), "Repressor") |
                                                                                   str_detect(colnames(tf.keywords), "KRAB")])) %>%
      mutate(type = case_when(
        activator_evidence > repressor_evidence ~ "activator",
        repressor_evidence > activator_evidence ~ "repressor",
        activator_evidence == repressor_evidence ~ "unknown"
      ))

    network.signed <- left_join(network.signed, tf.keywords.collapsed %>% select(TF, type), by = "TF")

    network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
      mutate(sign = case_when(
        type == "activator" ~ 1,
        type == "repressor" ~ -1,
        type == "unknown" ~ 0
      ))

    print("Interactions PMID + keywords:")
    print(table(network.signed$sign))
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

    print("Interactions PMID + keywords + TF regulon:")
    print(table(network.signed$sign))
  }

  if (unknown.to.pos){
    network.signed$sign[network.signed$sign == 0] <- 1
  }

  network.signed <- network.signed %>%
    filter(sign != 0)
  network.signed <- network.signed %>%
    select(c(TF.TG, sign, PMID))

  print("Interactions PMID final:")
  print(table(network.signed$sign))

  return(network.signed)
}


signed_network <- assign_sign(homogenized.table = collecTRI_homogenized,
                              tf.keywords = tf.keywords)

signed_dbTF_network <- assign_sign(homogenized.table = collecTRI_dbTF_homogenized,
                                   tf.keywords = tf.keywords)

write_csv(signed_network, "output/signed_CollecTRI_040722.csv")
write_csv(signed_dbTF_network, "output/signed_CollecTRI_dbTF_040722.csv")


