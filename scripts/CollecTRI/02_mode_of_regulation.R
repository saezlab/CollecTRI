#' In this script we assign a mode of regulation, either activation or
#' repression to each TF-gene interaction present. We use muliple sources
#' of information
#'
#' 1) Information from PMIDs
#' 2) TF classification into repressors or activators
#' 3) Other genes in the regulon

library(tidyverse)

## Load and prepare data ---------------------------
output.folder <- "output"
resource.file <- "output/aggregated_resources/CollecTRI_resources.csv"

collecTRI.resources <- read.csv(resource.file)

# load TF classification based on kexwords, GO terms and effector domain information from Soto
download.file("https://zenodo.org/record/7773985/files/CollecTRI_TFclassification.csv?download=1", file.path("data", "CollecTRI_TFclassification.csv"))
tf.class <- read.csv("data/CollecTRI_TFclassification.csv", sep = ";") %>%
  dplyr::rename("TF" = TFC2_Associated.Gene.Name,
                "strict" = STRICT_agreement..GO.UniProt.StructureFunction.)

## Create homogenized table ---------------------------
assign_sign <- function(aggregated.resources,
                        tf.classification = NULL, #needs to be provided if use.keywords is set to TRUE
                        use.PMIDs = T,
                        use.classification = T,
                        use.TFregulon = T,
                        unknown.to.pos = T,
                        save.decision = F,
                        regulon_end = T){

  # select sign based on more number of PMIDs
  network.signed <- aggregated.resources %>%
    mutate(sign = case_when(
      activation > repression ~ 1,
      repression > activation ~ -1,
      activation == repression ~ 0
    ))

  network.signed$TF <- unlist(lapply(str_split(network.signed$TF.TG, ":"), `[[`, 1))
  network.signed <- network.signed %>%
    mutate(sign.decision = case_when(
    sign != 0 ~ "PMID"
  ))

  if(!use.PMIDs){
    network.signed$sign <- 0
    network.signed$sign.decision <- NA
  }

  if(regulon_end){ #either run TF classification first and then regulon or other way around
    # for unknown interactions use keywords to select sign
    if (use.classification){
      network.signed <- left_join(network.signed, tf.classification %>%
                                    dplyr::select(TF, strict), by = "TF") %>%
        rename("type" = strict)

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "Act" ~ 1,
          type == "Repr" ~ -1,
          is.na(type) | type == "" ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(sign.decision = case_when(
          sign.decision == "PMID" ~ "PMID",
          sign != 0 & is.na(sign.decision) ~ "keywords"
        ))
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
        mutate(sign.decision = case_when(
          sign.decision == "PMID" ~ "PMID",
          sign.decision == "keywords" ~ "TF role",
          sign != 0 & is.na(sign.decision) ~ "regulon"
        ))

    }
  } else {
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
        mutate(sign.decision = case_when(
          sign.decision == "PMID" ~ "PMID",
          sign != 0 & is.na(sign.decision) ~ "regulon"
        ))

    }

    # for unknown interactions use keywords to select sign
    if (use.classification){
      network.signed <- left_join(network.signed, tf.classification %>%
                                    dplyr::select(TF, strict), by = "TF") %>%
        rename("type" = strict)

      network.signed[network.signed$sign == 0,] <- network.signed[network.signed$sign == 0,] %>%
        mutate(sign = case_when(
          type == "Act" ~ 1,
          type == "Repr" ~ -1,
          is.na(type) | type == "" ~ 0
        ))

      network.signed <- network.signed %>%
        mutate(sign.decision = case_when(
          sign.decision == "PMID" ~ "PMID",
          sign.decision == "regulon" ~ "regulon",
          sign != 0 & is.na(sign.decision) ~ "keywords"
        ))
    }
  }

  if (unknown.to.pos){
    network.signed$sign[network.signed$sign == 0] <- 1
    network.signed$sign.decision[is.na(network.signed$sign.decision)] <- "default activation"
  }

  network.signed <- network.signed %>%
    filter(sign != 0)

  if (save.decision){
    network.signed <- network.signed %>%
      dplyr::mutate(source = map_chr(str_split(TF.TG, ":"), 1)) %>%
      dplyr::mutate(target = map_chr(str_split(TF.TG, ":"), 2)) %>%
      dplyr::rename("weight" = sign) %>%
      dplyr::select(c(source, target, weight, TF.category, resources, PMID, sign.decision))
  }

  if (!save.decision){
    network.signed <- network.signed %>%
      dplyr::mutate(source = map_chr(str_split(TF.TG, ":"), 1)) %>%
      dplyr::mutate(target = map_chr(str_split(TF.TG, ":"), 2)) %>%
      dplyr::rename("weight" = sign) %>%
      dplyr::select(c(source, target, weight, TF.category, resources, PMID)) #sign.decision
  }


  return(network.signed)
}


## Construct  networks ---------------------------
signed.collecTRI <- assign_sign(aggregated.resources = collecTRI.resources,
                                tf.classification = tf.class,
                                use.PMIDs = T,
                                use.classification = F,
                                save.decision = T)

signed.collecTRI_TFrole <- assign_sign(aggregated.resources = collecTRI.resources,
                                tf.classification = tf.class,
                                use.PMIDs = F,
                                use.classification = T,
                                save.decision = T)
signed.collecTRI_PMID_TFrole <- assign_sign(aggregated.resources = collecTRI.resources,
                                            tf.classification = tf.class,
                                            use.PMIDs = T,
                                            use.classification = T,
                                            save.decision = T)
signed.collecTRI_PMID_regulon_TFrole <- assign_sign(aggregated.resources = collecTRI.resources,
                                                    tf.classification = tf.class,
                                                    use.PMIDs = T,
                                                    use.classification = T,
                                                    save.decision = T,
                                                    regulon_end = F)

## save networks ---------------------------
dir.create(file.path(output.folder, "CollecTRI"), showWarnings = FALSE)
dir.create(file.path(output.folder, "signed_networks"), showWarnings = FALSE)
readr::write_csv(signed.collecTRI, file.path(output.folder, "CollecTRI", "CollecTRI.csv"))
readr::write_csv(signed.collecTRI_TFrole, file.path(output.folder, "signed_networks", "CollecTRI_signed_TFrole.csv"))
readr::write_csv(signed.collecTRI_PMID_TFrole, file.path(output.folder, "signed_networks", "CollecTRI_signed_PMID_TFrole.csv"))
readr::write_csv(signed.collecTRI_PMID_regulon_TFrole, file.path(output.folder, "signed_networks", "CollecTRI_signed_PMID_regulon_TFrole.csv"))
