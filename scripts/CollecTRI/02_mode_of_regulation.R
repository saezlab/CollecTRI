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
tf.class <- read.csv("data/CollecTRI_TFclassification.csv", sep = ";") %>%
  dplyr::rename("TF" = TFC2_Associated.Gene.Name,
                "strict" = STRICT_agreement..GO.UniProt.StructureFunction.)

## Create homogenized table ---------------------------
assign_sign <- function(aggregated.resources,
                        tf.classification = NULL, #needs to be provided if use.keywords is set to TRUE
                        use.classification = T,
                        use.TFregulon = T,
                        unknown.to.pos = T){

  # select sign based on more number of PMIDs
  network.signed <- aggregated.resources %>%
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
        mutate(decision = case_when(
          decision == "PMID" ~ "PMID",
          sign != 0 & is.na(decision) ~ "keywords"
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
      mutate(decision = case_when(
        decision == "PMID" ~ "PMID",
        decision == "keywords" ~ "keywords",
        sign != 0 & is.na(decision) ~ "regulon"
      ))

  }

  if (unknown.to.pos){
    network.signed$sign[network.signed$sign == 0] <- 1
    network.signed$decision[is.na(network.signed$decision)] <- "unknown"
  }

  network.signed <- network.signed %>%
    filter(sign != 0)
  network.signed <- network.signed %>%
    dplyr::mutate(source = map_chr(str_split(TF.TG, ":"), 1)) %>%
    dplyr::mutate(target = map_chr(str_split(TF.TG, ":"), 2)) %>%
    dplyr::rename("weight" = sign) %>%
    dplyr::select(c(source, target, weight, PMID, TF.category)) #decision

  return(network.signed)
}


## Construct  networks ---------------------------
signed.collecTRI <- assign_sign(aggregated.resources = collecTRI.resources,
                              tf.classification = tf.class)

## save networks ---------------------------
dir.create(file.path(output.folder, "CollecTRI"), showWarnings = FALSE)
readr::write_csv(signed.collecTRI, file.path(output.folder, "CollecTRI", "CollecTRI.csv"))
