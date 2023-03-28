#' In this script we aggregate the compiled information of the resources
#' ExTRI, HTRI, TRRUST, TFactS, GOA, IntAct, SIGNOR, CytReg, GEREDB,
#' DoRotheEA_A, Pavlidis2021 and manual curations from the NTNU.
#'
#' 1) We gather the unique PubMed references for each TF-gene interaction and
#' removed those lacking any PubMed reference.
#'
#' 2) We only include TF-gene links from TFs classified as DNA-binding
#' transcription factors (dbTFs), co-regulatory transcription factors (coTFs)
#' or general initiation transcription factors (GTFs) based on criteria from
#' TFclass, Lambert et al., Lovering et al., or gene ontology annotations.

library(tidyverse)

## Load and prepare data ---------------------------
output.folder <- "output"
raw.file <- "data/CollecTRI.tsv"

collecTRI.raw <- read.table(raw.file,
                            sep = "\t",
                            header = TRUE) # load raw resource list

# Unify column names across resources
colnames(collecTRI.raw) <- colnames(collecTRI.raw) %>%
  str_replace(pattern = "\\.\\.", replacement = "_") %>%
  str_replace(pattern = "X.", replacement = "") %>%
  str_replace(pattern = "Sign", replacement = "Regulation") %>%
  str_replace(pattern = "Activation.Repression", replacement = "Regulation") %>%
  str_replace(pattern = "Mode.of.action", replacement = "Regulation")  %>%
  str_replace(pattern = "Effect", replacement = "Regulation")

## Aggregate resources ---------------------------
aggregate_resources <- function(collecTRI.raw){
  collecTRI.raw <- collecTRI.raw %>% column_to_rownames("TF.TG")

  # All resources with PMID information
  resources <- colnames(collecTRI.raw)[str_detect(colnames(collecTRI.raw), "PMID")] %>%
    str_remove("_PMID")

  # Construction of homogenized table
  homogenized_table <- map_dfr(1:nrow(collecTRI.raw), function(i){
    TF.TG.raw <- collecTRI.raw[i,]
    TF.TG.reg <- TF.TG.raw[, str_detect(colnames(TF.TG.raw), "PMID") |
                             str_detect(colnames(TF.TG.raw), "Regulation")]

    TF.TG.df <- map_dfr(resources, function(resource){
      # Add regulation column to resources without (add unknown regulation)
      if (sum(str_detect(colnames(TF.TG.reg), resource)) == 1){
        name <- paste0(resource, "_Regulation")
        TF.TG.reg[, name] <- "Unknown"
      }

      # Convert PMID with corresponding PMIDs into data.frame for each resource
      TF.TG.res <- TF.TG.reg[,str_detect(colnames(TF.TG.reg), resource)]
      PMID <- TF.TG.res[,str_detect(colnames(TF.TG.res), "PMID")]
      Regulation <- TF.TG.res[,str_detect(colnames(TF.TG.res), "Regulation")]

      PMID_split <- unlist(str_split(PMID, "\\|")) %>% str_remove(" ")
      Regulation_split <- unlist(str_split(Regulation, "\\|")) %>% str_remove(" ")

      map_dfr(1:length(PMID_split), function(i){
        data.frame(PMID = unlist(str_split(PMID_split[i], ";")),
                   Regulation = Regulation_split[i],
                   Resource = resource)
      })

    })

    # Remove rows without PMID information
    TF.TG.df <- TF.TG.df %>%
      filter(PMID != "")

    # Consistent naming for activation, repression, unknown
    TF.TG.df$Regulation[TF.TG.df$Regulation %in% c("UP", "Activation", "positive", "activation", "+", "Stimulate")] <- "activation"
    TF.TG.df$Regulation[TF.TG.df$Regulation %in% c("DOWN", "Repression", "negative", "-", "Inhibit")] <- "repression"
    TF.TG.df$Regulation[TF.TG.df$Regulation %in% c("Unknown", "unknown", "", "not_applicable")] <- "unknown"
    TF.TG.df$Regulation[is.na(TF.TG.df$Regulation)] <- "unknown"

    # Remove duplicated rows (same PMID and same type of regulation)
    filtered_TF.TG <- TF.TG.df %>% filter(!duplicated(TF.TG.df %>% select(PMID, Regulation)))

    # Table data.frame to find PMID with multiple information about regulation
    table_TF.TG <- table(filtered_TF.TG$PMID, filtered_TF.TG$Regulation) %>%
      as.data.frame.matrix()

    if (!"repression" %in% colnames(table_TF.TG)){
      table_TF.TG <- table_TF.TG %>% add_column("repression" = 0)
    }

    if (!"activation" %in% colnames(table_TF.TG)){
      table_TF.TG <- table_TF.TG %>% add_column("activation" = 0)
    }

    # change PMID with activation + unknown or repression + unknown to activation
    # or repression, respectively

    final_TF.TG <- table_TF.TG %>%
      mutate(Sign = case_when(
        activation >= 1 & repression == 0 ~ "activation",
        activation == 0 & repression >= 1 ~ "repression",
        activation >= 1 & repression >= 1 ~ "activation+repression"
      ))

    final_TF.TG$Sign[is.na(final_TF.TG$Sign)] <- "unknown"

    # create final data frame with number of unique PMID for activation, repression
    # and unknown + all PMIDs
    data.frame(activation = sum(str_detect(final_TF.TG$Sign, "activation")),
               repression = sum(str_detect(final_TF.TG$Sign, "repression")),
               unknown = sum(final_TF.TG$Sign == "unknown"),
               PMID = paste(rownames(final_TF.TG), collapse = ";"),
               resources = paste(unique(TF.TG.df$Resource), collapse = ";"))
  })
  homogenized_table <- homogenized_table %>% add_column(TF.TG = rownames(collecTRI.raw),
                                                        .before = "activation") %>%
    add_column(TF.category = collecTRI.raw$TF.Category)

  # final clean up to be sure to only have connections with associated PMIDs
  homogenized_table <- homogenized_table %>% filter(PMID != "")

  return(homogenized_table)
}

resource.table.raw <- aggregate_resources(collecTRI.raw)

## Filter proteins not classified as TFs ---------------------------
# Add information about TF category
resource.table <- resource.table.raw %>%
  left_join(collecTRI.raw %>%
              dplyr::select(TF.TG, Lambert, Lovering, TFClass,
                     GO.0003700, GO.0140223, GO.0003712),
            by = "TF.TG") %>%
  dplyr::mutate(TF.category = case_when(
    Lambert != "" | Lovering != "" | TFClass != "" | GO.0003700 != "" ~ "dbTF",
    GO.0003712 != "" ~ "coTF",
    GO.0140223 != "" ~ "GTF"
  )) %>%
  dplyr::select(!c(Lambert, Lovering, TFClass, GO.0003700, GO.0140223, GO.0003712))

# Filter out interactions from non-TFs
resource.table <- resource.table[!is.na(resource.table$TF.category),]

## Save data ---------------------------
dir.create(file.path(output.folder, "aggregated_resources"), showWarnings = FALSE)

write_csv(resource.table.raw, file.path(output.folder, "aggregated_resources", "CollecTRI_resources_raw.csv"))
write_csv(resource.table, file.path(output.folder, "aggregated_resources", "CollecTRI_resources.csv"))
