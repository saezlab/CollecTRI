# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

library(tidyverse)

## define version and path for collected TRI pairs ---------------------------
file.version <- "040722"
raw.file <- "data/CollecTRI_alias_pairs_final_from_Miguel_040722.tsv"
#Code is quite specific to the structure of this file.
#If the column names and overall structure change in a never version the code will need to be checked again.

## Create homogenized table ---------------------------
NTNU_compiled_raw <- read.table(raw.file, sep = "\t", header = TRUE)

get_homogenized_table <- function(collecTRI.raw){
  # Rename columns to have a consistent naming
  colnames(collecTRI.raw) <- colnames(collecTRI.raw) %>%
    str_replace(pattern = "\\.\\.", replacement = "_") %>%
    str_replace(pattern = "X.", replacement = "") %>%
    str_replace(pattern = "Sign", replacement = "Regulation") %>%
    str_replace(pattern = "Activation.Repression", replacement = "Regulation") %>%
    str_replace(pattern = "Mode.of.action", replacement = "Regulation")

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

      PMID_split <- unlist(str_split(PMID, ",")) %>% str_remove(" ")
      Regulation_split <- unlist(str_split(Regulation, ",")) %>% str_remove(" ")

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
    TF.TG.df$Regulation[TF.TG.df$Regulation %in% c("UP", "Activation", "positive", "activation", "+")] <- "activation"
    TF.TG.df$Regulation[TF.TG.df$Regulation %in% c("DOWN", "Repression", "negative", "-")] <- "repression"
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
    # and respression, respectively. Only unknown or repression + unknown will be
    # unknown type of regulation

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
               PMID = paste(rownames(final_TF.TG), collapse = ","))
    })
  homogenized_table <- homogenized_table %>% add_column(TF.TG = rownames(collecTRI.raw),
                                                        .before = "activation") %>%
    add_column(TF.category = collecTRI.raw$TF.Category)

  # final clean up to be sure to only have connections with associated PMIDs
  homogenized_table <- homogenized_table %>% filter(PMID != "")

  return(homogenized_table)
}

homogenized_table <- get_homogenized_table(NTNU_compiled_raw)


## Remove TFs with uncertain role ---------------------------
all.keywords <- readxl::read_excel("data/CollecTRI_TF-role_new_draft_211122.xlsx") %>%
  rename("TF" = "TFC2_Associated.Gene.Name",
         "strict" = "STRICT_agreement (GO/UniProt-StructureFunction)",
         "relaxed_KB" = "RELAXED - KB (STRICT_add aggregated GO/UniProt)",
         "relaxed_SF" =  "RELAXED - StructureFunction (STRICT_add KRAB/Soto)",
         "relaxed" = "RELAXED_KB OR StructureFunction")

rm.TFs <- all.keywords %>%
  filter(TF_category %in% c("uncertain TF role", "likely coTF")) %>%
  pull(TF)

homogenized_table <- homogenized_table %>%
  mutate(TF = map_chr(str_split(TF.TG, ":"), 1)) %>%
  filter(!TF %in% rm.TFs) %>%
  select(-TF)

## Save homogenized table ---------------------------
dir.create(file.path("output", file.version), showWarnings = FALSE)
dir.create(file.path("output", file.version, "01_homogenized_table"), showWarnings = FALSE)
write_csv(homogenized_table, file.path("output", file.version, "01_homogenized_table", "PMID_CollecTRI.csv"))
