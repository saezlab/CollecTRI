NTNU_compiled_raw <- read.table("data/NTNU_BSC_compiled.tsv", sep = "\t", header = TRUE)

colnames(NTNU_compiled_raw) <- c(
  "TF.TG", "TF", "TG", "ExTRI_Confidence", "ExTRI_PMID", "ExTRI_present",
  "HTRI_present", "HTRI_Technique", "HTRI_PMID", "HTRI_Confidence",
  "TRRUST_present", "TRRUST_Regulation", "TRRUST_PMID",
  "TFactS_present", "TFactS_Regulation", "TFactS_Species", "TFactS_Source", "TFactS_PMID", "TFactS_Confidence",
  "GOA_present", "GOA_Sign",
  "IntAct_present", "IntAct_PMID", "IntAct_Method.ID",
  "SIGNOR_present", "SIGNOR_Effect", "SIGNOR_Regulation", "SIGNOR_PMID",
  "CytReg_present", "CytReg_Assay.type", "CytReg_species", "CytReg_Regulation", "CytReg_PMID", "CytReg_Year.of.publication",
  "GEREDB_present", "GEREDB_Regulation", "GEREDB_PMID",
  "NTNU.Curated_present", "NTNU.Curated_Regulation", "NTNU.Curated_PMID", "TFClass", "Auto.regulation"
)

NTNU_raw <- NTNU_compiled_raw[str_detect(colnames(NTNU_compiled_raw), "PMID") | str_detect(colnames(NTNU_compiled_raw), "Regulation")]
rownames(NTNU_raw) <- NTNU_compiled_raw$TF.TG

head(NTNU_raw)
ressources <- unique(sapply(str_split(colnames(NTNU_raw), "_"), "[[", 1))

summarized_tables <- map(ressources, function(x) {
  counts_ressource <- NTNU_raw[str_detect(colnames(NTNU_raw), x)]
  colnames(counts_ressource) <- str_remove(colnames(counts_ressource), paste0(x, "_"))

  if (!any(colnames(counts_ressource) == "Regulation")) {
    counts_ressource <- cbind(counts_ressource, Regulation = "Unknown")
  }

  counts_ressource <- counts_ressource %>% filter(!PMID == "")
  PMID_ressource <- counts_ressource[colnames(counts_ressource) == "PMID"]

  Regulation_ressource <- counts_ressource[colnames(counts_ressource) == "Regulation"]

  summarized_matrix <- matrix(ncol = 3, nrow = nrow(counts_ressource))
  colnames(summarized_matrix) <- c("Positive", "Negative", "Unknown")
  rownames(summarized_matrix) <- rownames(counts_ressource)

  for (i in 1:nrow(PMID_ressource)) {
    counts <- unlist(str_split(PMID_ressource[i, ], "\\|"))
    names(counts) <- unlist(str_split(Regulation_ressource[i, ], "\\|"))
    counts[counts == ""] <- paste0("Unknown_", 1:sum(counts == ""))

    summarized_matrix[i, 1] <- length(unique(unlist(str_split(counts[names(counts) %in% c("UP", "Activation", "positive")], ";"))))
    summarized_matrix[i, 2] <- length(unique(unlist(str_split(counts[names(counts) %in% c("DOWN", "Repression", "negative")], ";"))))
    summarized_matrix[i, 3] <- length(unique(unlist(str_split(counts[names(counts) %in% c("Unknown", "unknown", "")], ";"))))
  }

  colnames(summarized_matrix) <- paste0(x, "_", colnames(summarized_matrix))

  summarized_df <- data.frame(summarized_matrix) %>% rownames_to_column("TF.TG")

  return(summarized_df)
})

homogenized_table <- summarized_tables %>% reduce(full_join, by = "TF.TG") %>% column_to_rownames("TF.TG")
homogenized_table[is.na(homogenized_table)] <- 0
homogenized_table <- homogenized_table %>%
                      add_column(totNeg = rowSums(homogenized_table[str_detect(colnames(homogenized_table), "Positive")])) %>%
                      add_column(totPositive = rowSums(homogenized_table[str_detect(colnames(homogenized_table), "Negative")])) %>%
                      add_column(totUnknown = rowSums(homogenized_table[str_detect(colnames(homogenized_table), "Unknown")])) %>%
                      rownames_to_column("TF.TG")
write_csv(homogenized_table, "data/homogenized_ressource_ExTRI.csv")
