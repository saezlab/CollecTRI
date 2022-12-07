# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we filter the networks
#' to test the effects of removing dimers (NFkB and AP1)

library(tidyverse)

## load network ---------------------------
file.version <- "040722"

collecTRI <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))


## Removal of dimers ---------------------------
dimers_to_remove <- c("NFKB", "AP1")

collecTRI_filtered <- collecTRI %>%
  filter(!source %in% dimers_to_remove)

print(paste0("Remove of dimers: ", nrow(collecTRI) - nrow(collecTRI_filtered), " edges removed"))

dir.create(file.path("output", file.version, "05_filtered_networks"), showWarnings = FALSE)
write_csv(collecTRI_filtered, file.path("output", file.version, "05_filtered_networks", "collecTRI_dimers_removed.csv"))
