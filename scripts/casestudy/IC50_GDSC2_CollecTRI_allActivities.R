library(dplyr)
library(data.table)

## Load files
CollecTRI <- read.csv("CollecTRI.csv")
IC50_ulm_CollecTRI <- read.csv("IC50_ulm_CollecTRI.csv")

GDSC2_data <- read.delim("GDSC2_fitted_dose_response.csv")
model_list <- read.csv("model_list_20230608.csv")

# Select metadata for cell lines in GDSC2 cohort
GDSC2_cell_lines <- subset(model_list, model_list$model_id %in% GDSC2_data$SANGER_MODEL_ID) %>%
  select(model_id, tissue, msi_status,mutational_burden) %>%
  filter(msi_status != "" & mutational_burden != "")

GDSC2_data_flt <- GDSC2_data %>%
  select(SANGER_MODEL_ID, DRUG_ID, LN_IC50) %>%
  inner_join(GDSC2_cell_lines, by = c("SANGER_MODEL_ID" = "model_id")) %>%
  distinct()

IC50_ulm_CollecTRI <- IC50_ulm_CollecTRI %>%
  select(-Pval) %>%
  rename(SANGER_MODEL_ID = "Sample")
#  spread(TF, Activity)

names(IC50_ulm_CollecTRI) <- gsub("-","_", names(IC50_ulm_CollecTRI)) #substitute symbols that return an error with the Anova function

GDSC2_data_flt <- as.data.table(GDSC2_data_flt)
IC50_ulm_CollecTRI <- as.data.table(IC50_ulm_CollecTRI)

# Set key columns for joining
setkey(GDSC2_data_flt, SANGER_MODEL_ID)
setkey(IC50_ulm_CollecTRI, SANGER_MODEL_ID)

# allow.cartesian=TRUE for many-to-many merging
model_inp <- GDSC2_data_flt[IC50_ulm_CollecTRI, allow.cartesian = TRUE]

model_inp[, split_ID := paste(DRUG_ID, TF, sep = "_")]

# Removing SANGER_MODEL_ID, DRUG_ID and TF columns
model_inp[, c("SANGER_MODEL_ID","DRUG_ID", "TF") := NULL]

test[, fwrite(copy(.SD)[, split_ID := split_ID], paste0("input", var1,".csv")), by = split_ID]

# Grouping the model_inp data table by split_ID
setkey(model_inp, split_ID)

# Define a function to run linear model on each group
run_linear_model <- function(dt) {
  lm(LN_IC50 ~ tissue + msi_status + mutational_burden + Activity, data = dt)
}

# Run linear model on each split_ID group
model_results <- model_inp[, .(lm_result = run_linear_model(.SD)), by = split_ID]

