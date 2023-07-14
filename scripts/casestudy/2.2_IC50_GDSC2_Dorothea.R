library(tidyverse)
library(car)
library(parallel)

## Load files
Dorothea <- read.csv("data/Networks/dorothea_ABC.csv")
IC50_ulm_Dorothea <- read.csv("IC50_ulm_Dorothea.csv")

GDSC2_data <- read.delim("data/GDSC/GDSC2_fitted_dose_response.csv")
model_list <- read.csv("data/GDSC/model_list_20230608.csv")

# Select metadata for cell lines in GDSC2 cohort
GDSC2_cell_lines <- subset(model_list, model_list$model_id %in% GDSC2_data$SANGER_MODEL_ID) %>%
  select(model_id, tissue, msi_status,mutational_burden) %>%
  filter(msi_status != "" & mutational_burden != "")

GDSC2_data_flt <- GDSC2_data %>%
  select(SANGER_MODEL_ID, DRUG_ID, LN_IC50) %>%
  inner_join(GDSC2_cell_lines, by = c("SANGER_MODEL_ID" = "model_id")) %>%
  distinct()

ulm_Dorothea <- IC50_ulm_Dorothea %>%
  filter(Pval < 0.01)

## Quantify the number of TFs per cell line
ulm_Dorothea %>%
  group_by(Sample) %>%
  summarise(n())

ulm_Dorothea <- ulm_Dorothea %>%
  select(-Pval)
#  spread(TF, Activity)

names(ulm_Dorothea) <- gsub("-","_", names(ulm_Dorothea))

model_inp <- left_join(GDSC2_data_flt, ulm_Dorothea,
                       by = c("SANGER_MODEL_ID"="Sample"), relationship = "many-to-many")

model_inp_t <- model_inp %>%
  unite(split_ID, DRUG_ID, TF, remove = F) %>%
  split(~split_ID)

# Remove those TFs with an activity in less than 20% of the cell lines
model_inp_flt <- model_inp_t[sapply(model_inp_t, nrow)> 0.2 * nrow(GDSC2_cell_lines)]

#model_inp[is.na(model_inp)] <- 0

model_inp_flt <- lapply(model_inp_flt, function(x) x %>%
                          #    select(-where(function(y) all(y == 0))) %>%
                          select(-SANGER_MODEL_ID, - DRUG_ID, -split_ID,-TF))

save(model_inp_flt, file = "IC50_input_GDSC2_Dorothea.RData")

# Fit the regression model for each group and store the models
IC50_model <- function(df){
  tryCatch((model = lm(LN_IC50 ~ ., data = df)), error = function(e) NULL)
}

lm(LN_IC50 ~ ., data = model_inp_flt[[1]])
model <- IC50_model(model_inp_flt[[1]])
summary(model)

#Run models and Anova in parallel

cl <- makeCluster(detectCores())

# Load necessary packages in parallel workers
clusterEvalQ(cl, {
  library(car)
})

# Export necessary objects to parallel workers
clusterExport(cl, "Anova")

models_GDSC2_Dorothea <- parLapply(cl, model_inp_flt, function(x) lm(LN_IC50 ~ ., data = x))

anova_GDSC2_Dorothea <- parLapply(cl, models_GDSC2_Dorothea, function(x) Anova(x, type = "2"))

stopCluster(cl)

save(models_GDSC2_Dorothea, file = "IC50_models_GDSC2_Dorothea.RData")
save(anova_GDSC2_Dorothea, file = "IC50_anova_GDSC2_Dorothea.RData")

anova_adj <- lapply(anova_GDSC2_Dorothea, function(x) x %>% mutate(Adjusted_P_Value = p.adjust(.$"Pr(>F)", method = "BH")) %>% rownames_to_column("Factor"))

anova_adj_005 <- lapply(anova_adj, function(x) x %>% filter(Factor == "Activity" & Adjusted_P_Value < 0.05) )
anova_adj_005 <- anova_adj_005[sapply(anova_adj_005, nrow)>0]
save(anova_adj_005, file = "IC50_anova_005_GDSC2_Dorothea.RData")

# Merge model coefficients and Anova p-vals
merge_dataframes <- function(list1, list2) {
  merged_list <- map(names(list1), function(name) {
    inner_join(list1[[name]], list2[[name]], by = "Factor", )
  })

  return(merged_list)
}
models_coeff <- lapply(models_GDSC2_Dorothea, function(x) (x$effects) %>% stack(.) %>%
                         rename(Effect = 'values', Factor = 'ind') %>%
                         mutate(Factor = as.character(Factor)) %>%
                         filter(Factor != ""))

models_coeff <- Map(cbind, models_coeff, Combo = names(models_coeff))

res_all <- merge_dataframes(anova_adj_005, models_coeff)
res_all <- res_all %>%
  purrr::reduce(rbind)
res_all_GDSC2_Dorothea <- res_all
save(res_all_GDSC2_Dorothea, file ="IC50_GDSC2_Dorothea_resAll.RData")
res_all_flt <- res_all %>%
  separate(Combo, c("DrugA","TF"), remove = F, sep = "_")

length(unique(res_all_flt$TF))
length(unique(res_all_flt$DrugA))
