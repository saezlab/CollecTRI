library(tidyverse)
library(car)
library(parallel)
library(janitor)

## Load files
IC50_ulm_CollecTRI <- read.csv("IC50_ulm_CollecTRI.csv") %>%
  filter(Pval < 0.01)

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

gene_counts <- read.csv("data/GDSC/rnaseq_fpkm_20220624_flt.csv")
gene_counts <- gene_counts[complete.cases(gene_counts),]
gene_counts[,2:ncol(gene_counts)] <- lapply(gene_counts[,2:ncol(gene_counts)], as.numeric)

GDSC2_expr <- select(gene_counts, X, any_of(GDSC2_cell_lines$model_id)) %>%
  remove_rownames() %>%
  distinct(X, .keep_all = T) %>%
  filter(X %in% unique(IC50_ulm_CollecTRI$TF)) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate( across( where(is.character), ~ as.numeric(.x) ) ) %>%
  rownames_to_column("SIDM")

GDSC2_expr[GDSC2_expr < 10] <- 0
i <- colSums(GDSC2_expr == 0, na.rm=TRUE) < 0.2 * ncol(GDSC2_expr)
GDSC2_expr <- GDSC2_expr[, i, drop=FALSE]

names(GDSC2_expr) <- gsub("-","_", names(GDSC2_expr))

GDSC2_expr <- gather(GDSC2_expr, TF, Expression,AIP:YY1)
GDSC2_expr_log <- GDSC2_expr %>% mutate(Expression = log2(Expression +1))

model_inp <- inner_join(GDSC2_data_flt, GDSC2_expr_log,
                       by = c("SANGER_MODEL_ID"="SIDM"), relationship = "many-to-many") %>%
  filter(!is.na(TF))

model_inp <- model_inp %>%
  unite(split_ID, DRUG_ID, TF, remove = F) %>%
  split(~split_ID)

# Remove those TFs with expression in less than 20% of the cell lines
model_inp <- model_inp[sapply(model_inp, nrow) > 0.2 * nrow(GDSC2_cell_lines)]

#model_inp[is.na(model_inp)] <- 0

model_inp_t <- map(model_inp, function(x) {
  x[, -c("SANGER_MODEL_ID", "DRUG_ID", "split_ID", "TF")]
})

model_inp <- lapply(model_inp, function(x)  subset(x, select = -c(SANGER_MODEL_ID, DRUG_ID, split_ID, TF)))

save(model_inp, file = "IC50_input_GDSC2_expr_tissue.RData")

lm(LN_IC50 ~ ., data = model_inp[[1]])

load("IC50_input_GDSC2_expr_tissue.RData")

#Run models and Anova in parallel

cl <- makeCluster(detectCores())

# Load necessary packages in parallel workers
clusterEvalQ(cl, {
  library(car)
})

# Export necessary objects to parallel workers
clusterExport(cl, "Anova")

load("IC50_models_GDSC2_Expression.RData")
models_GDSC2_expr <- parLapply(cl, model_inp, function(x) lm(LN_IC50 ~ ., data = x))
anova_GDSC2_expr <- parLapply(cl, models_GDSC2_expr, function(x) Anova(x, type = "2"))

stopCluster(cl)

save(models_GDSC2_expr, file = "IC50_models_GDSC2_Expression.RData")
save(anova_GDSC2_expr, file = "IC50_anova_GDSC2_Expression.RData")

anova_adj <- lapply(anova_GDSC2_expr, function(x) x %>% mutate(Adjusted_P_Value = p.adjust(.$"Pr(>F)", method = "BH")) %>% rownames_to_column("Factor"))

anova_adj_005 <- lapply(anova_adj, function(x) x %>% filter(Factor == "Expression" & Adjusted_P_Value < 0.05) )
anova_adj_005 <- anova_adj_005[sapply(anova_adj_005, nrow)>0]
save(anova_adj_005, file = "IC50_anova_005_GDSC2_expr.RData")


# Merge model coefficients and Anova p-vals
merge_dataframes <- function(list1, list2) {
  merged_list <- map(names(list1), function(name) {
    inner_join(list1[[name]], list2[[name]], by = "Factor", )
  })

  return(merged_list)
}

models_coeff <- lapply(models_GDSC2_expr, function(x) (x$effects) %>% stack(.) %>%
                         rename(Effect = 'values', Factor = 'ind') %>%
                         mutate(Factor = as.character(Factor)) %>%
                         filter(Factor != ""))

models_coeff <- Map(cbind, models_coeff, Combo = names(models_coeff))

res_all <- merge_dataframes(anova_adj_005, models_coeff)
res_all <- res_all %>%
  purrr::reduce(rbind)

res_all_GDSC2_Expression <- res_all
save(res_all_GDSC2_Expression, file ="IC50_GDSC2_Expression_resAll.RData")

res_all_flt <- res_all %>%
  separate(Combo, c("DrugA","TF"), remove = F, sep = "_")

length(unique(res_all_flt$TF))
length(unique(res_all_flt$DrugA))
