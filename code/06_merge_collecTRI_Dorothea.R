# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we merge the top
#' performing methods
#' 1) CollecTRI
#' 2) Dorothea
#' We try different merging strategies and benchmark them


library(tidyverse)
library(ggbreak)

## Load data ---------------------------
set.seed(123)
file.version <- "040722"

collecTRI <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))
dorothea <- read_csv("data/raw/dorothea_ABC.csv")[2:5]
mat <- read_csv('data/knockTF_expr.csv')
obs <- read_csv('data/knockTF_meta.csv')
msk <- obs$logFC < -1
obs <- obs[msk,]

TFs_bench <- obs$TF %>% unique()

all.keywords <- readxl::read_excel("data/CollecTRI_TF-role_new_draft_211122.xlsx") %>%
  rename("source" = "TFC2_Associated.Gene.Name",
         "strict" = "STRICT_agreement (GO/UniProt-StructureFunction)",
         "relaxed_KB" = "RELAXED - KB (STRICT_add aggregated GO/UniProt)",
         "relaxed_SF" =  "RELAXED - StructureFunction (STRICT_add KRAB/Soto)",
         "relaxed" = "RELAXED_KB OR StructureFunction") %>%
  select(source, strict)


## Merge networks ---------------------------
merging_strategies <- c("raw", "redo_unknowns", "redo_dorothea", "random")

collecTRI <- collecTRI %>%
  mutate(edge = paste(source, target, sep = ":"))

dorothea <- dorothea %>%
  mutate(edge = paste(source, target, sep = ":")) %>%
  mutate(weight = case_when(
    weight > 0 ~ 1,
    weight < 0 ~ -1
  ))

merged_networks <- map(merging_strategies, function(merge_strategy){
  merged <- full_join(collecTRI, dorothea, by = c("edge", "source", "target"))
  if (merge_strategy == "raw"){
    merged <- merged %>%
      mutate(weight = case_when(
        !is.na(weight.x) ~ weight.x,
        is.na(weight.x) ~  weight.y
      ))
  } else if (merge_strategy == "redo_unknowns"){
    merged <- merged %>%
      mutate(weight = case_when(
        !is.na(weight.x) & decision != "unknown" ~ weight.x,
        !is.na(weight.x) & is.na(weight.y) & decision == "unknown" ~ weight.x,
        is.na(weight.x) ~  weight.y,
        decision == "unknown" & !is.na(weight.y) ~ weight.y
      ))
  } else if (merge_strategy == "redo_dorothea"){
    merged <- left_join(merged, all.keywords) %>%
      mutate(weight.y = case_when(
        is.na(strict) ~ weight.y,
        strict == "Act" ~ 1,
        strict == "Repr" ~ -1
      )) %>%
      mutate(weight = case_when(
        !is.na(weight.x) ~ weight.x,
        is.na(weight.x) & is.na(strict) ~ 0,
        is.na(weight.x) & !is.na(strict) ~ weight.y
      ))
    TF_regulon <- merged %>%
      group_by(source) %>%
      summarize(ratio_act = mean(weight)) %>%
      mutate(TF_role = case_when(
        ratio_act > 0 ~ "activator",
        ratio_act < 0 ~ "repressor"
      ))

    merged <- left_join(merged, TF_regulon) %>%
      mutate(weight = case_when(
        weight != 0 ~ weight,
        weight == 0 & TF_role == "activator" ~ 1,
        weight == 0 & TF_role == "reproessor" ~ -1
      ))
    merged$weight[is.na(merged$weight)] <- 1
  } else if (merge_strategy == "random"){
    doro <- merged %>%
      filter(is.na(weight.x))

    doro_new <- map_dfr(unique(doro$source), function(tf){
      genes_tf <- merged %>% filter(source == tf) %>% pull(target)
      genes_accepted <- setdiff(unique(merged$target), genes_tf)
      doro_sub <- doro %>% filter(source == tf)

      doro_sub$target <- sample(genes_accepted, nrow(doro_sub))
      doro_sub$edge <- paste(doro_sub$source, doro_sub$target, sep = ":")

      doro_sub
    })

    doro_new <- doro_new %>% select(source, target, edge, confidence, weight.y)

    collectri <- merged %>% filter(!is.na(weight.x)) %>%
      select(source, target, weight.x, PMID, TF.category, decision, edge)

    merged <- full_join(collectri, doro_new) %>%
      mutate(weight = case_when(
        !is.na(weight.x) & decision != "unknown" ~ weight.x,
        !is.na(weight.x) & is.na(weight.y) & decision == "unknown" ~ weight.x,
        is.na(weight.x) ~  weight.y,
        decision == "unknown" & !is.na(weight.y) ~ weight.y
      ))
  }

 merged <- merged %>%
   select(source, target, weight, PMID, TF.category, decision, confidence) %>%
   rename("confidence.dorothea" = "confidence")

  dir.create(file.path("output", file.version, "06_merged_network"), showWarnings = FALSE)
  write_csv(merged, file.path("output", file.version, "06_merged_network", paste0(merge_strategy, "_doro_collecTRI.csv")))

  merged
})

names(merged_networks) <- merging_strategies
merged <- merged_networks$redo_unknowns
shared_TFs <- unique(collecTRI$source)[unique(collecTRI$source) %in% unique(dorothea$source)]

n_shared <- sum(rowSums(is.na(merged)) == 0)
n_doro <- nrow(dorothea)-n_shared
n_collecTRI <- nrow(collecTRI)-n_shared

n_shared_TF <- sum(unique(collecTRI$source) %in% unique(dorothea$source))
n_doro_TF <- length(unique(dorothea$source))-n_shared_TF
n_collecTRI_TF <- length(unique(collecTRI$source))-n_shared_TF
n_shared_TF_bench <- sum(shared_TFs %in% TFs_bench)
n_doro_TF_bench <- sum(unique(dorothea$source)[!unique(dorothea$source) %in% shared_TFs] %in% TFs_bench)
n_collecTRI_TF_bench <- sum(unique(collecTRI$source)[!unique(collecTRI$source) %in% shared_TFs] %in% TFs_bench)



plot.df <- data.frame(x = c(n_shared,
                            n_doro,
                            n_collecTRI,
                            n_shared_TF,
                            n_doro_TF,
                            n_collecTRI_TF,
                            n_shared_TF_bench,
                            n_doro_TF_bench,
                            n_collecTRI_TF_bench),
              network = rep(c("shared", "dorothea", "collecTRI"), times = 3),
              label = rep(c("edge", "TF", "TF in bench"), each = 3))

ggplot(plot.df %>% filter(!label == "edge")) +
  aes(fill=network, y=x, x=label) +
  geom_bar(stat="identity", color="black", position="stack", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  theme(text = element_text(size = 12)) + theme(legend.title = element_blank())

ggplot(plot.df %>% filter(label == "edge")) +
  aes(fill=network, y=x, x=label) +
  geom_bar(stat="identity", color="black", position="stack", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  theme(text = element_text(size = 12)) + theme(legend.title = element_blank())

