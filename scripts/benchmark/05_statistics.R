# Copyright (c) Sophia Müller-Dott [2023]
# sophia.mueller-dott@uni-heidelberg.de

# In this script we perform the statistical analysis of our results

library(tidyverse)

## Comparison of benchmark results---------------------------
### load data
bench_agnositc_res <- read_csv("output/benchmark/benchmark_res.csv")  %>%
  arrange(factor(net, levels = c("collecTRI", "ABC", "regnet",
                                 "pathComp", "chea3_archs4", "chea3_GTEx",
                                 "chea3_enrich", "rand", "chea3_remap",
                                 "chea3_encode", "chea3_lit"))) %>%
  mutate(net = recode(net,
                       chea3_archs4 = "ChEA3 ARCHS4",
                       chea3_GTEx = "ChEA3 GTEx",
                       chea3_enrich = "ChEA3 Enrichr",
                       regnet = "RegNetwork",
                       ABC = "DoRothEA ABC",
                       collecTRI = "CollecTRI",
                       chea3_remap = "ChEA3 ReMap",
                       chea3_lit = "ChEA3 Literature",
                       chea3_encode = "ChEA3 ENCODE",
                       rand = "shuffled CollecTRI",
                       pathComp = "Pathway Commmons"))

### change format into matrix for AUROC and AUPRC
auroc_mat <- bench_agnositc_res %>%
  filter(metric == "mcauroc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(bench_agnositc_res$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

auroc_mat$CollecTRI %>% median()

auprc_mat <- bench_agnositc_res %>%
  filter(metric == "mcauprc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(bench_agnositc_res$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

auprc_mat$CollecTRI %>% median()

multi.ttest <- function(mat, pVal = T, alternative = "two.sided") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- t.test(mat[, i], mat[, j], alternative = "two.sided")
      if(pVal){
        p.mat[i, j] <- p.mat[j, i] <- test$p.value
      } else {
        p.mat[i, j] <- p.mat[j, i] <- test$statistic
      }

    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

### Perform t-test
perform.multi.ttest <- function(mat){

  # extract p_values
  ttest.p <- multi.ttest(mat) %>%
    as.data.frame() %>%
    rownames_to_column("comp1") %>%
    pivot_longer(!comp1, names_to = "comp2", values_to = "p.value")

  # extract t_values
  ttest.t <- multi.ttest(mat, pVal = F) %>%
    as.data.frame() %>%
    rownames_to_column("comp1") %>%
    pivot_longer(!comp1, names_to = "comp2", values_to = "t.value")

  # Adjust t-values direction
  # Identify positions were direction needs to be adjusted
  idx <- c()
  for (i in 1:(length(unique(ttest.p$comp1))-1)){
    new_idx <- length(unique(ttest.t$comp1))*i + rep(1:i)
    idx <- append(idx, new_idx)
  }
  ttest.t$t.value[idx] <- ttest.t$t.value[idx]*-1

  ttest <- ttest.p %>%
    add_column(p.adj = p.adjust(ttest.p$p.value, method = "BH")) %>%
    left_join(ttest.t, by = c("comp1", "comp2")) %>%
    filter(!comp1 == comp2)

  # remove duplicated comparisons
  idx_rm <- c()
  for (i in 1:(length(unique(ttest.t$comp1))-1)){
    new_idx <- (length(unique(ttest.t$comp1))-1)*i + rep(1:i)
    idx_rm <- append(idx_rm, new_idx)
  }

  ttest <- ttest %>%
    slice(-idx_rm)

  return(ttest)
}

# AUROC
auroc.ttest <- perform.multi.ttest(auroc_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2)) %>%
  select(-p.value)

auroc.ttest <- auroc.ttest %>%
  mutate(comp1_new = case_when(
  t.value > 0 ~ comp1,
  t.value < 0 ~ comp2))  %>%
  mutate(comp2_new = case_when(
    t.value > 0 ~ comp2,
    t.value < 0 ~ comp1))  %>%
  mutate(t.value_new = abs(t.value)) %>%
  select(comp1_new, comp2_new, p.adj, t.value_new) %>%
  rename("comp1" = comp1_new, "comp2" = comp2_new, "t.value" = t.value_new)

auroc.ttest <- auroc.ttest %>% arrange(comp1, comp2)

auroc.ttest %>%
  filter(comp1 == "CollecTRI") %>%
  pull(t.value) %>%
  mean()

auroc.ttest %>%
  filter(comp1 == "shuffled CollecTRI") %>%
  filter(t.value > 0) %>%
  filter(p.adj < 0.05) %>%
  pull(t.value) %>%
  mean()

# AUPRC
auprc.ttest <- perform.multi.ttest(auprc_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2)) %>%
  select(-p.value)

auprc.ttest <- auprc.ttest %>%
  mutate(comp1_new = case_when(
    t.value > 0 ~ comp1,
    t.value < 0 ~ comp2))  %>%
  mutate(comp2_new = case_when(
    t.value > 0 ~ comp2,
    t.value < 0 ~ comp1))  %>%
  mutate(t.value_new = abs(t.value)) %>%
  select(comp1_new, comp2_new, p.adj, t.value_new) %>%
  rename("comp1" = comp1_new, "comp2" = comp2_new, "t.value" = t.value_new)

auprc.ttest <- auprc.ttest %>% arrange(comp1, comp2)

auprc.ttest %>%
  filter(comp1 == "CollecTRI") %>%
  pull(t.value) %>%
  mean()

auprc.ttest %>%
  filter(comp1 == "shuffled CollecTRI") %>%
  filter(t.value > 0) %>%
  filter(p.adj < 0.05) %>%
  pull(t.value) %>%
  mean()

# Merge AUROC and AUPRC results into one table
statistics_bench <- full_join(auroc.ttest, auprc.ttest, by = c("comp1", "comp2"))
colnames(statistics_bench) <- c("GRN 1", "GRN 2", "AUROC adjusted p value", "AUROC t value", "AUPRC adjusted p value", "AUPRC t value")
#as Wolfgang Huber recommended we report that the values are below detection limit (https://www.sciencedirect.com/science/article/pii/S2405471219300717?via%3Dihub)
statistics_bench$`AUROC adjusted p value`[statistics_bench$`AUROC adjusted p value` == 0] <- "<2.2x10^-16"
statistics_bench$`AUPRC adjusted p value`[statistics_bench$`AUPRC adjusted p value` == 0] <- "<2.2x10^-16"

write_csv(statistics_bench, "output/benchmark/benchmark_ttest.csv")


## Comparison of source benchmark results---------------------------
bench_source_res <- read_csv("output/benchmark/benchmark_source_res.csv")  %>%
  arrange(factor(net, levels = c("collecTRI", "ABC", "regnet"))) %>%
  mutate(net = recode(net,
                      regnet = "RegNetwork",
                      ABC = "DoRothEA ABC",
                      collecTRI = "CollecTRI"))

### change format into matrix for AUROC and AUPRC
source_auroc_mat <- bench_source_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(source, net) %>%
  summarise(median = median(score))  %>%
  pivot_wider(names_from = net, values_from = median) %>%
  column_to_rownames("source")

source_auprc_mat <- bench_source_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  group_by(source, net) %>%
  summarise(median = median(score))  %>%
  pivot_wider(names_from = net, values_from = median) %>%
  column_to_rownames("source")


### Perform t.test
source.auroc.ttest <- perform.multi.ttest(source_auroc_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2))

source.auroc.ttest

source.auprc.ttest <- perform.multi.ttest(source_auprc_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2))

source.auprc.ttest


## Repeat with full comparison
full_source_auroc_mat <- bench_source_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  filter(!source == "REST") %>%
  mutate(comb = paste(net, source, sep = ":"))  %>%
  add_column(counter = rep(c(1:1000), times = length(unique(paste(bench_source_res$net, bench_source_res$source, sep = ":")))-3)) %>%
  select(score, comb, counter) %>%
  pivot_wider(names_from = comb, values_from = score) %>%
  column_to_rownames("counter")

full_source_auprc_mat <- bench_source_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  filter(!source == "REST") %>%
  mutate(comb = paste(net, source, sep = ":"))  %>%
  add_column(counter = rep(c(1:1000), times = length(unique(paste(bench_source_res$net, bench_source_res$source, sep = ":")))-3)) %>%
  select(score, comb, counter) %>%
  pivot_wider(names_from = comb, values_from = score) %>%
  column_to_rownames("counter")

### Perform t.test
full.source.auroc.ttest <- perform.multi.ttest(full_source_auroc_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(net1 = map_chr(str_split(comp1, ":"), 1)) %>%
  mutate(TF1 = map_chr(str_split(comp1, ":"), 2)) %>%
  mutate(net2 = map_chr(str_split(comp2, ":"), 1)) %>%
  mutate(TF2 = map_chr(str_split(comp2, ":"), 2)) %>%
  filter(TF1 == TF2)

full.source.auroc.ttest %>%
  filter(net1 == "CollecTRI") %>%
  filter(t.value > 0) %>%
  filter(p.adj < 0.05) %>%
  filter(TF1 %in% c("TP53", "FLI1", "NR2F2", "SOX2")) %>%
  pull(p.adj) %>%
  range

full.source.auroc.ttest %>%
  filter(net1 == "CollecTRI") %>%
  filter(t.value > 0) %>%
  filter(p.adj < 0.05) %>%
  filter(TF1 %in% c("TP53", "FLI1", "NR2F2", "SOX2")) %>%
  pull(t.value) %>%
  mean

mean(source_auroc_mat$CollecTRI[rownames(source_auroc_mat) %in% c("TP53", "FLI1", "NR2F2", "SOX2", "REST")])

full.source.auprc.ttest <- perform.multi.ttest(full_source_auprc_mat)%>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(net1 = map_chr(str_split(comp1, ":"), 1)) %>%
  mutate(TF1 = map_chr(str_split(comp1, ":"), 2)) %>%
  mutate(net2 = map_chr(str_split(comp2, ":"), 1)) %>%
  mutate(TF2 = map_chr(str_split(comp2, ":"), 2)) %>%
  filter(TF1 == TF2)

full.source.auprc.ttest %>%
  filter(net1 == "CollecTRI") %>%
  filter(t.value > 0) %>%
  filter(p.adj < 0.05) %>%
  filter(TF1 %in% c("TP53", "FLI1", "NR2F2", "SOX2")) %>%
  pull(p.adj) %>%
  range

full.source.auprc.ttest %>%
  filter(net1 == "CollecTRI") %>%
  filter(t.value > 0) %>%
  filter(p.adj < 0.05) %>%
  filter(TF1 %in% c("TP53", "FLI1", "NR2F2", "SOX2")) %>%
  pull(t.value) %>%
  mean

mean(source_auprc_mat$CollecTRI[rownames(source_auroc_mat) %in% c("TP53", "FLI1", "NR2F2", "SOX2", "REST")])



## Weights ---------------------------
benchmark_weights <- read_csv("output/benchmark/benchmark_weights_res.csv")

auroc_mat_weights <- benchmark_weights %>%
  filter(metric == "mcauroc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(benchmark_weights$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")
auprc_mat_weights <- benchmark_weights %>%
  filter(metric == "mcauprc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(benchmark_weights$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")
auroc.ttest.p_weights <- multi.ttest(auroc_mat_weights, alternative = "less") %>%
  as.data.frame() %>%
  rownames_to_column("net1") %>%
  pivot_longer(!net1, names_to = "net2", values_to = "p.value")
auroc.ttest.t_weights <- multi.ttest(auroc_mat_weights, pVal = F, alternative = "less") %>%
  as.data.frame() %>%
  rownames_to_column("net1") %>%
  pivot_longer(!net1, names_to = "net2", values_to = "t.value")

# Adjust t-values direction
# Identify positions were direction needs to be adjusted
idx <- c()
for (i in 1:(length(unique(auroc.ttest.t_weights$net1))-1)){
  new_idx <- length(unique(auroc.ttest.t_weights$net1))*i + rep(1:i)
  idx <- append(idx, new_idx)
}
auroc.ttest.t_weights$t.value[idx] <- auroc.ttest.t_weights$t.value[idx]*-1

auroc.ttest_weights <- auroc.ttest.p_weights %>%
  add_column(p.adj = p.adjust(auroc.ttest.p_weights$p.value, method = "BH")) %>%
  left_join(auroc.ttest.t_weights, by = c("net1", "net2")) %>%
  filter(!net1 == net2)
auroc.ttest_weights %>% filter(net1 == "collecTRI")

comp_full_coverage <- cbind(auprc_mat_weights[c("FIMO10gene", "matRid10gene")],
      auprc_mat[c("CollecTRI")])
auroc.ttest.p_coverage <- multi.ttest(comp_full_coverage) %>%
  as.data.frame() %>%
  rownames_to_column("net1") %>%
  pivot_longer(!net1, names_to = "net2", values_to = "p.value")
auroc.ttest.t_coverage <- multi.ttest(comp_full_coverage, pVal = F) %>%
  as.data.frame() %>%
  rownames_to_column("net1") %>%
  pivot_longer(!net1, names_to = "net2", values_to = "t.value")

idx <- c()
for (i in 1:(length(unique(auroc.ttest.t_coverage$net1))-1)){
  new_idx <- length(unique(auroc.ttest.t_coverage$net1))*i + rep(1:i)
  idx <- append(idx, new_idx)
}
auroc.ttest.t_coverage$t.value[idx] <- auroc.ttest.t_coverage$t.value[idx]*-1

auroc.ttest_coverage <- auroc.ttest.p_coverage %>%
  add_column(p.adj = p.adjust(auroc.ttest.p_coverage$p.value, method = "BH")) %>%
  left_join(auroc.ttest.t_coverage, by = c("net1", "net2")) %>%
  filter(!net1 == net2)
auroc.ttest_coverage




## Sign ---------------------------
tf.class <- read.csv("data/CollecTRI_TFclassification.csv", sep = ";") %>%
  dplyr::rename("TF" = TFC2_Associated.Gene.Name,
                "strict" = STRICT_agreement..GO.UniProt.StructureFunction.)
collectri <- read_csv("output/CollecTRI/CollecTRI_GRN.csv")

pmid_sign <- collectri %>%
  filter(sign_decision == "PMID") %>%
  pull(weight) %>%
  table()

pmid_sign
pmid_sign["1"]/(pmid_sign["1"] + pmid_sign["-1"])

act <- tf.class %>%
  filter(strict == "Act") %>%
  pull(TF)

rep <- tf.class %>%
  filter(strict == "Repr") %>%
  pull(TF)

pos <- collectri %>% filter(source %in% act) %>% nrow()
neg <- collectri %>% filter(source %in% rep) %>% nrow()

c(pos, neg)
pos/(pos+neg)

## benchmarl sign
benchmark_sign <- read_csv("output/benchmark/benchmark_sign_res.csv")

auroc_mat_sign <- benchmark_sign %>%
  filter(metric == "mcauroc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(benchmark_sign$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")
auprc_mat_sign <- benchmark_sign %>%
  filter(metric == "mcauprc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(benchmark_sign$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

### Perform t.test
source.auroc.ttest <- perform.multi.ttest(auroc_mat_sign) %>%
  mutate(t.value = round(t.value, digits = 1))

source.auroc.ttest

source.auprc.ttest <- perform.multi.ttest(auprc_mat_sign) %>%
  mutate(t.value = round(t.value, digits = 1))

source.auprc.ttest


## Size effect ----
# The effect of the number of targets and the estimated activities was investigated
# directly in the figures_manuscript.R script. For more details please check
# out the Section Supp 1 Bias, S1.1 Size difference between TFs in benchmark and background
