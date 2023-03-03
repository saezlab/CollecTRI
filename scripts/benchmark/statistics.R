# Copyright (c) Sophia Müller-Dott [2023]
# sophia.mueller-dott@uni-heidelberg.de

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

auprc_mat <- bench_agnositc_res %>%
  filter(metric == "mcauprc") %>%
  filter(method == "consensus_estimate") %>%
  add_column(counter = rep(c(1:1000), times = length(unique(bench_agnositc_res$net)))) %>%
  select(score, net, counter) %>%
  pivot_wider(names_from = net, values_from = score) %>%
  column_to_rownames("counter")

### Perform t-test
perform.multi.ttest <- function(mat){
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

write_csv(auroc.ttest, "output/benchmark/bench_comp_AUROC.csv")

auroc.ttest %>%
  filter(comp1 == "CollecTRI") %>%
  pull(t.value) %>%
  mean()

auroc.ttest %>%
  filter(comp1 == "shuffled CollecTRI") %>%
  filter(t.value > 0) %>%
  pull(t.value) %>%
  mean()

# AUPRC
auprc.ttest <- perform.multi.ttest(auprc_mat) %>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2)) %>%
  select(-p.value)

write_csv(auprc.ttest, "output/benchmark/bench_comp_AUPRC.csv")

auprc.ttest %>%
  filter(comp1 == "CollecTRI") %>%
  pull(t.value) %>%
  mean()

auprc.ttest %>%
  filter(comp1 == "shuffled CollecTRI") %>%
  filter(t.value > 0) %>%
  pull(t.value) %>%
  mean()



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
  mutate(p.adj = round(p.adj, digits = 2)) %>%
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
  pull(t.value) %>%
  mean

full.source.auprc.ttest <- perform.multi.ttest(full_source_auprc_mat)%>%
  mutate(t.value = round(t.value, digits = 1)) %>%
  mutate(p.adj = round(p.adj, digits = 2)) %>%
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
  pull(t.value) %>%
  mean


## Weights ---------------------------
benchmark_weights <- read_csv("output/040722/benchmark/weights_res.csv")



## Weights ---------------------------
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
      auprc_mat[c("collecTRI_signed")])
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



