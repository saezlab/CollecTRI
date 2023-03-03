# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will construct the figures for the final manuscript

library(tidyverse)
library(UpSetR)
library(ggsignif)
library(rstatix)
library(ggpubr)
library(MLeval)
library(caret)
library(RColorBrewer)

## Load data---------------------------
### networks
doro <- read_csv("data/networks/filtered_dorothea_ABC.csv")
regnet <- read_csv("data/networks/filtered_regnetwork.csv")
pathComp <- read_csv("data/networks/filtered_pathwayCommons.csv")
chea <- read_csv("data/networks/filtered_chea3.csv")
chea_arch <- chea %>%
  filter(confidence == "ARCHS4_Coexpression")
chea_encode <- chea %>%
  filter(confidence == "ENCODE_ChIP-seq")
chea_enrichr <- chea %>%
  filter(confidence == "Enrichr_Queries")
chea_GTEx <- chea %>%
  filter(confidence == "GTEx_Coexpression")
chea_lit <- chea %>%
  filter(confidence == "Literature_ChIP-seq")
chea_remap <- chea %>%
  filter(confidence == "ReMap_ChIP-seq")

collecTRI <- read_csv("output/CollecTRI/CollecTRI.csv")

#merged <- read_csv("output/230223/06_merged_network/redo_unknowns_doro_collecTRI.csv")

networks <- list(doro = doro,
                 regnet = regnet,
                 pathComp = pathComp,
                 chea_arch = chea_arch,
                 chea_encode = chea_encode,
                 chea_enrichr = chea_enrichr,
                 chea_GTEx = chea_GTEx,
                 chea_lit = chea_lit,
                 chea_remap = chea_remap,
                 collecTRI = collecTRI)

### meta data benchmark
obs <- read_csv('data/knockTF_meta.csv')
msk <- obs$logFC < -1
obs_filtered <- obs[msk,]

### TF activities
decoupler_path <- list.files(file.path("output", 'benchmark'), pattern = 'activity', full.names = T)
act <-  map(decoupler_path, read_csv)
names(act) <-  str_remove(list.files(file.path("output", 'benchmark'),  pattern = 'activity'), "_activity.csv")

## Figure 1 construction of collecTRI ---------------------------
### 1.1 Overview construction of collecTRI

### 1.2 CollecTRI composition
sign_comp_df <- table(collecTRI$weight) %>%
  as.data.frame() %>%
  mutate(regulation = case_when(
    Var1 == "-1" ~ "inhibition",
    Var1 == "1" ~ "activation"
  ))
100/sum(sign_comp_df$Freq) * sign_comp_df$Freq[1]
100/sum(sign_comp_df$Freq) * sign_comp_df$Freq[2]

p_1.2.1 <- ggplot(data=sign_comp_df, aes(x=Freq, y=regulation, fill=regulation)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.85), width = 0.7) +
  theme_minimal() +
  scale_fill_manual(values=c('#91b57c','#a36a69')) +
  ylab("") +
  xlab("Total size") +
  theme(text = element_text(size = 9),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip() +
  guides(fill=guide_legend(title=""))


mor_TFs_df <- collecTRI %>%
  group_by(source) %>%
  summarize(n_pos = sum(weight > 0),
            n_neg = sum(weight < 0)) %>%
  mutate(mor = case_when(
    n_pos > 0 & n_neg == 0 ~ "activator",
    n_pos > 0 & n_neg > 0 ~ "dual",
    n_pos == 0 & n_neg > 0 ~ "repressor"
  )) %>%
  group_by(mor) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

100/sum(mor_TFs_df$n) * mor_TFs_df$n[1]
100/sum(mor_TFs_df$n) * mor_TFs_df$n[2]
100/sum(mor_TFs_df$n) * mor_TFs_df$n[3]

mor_TFs_df$mor <- factor(mor_TFs_df$mor, levels = unique(mor_TFs_df$mor))

p_1.2.2 <- ggplot(mor_TFs_df, aes(x="", y=n, fill=mor)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = c("#c0b372", "#91b57c", "#a36a69"),
                    name = NULL) +
  theme(legend.position = "left",
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"))


pdf("figures/manuscript/p1.2.1.pdf", width = 1.5, height = 2)
p_1.2.1
dev.off()

pdf("figures/manuscript/p1.2.2.pdf", width = 2.8, height = 1.5)
p_1.2.2
dev.off()

### 1.3 coverage compared to other networks
TFs <- map(networks, function(x){table(x$source) %>% as.data.frame() %>% mutate(Freq = 1)})
TF_df <- TFs %>% reduce(full_join, by = "Var1")
TF_df[is.na(TF_df)] <- 0
colnames(TF_df) <- c("source", names(networks))

interactions <- map(networks, function(x){table(paste0(x$source, "_", x$target)) %>% as.data.frame()  %>% mutate(Freq = 1)})

interactions_df <- interactions %>% reduce(full_join, by = "Var1")
interactions_df[is.na(interactions_df)] <- 0

colnames(interactions_df) <- c("interactions", names(networks))

tmp_1 <- interactions_df %>% column_to_rownames("interactions") %>%
  add_column(total =  rowSums(interactions_df %>% column_to_rownames("interactions"))) %>%
  mutate(unique = case_when(
    total == 1 & doro == 1 ~ "doro",
    total == 1 & regnet == 1 ~ "regnet",
    total == 1 & pathComp == 1 ~ "pathComp",
    total == 1 & chea_arch == 1 ~ "chea_arch",
    total == 1 & chea_encode == 1 ~ "chea_encode",
    total == 1 & chea_enrichr == 1 ~ "chea_enrichr",
    total == 1 & chea_GTEx == 1 ~ "chea_GTEx",
    total == 1 & chea_lit == 1 ~ "chea_lit",
    total == 1 & chea_remap == 1 ~ "chea_remap",
    total == 1 & collecTRI == 1 ~ "collecTRI",
    total > 1 ~ "shared"
  ))

plot_df_1 <- data.frame(network = colnames(tmp_1)[1:length(networks)],
                        total = colSums(tmp_1[1:length(networks)])
) %>%
  left_join(table(tmp_1$unique) %>% as.data.frame() %>% rename("network" = "Var1"), by = "network") %>%
  rename("unique" = "Freq") %>%
  replace_na(list(unique = 0)) %>%
  mutate("shared" = total-unique) %>%
  pivot_longer(!network, names_to = "unique", values_to = "count") %>%
  arrange(desc(count)) %>%
  filter(!unique == "total") %>%
  add_column(type = "Interactions") %>%
  mutate(network = recode(network,
                          chea_arch = "ChEA3 ARCHS4",
                          chea_GTEx = "ChEA3 GTEx",
                          chea_enrichr = "ChEA3 Enrichr",
                          regnet = "RegNetwork",
                          doro = "DoRothEA ABC",
                          collecTRI = "CollecTRI",
                          chea_remap = "ChEA3 ReMap",
                          chea_lit = "ChEA3 Literature",
                          chea_encode = "ChEA3 ENCODE",
                          pathComp = "Pathway Commons"))

plot_df_1 %>% group_by(network) %>%
  reframe(total = sum(count),
          individual_counts = count,
          labels = unique) %>%
  mutate(perc = 100/total * individual_counts) %>%
  filter(labels == "unique")

tmp <- TF_df %>% column_to_rownames("source") %>%
  add_column(total =  rowSums(TF_df %>% column_to_rownames("source"))) %>%
  mutate(unique = case_when(
    total == 1 & doro == 1 ~ "doro",
    total == 1 & regnet == 1 ~ "regnet",
    total == 1 & chea_arch == 1 ~ "chea_arch",
    total == 1 & chea_encode == 1 ~ "chea_encode",
    total == 1 & chea_enrichr == 1 ~ "chea_enrichr",
    total == 1 & chea_GTEx == 1 ~ "chea_GTEx",
    total == 1 & pathComp == 1 ~ "pathComp",
    total == 1 & chea_lit == 1 ~ "chea_lit",
    total == 1 & chea_remap == 1 ~ "chea_remap",
    total == 1 & collecTRI == 1 ~ "collecTRI",
    total > 1 ~ "shared"
  ))

plot_df <- data.frame(network = colnames(tmp)[1:length(networks)],
                      total = colSums(tmp[1:length(networks)])
) %>%
  left_join(table(tmp$unique) %>% as.data.frame() %>% rename("network" = "Var1"), by = "network") %>%
  rename("unique" = "Freq") %>%
  replace_na(list(unique = 0)) %>%
  mutate("shared" = total-unique) %>%
  pivot_longer(!network, names_to = "unique", values_to = "count") %>%
  arrange(desc(count)) %>%
  filter(!unique == "total") %>%
  add_column(type = "TFs") %>%
  mutate(network = recode(network,
                          chea_arch = "ChEA3 ARCHS4",
                          chea_GTEx = "ChEA3 GTEx",
                          chea_enrichr = "ChEA3 Enrichr",
                          regnet = "RegNetwork",
                          doro = "DoRothEA ABC",
                          collecTRI = "CollecTRI",
                          chea_remap = "ChEA3 ReMap",
                          chea_lit = "ChEA3 Literature",
                          chea_encode = "ChEA3 ENCODE",
                          pathComp = "Pathway Commons"))

plot_df_final <- rbind(plot_df, plot_df_1)

net_order <- plot_df_final %>%
  filter(type == "TFs") %>%
  group_by(network) %>%
  summarise(all_counts = sum(count)) %>%
  arrange(desc(all_counts)) %>%
  pull(network)

plot_df_final$network <- factor(plot_df_final$network, levels = net_order)
plot_df_final$unique <- factor(plot_df_final$unique, levels = unique(plot_df_final$unique))
plot_df_final$type <- factor(plot_df_final$type, levels = unique(plot_df_final$type))


p_1.3 <- ggplot(plot_df_final, aes(fill=unique, y=count, x=network)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#b2c6e8", "#496bac")) +
  facet_grid(type ~ ., scales='free') +
  #coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines")) +
  ylab("Total Size") +
  xlab("Networks")


pdf("figures/manuscript/p1.3.pdf", width = 2.8, height = 4)
p_1.3
dev.off()

## Figure 2 systematic comparison ---------------------------
### 2.1 Overview benchmark
# load data and run Caret
data(Sonar)
ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                     savePredictions = T)
fit1 <- train(Class ~ .,data=Sonar,method="rf",trControl=ctrl)

# run MLeval
res <- evalm(list(fit1),gnames=c('rf'))

pdf("figures/manuscript/p2.1.1.pdf", width = 6, height = 3)
plot(res$roc)
dev.off()

pdf("figures/manuscript/p2.1.2.pdf", width = 6, height = 3)
plot(res$prg)
dev.off()


### 2.2 Benchmark
bench_agnositc_res <- read_csv("output/benchmark/benchmark_res.csv")
bench_agnositc_res <- bench_agnositc_res %>%
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
order_net <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_agnositc_res$net <- factor(bench_agnositc_res$net, levels = rev(order_net$net))

p_2.2.1 <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  ylim(0.35, 0.79) +
  geom_signif(y_position = c(0.77), xmin = c(10), xmax = c(11),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              margin_top = 0) +
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUROC") +
  xlab("") + coord_flip()

p_2.2.2 <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  ylim(0.4, 0.82) +
  geom_signif(y_position = c(0.81), xmin = c(10), xmax = c(11),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7 ,
              margin_top = 0) +
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUPRC") +
  xlab("") + coord_flip()


pdf("figures/manuscript/p2.2.1.pdf", width = 3.7, height = 2.4)
p_2.2.1
dev.off()

pdf("figures/manuscript/p2.2.2.pdf", width = 3.7, height = 2.4)
p_2.2.2
dev.off()


## Supp 1 Bias ---------------------------
### S1.1 Size difference between TFs in benchmark and background
TFs_bench <- obs_filtered$TF %>% unique()

TFs_collecTRI <- table(collecTRI$source) %>% as.data.frame() %>%
  add_column(network = "collecTRI")
TFs_doro <- table(doro$source) %>% as.data.frame() %>%
  add_column(network = "dorothea")
TFs_regnet <- table(regnet$source) %>% as.data.frame() %>%
  add_column(network = "regnet")


TF_df <- rbind(TFs_collecTRI, TFs_doro, TFs_regnet) %>%
  mutate(benchmark = case_when(
    Var1 %in% TFs_bench ~ "TF in benchmark",
    !Var1 %in% TFs_bench ~ "TF not in benchmark"
  )) %>%
  rename("source" = "Var1") %>%
  rename("nTargets" = "Freq") %>%
  filter(nTargets >= 5) %>%
  mutate(network = recode(network,
                          dorothea = "DoRothEA ABC",
                          collecTRI = "CollecTRI",
                          regnet = "RegNetwork",
                          chea_arch = "ChEA3 ARCHS4"))

TF_df$network <- factor(TF_df$network, levels = c("CollecTRI", "DoRothEA ABC", "RegNetwork"))

stat.test <- TF_df %>%
  group_by(network) %>%
  t_test(nTargets ~ benchmark) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "network", dodge = 0.8) %>%
  mutate(p.adj.signif = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01 & p.adj >= 0.001  ~ "**",
    p.adj < 0.05 & p.adj >= 0.01 ~ "*",
  ))

p_S1.1 <- ggboxplot(TF_df, x = "network", y = "nTargets", fill = "benchmark", outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif",
                     tip.length = 0.01, bracket.nudge.y = -2) +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm")) +
  xlab("TFs") +
  ylab("Number of targets")


pdf("figures/manuscript/pS1.1.pdf", width = 6, height = 2.8)
p_S1.1
dev.off()

### S1.2 Correlation per exp
# correlation per experiment
cor_perExp <- map_dfr(names(act), function(act_i){
  act_df <- act[[act_i]] %>%
    pivot_longer(!...1,
                 names_to = "source",
                 values_to = "act"
    )

  if(str_detect(string = act_i, pattern = "collecTRI")){
    act_df <- act_df %>% left_join(TF_df %>% filter(network == "CollecTRI"), by = "source")
  } else if (str_detect(string = act_i, pattern = "dorothea")){
    act_df <- act_df %>% left_join(TF_df %>% filter(network == "DoRothEA ABC"), by = "source")
  } else if (str_detect(string = act_i, pattern = "regnet")){
    act_df <- act_df %>% left_join(TF_df %>% filter(network == "RegNetwork"), by = "source")
  }

  act_df <- act_df %>%
    filter(!is.na(act))
  map_dfr(unique(act_df$...1), function(exp){
    df_exp <- act_df %>% filter(...1 == exp)
    data.frame(experiment = exp,
               pearson.cor = cor(abs(df_exp$act), df_exp$nTargets, method = "pearson"),
               network = act_i)
  })
})

cor_perExp <- cor_perExp %>%
  mutate(network = recode(network,
                          collecTRI = "CollecTRI",
                          dorothea = "DoRothEA ABC",
                          regnet = "RegNetwork")) %>%
  filter(network %in% c("CollecTRI", "DoRothEA ABC", "RegNetwork"))
# %>% filter(network != "collecTRI_rand")
labels_plot <- cor_perExp %>%
  group_by(network) %>%
  summarise(mean_cor = round(mean(pearson.cor), digits = 2)) %>%
  pull(mean_cor) %>%
  paste0("mean(r) = ", .)
names(labels_plot) <- unique(cor_perExp %>% pull(network))

p_S1.2 <- ggplot(cor_perExp ,
                 aes(x = pearson.cor, fill = network) ) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  xlab("pearson correlation (r)") +
  ylab("Number of experiments") +
  facet_grid(. ~ network, labeller = as_labeller(labels_plot))


pdf("figures/manuscript/pS1.2.pdf", width = 6, height = 2.7)
p_S1.2
dev.off()

## Supp 2 Benchmark per source ---------------------------
bench_signed_source_res <- read_csv("output/benchmark/benchmark_source_res.csv")
bench_signed_source_res <- bench_signed_source_res %>%
  mutate(net = recode(net,
                      ABC = "DoRothEA ABC",
                      collecTRI = "CollecTRI",
                      regnet = "RegNetwork"))

TFs_to_keep <- bench_signed_source_res %>%
  filter(net == "CollecTRI") %>%
  pull(source) %>%
  unique()
bench_signed_source_res <- bench_signed_source_res %>%
  filter(source %in% TFs_to_keep)

signed_source_summarized <- bench_signed_source_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric %in% c("mcauprc", "mcauroc")) %>%
  group_by(net, source, metric) %>%
  summarise(median = median(score)) %>%
  pivot_wider(names_from = metric, values_from = median) %>%
  rename("Network" = "net")

signed_source_summarized$Network <- factor(signed_source_summarized$Network, levels = rev(order_net$net))

p_S2.1 <- ggplot(signed_source_summarized, aes(color=Network, y=mcauprc, x=mcauroc)) +
  geom_point(position=position_jitter(h=0.02,w=0.02), size = 1) +
  geom_vline(xintercept=0.5, linetype="dashed", color = "black", alpha = 0.4) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.4) +
  theme_bw() +
  facet_wrap(~source, strip.position = "top", nrow = 3) +
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("median AUPRC") +
  xlab("median AUROC")  +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing=unit(0.1,'cm')) + theme(legend.position = "bottom")

top_auroc <- signed_source_summarized %>%
  group_by(source) %>%
  filter(mcauroc == max(mcauroc))
top_auroc <- table(top_auroc$Network) %>% as.data.frame() %>% add_column(metric = "meadian auroc")
top_auprc <- signed_source_summarized %>%
  group_by(source) %>%
  filter(mcauprc == max(mcauprc))
top_auprc <- table(top_auprc$Network) %>% as.data.frame() %>% add_column(metric = "meadian auprc")
top_performing_method <- data.frame(network = signed_source_summarized$Network,
                                    source = signed_source_summarized$source,
                                    score = c(signed_source_summarized$mcauroc, signed_source_summarized$mcauprc),
                                    metric = rep(c("mcauroc", "mcauprc"), each = nrow(signed_source_summarized)))

p_S2.2.1 <- top_performing_method %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=network, y=score, fill=network)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  ylim(0.35, 0.81) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("median AUROC") +
  xlab("")

p_S2.2.2 <- top_performing_method %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=network, y=score, fill=network)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  ylim(0.35, 0.81) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("median AUPRC") +
  xlab("")


pdf("figures/manuscript/pS2.1.pdf", width = 4.4, height = 4.2)
p_S2.1
dev.off()

pdf("figures/manuscript/pS2.2.1.pdf", width = 1.5, height = 2.5)
p_S2.2.1
dev.off()

pdf("figures/manuscript/pS2.2.2.pdf", width = 1.5, height = 2.5)
p_S2.2.2
dev.off()


## Supp 1 weights ---------------------------
### S1.1 Weighting strategy
bench_weights_res <- read_csv("output/040722/benchmark/weights_res.csv")
bench_weights_res <- bench_weights_res %>%
  mutate(normalisation = case_when(
    str_detect(net, "raw") ~ "raw",
    str_detect(net, "gene") ~ "per gene",
    str_detect(net, "tf") ~ "per TF",
    str_detect(net, "collecTRI") ~ "raw",
  )) %>%
  mutate(weightingMethod = case_when(
    str_detect(net, "collecTRI") ~ "none",
    str_detect(net, "FIMO") ~ "FIMO",
    str_detect(net, "matRid") ~ "matrixRider"
  ))

order_net <- bench_weights_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_weights_res$net <- factor(bench_weights_res$net, levels = rev(order_net$net))

bench_weights_res_1kb <- bench_weights_res %>%
  filter(!str_detect(net, "10") | str_detect(net, "collecTRI"))

bench_weights_res_10kb <- bench_weights_res %>%
  filter(str_detect(net, "10") | str_detect(net, "collecTRI"))

bench_weights_res

p_S1.1.1 <- bench_weights_res_1kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=normalisation)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUROC") +
  xlab("Weighting method")

p_S1.1.2 <- bench_weights_res_1kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=normalisation)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUPRC") +
  xlab("Weighting method")


p_S1.2.1 <- bench_weights_res_10kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=normalisation)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUROC") +
  xlab("Weighting method")

p_S1.2.2 <- bench_weights_res_10kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=normalisation)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUPRC") +
  xlab("Weighting method")


pdf("figures/manuscript/pS1.1.1.pdf", width = 3.7, height = 2)
p_S1.1.1
dev.off()

pdf("figures/manuscript/pS1.1.2.pdf", width = 3.7, height = 2)
p_S1.1.2
dev.off()

pdf("figures/manuscript/pS1.2.1.pdf", width = 3.7, height = 2)
p_S1.2.1
dev.off()

pdf("figures/manuscript/pS1.2.2.pdf", width = 3.7, height = 2)
p_S1.2.2
dev.off()


## Use weights as filtering
bench_filtered_res <- read_csv("output/040722/benchmark/weights_filtered_res.csv")
bench_filtered_res <- bench_filtered_res %>%
  mutate(quantile = case_when(
    str_detect(net, "raw") ~ "full",
    net == "collecTRI" ~ "full",
    str_detect(net, "_10") ~ "10%",
    str_detect(net, "_20") ~ "20%",
    str_detect(net, "_30") ~ "30%",
  )) %>%
  mutate(weightingMethod = case_when(
    str_detect(net, "collecTRI") ~ "none",
    str_detect(net, "FIMO") ~ "FIMO",
    str_detect(net, "matRid") ~ "matrixRider"
  ))

order_net <- bench_filtered_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_filtered_res$net <- factor(bench_filtered_res$net, levels = rev(order_net$net))
bench_filtered_res$quantile <- factor(bench_filtered_res$quantile, levels = c("full", "10%", "20%", "30%"))

bench_filtered_res_1kb <- bench_filtered_res %>%
  filter(str_detect(net, "matRid1_") | str_detect(net, "FIMO1_") | str_detect(net, "matRid1r") | str_detect(net, "FIMO1r") | str_detect(net, "collecTRI"))

bench_filtered_res_10kb <- bench_filtered_res %>%
  filter(str_detect(net, "matRid10") | str_detect(net, "FIMO10") | str_detect(net, "collecTRI"))


p_S1.3.1 <- bench_filtered_res_1kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=quantile)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUROC") +
  xlab("Weighting method")

p_S1.3.2 <- bench_filtered_res_10kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=quantile)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUPRC") +
  xlab("Weighting method")


p_S1.4.1 <- bench_filtered_res_10kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=quantile)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUROC") +
  xlab("Weighting method")

p_S1.4.2 <- bench_filtered_res_10kb %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=quantile)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUPRC") +
  xlab("Weighting method")


pdf("figures/manuscript/pS1.3.1.pdf", width = 3.6, height = 2)
p_S1.3.1
dev.off()

pdf("figures/manuscript/pS1.3.2.pdf", width = 3.6, height = 2)
p_S1.3.2
dev.off()

pdf("figures/manuscript/pS1.4.1.pdf", width = 3.6, height = 2)
p_S1.4.1
dev.off()

pdf("figures/manuscript/pS1.4.2.pdf", width = 3.6, height = 2)
p_S1.4.2
dev.off()



## Supp Sign decision ---------------------------


## Methods TFs covered in benchmark ---------------------------
TFs_in_bench <- map_dfr(names(networks), function(name_net){
  net <- networks[[name_net]]
  data.frame(network = name_net,
             n_TFs = sum(unique(net$source) %in% obs_filtered$TF))
})

write_csv(TFs_in_bench %>% arrange(network), "figures/manuscript/table_TFs_in_benchmark.csv")
