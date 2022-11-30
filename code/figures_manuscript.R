# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will construct the figures for the final manuscript

library(tidyverse)
library(UpSetR)

## Load networks---------------------------
doro <- read_csv("data/raw/dorothea_ABC.csv") %>%
  select(-...1)
regnet <- read_csv("data/raw/regnetwork.csv")
chea <- read_csv("data/raw/chea3.csv")
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

collecTRI <- read_csv("output/040722/final_dbTF/signed_collecTRI.csv")

networks <- list(doro = doro,
                 regnet = regnet,
                 chea_ = chea,
                 collecTRI = collecTRI)

## Figure 1 coverage ---------------------------
### 1.1 coverage compared to other networks
TFs <- map(networks, function(x){table(x$source) %>% as.data.frame() %>% mutate(Freq = 1)})

TF_df <- TFs %>% reduce(full_join, by = "Var1")
TF_df[is.na(TF_df)] <- 0

colnames(TF_df) <- c("source", names(networks))

TF_upset_p <- UpSetR::upset(TF_df, sets = c(names(networks)), order.by = "freq", text.scale = 1.5)


Genes <- map(networks, function(x){table(x$target) %>% as.data.frame()  %>% mutate(Freq = 1)})

Genes_df <- Genes %>% reduce(full_join, by = "Var1")
Genes_df[is.na(Genes_df)] <- 0

colnames(Genes_df) <- c("targets", names(networks))

Genes_upset_p <- UpSetR::upset(Genes_df, sets = c(names(networks)), order.by = "freq", text.scale = 1.5)


interactions <- map(networks, function(x){table(paste0(x$source, "_", x$target)) %>% as.data.frame()  %>% mutate(Freq = 1)})

interactions_df <- interactions %>% reduce(full_join, by = "Var1")
interactions_df[is.na(interactions_df)] <- 0

colnames(interactions_df) <- c("interactions", names(networks))

interactions_upset_p <- UpSetR::upset(interactions_df, sets = c(names(networks)),
                                      order.by = "freq", text.scale = 1.5)

### 1.2 Benchmark agnostic
bench_agnositc_res <- read_csv("output/040722/benchmark/agnositc_res.csv")
bench_agnositc_res <- bench_agnositc_res %>%
  filter(!net == "collecTRI") %>%
  mutate(net = recode(net, collecTRI_signed = "collecTRI",
                      ABC = "Dorothea ABC", rand = "random network"))
order_net <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_agnositc_res$net <- factor(bench_agnositc_res$net, levels = rev(order_net$net))
p_1.2 <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 12)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUROC") +
  xlab("")

pdf("figures/manuscript/p1.2.pdf", width = 5, height = 5)
p_1.2
dev.off()
## Figure 2 effect of signs ---------------------------
### 2.1 Overview sign decision


### 2.2 Number positive and negative interactions after each step
sign_collecTRI <- read_csv("output/040722/signed_CollecTRI_dbTF_040722.csv")

sign_decision_df <- map_dfr(1:length(unique(sign_collecTRI$decision)), function(i){
  typ <- c("PMID", "keywords", "regulon", "unknown")[1:i]
  typ_net <- sign_collecTRI %>%
    filter(decision %in% typ)

  data.frame(decision = rep(paste(typ, collapse = "\n"), times = 3),
             regulation = c("pos", "neg", "unknown"),
             n = c(sum(typ_net$sign > 0),
                   sum(typ_net$sign < 0),
                   nrow(sign_collecTRI) - nrow(typ_net)))



})
sign_decision_df$decision <- factor(sign_decision_df$decision, levels = unique(sign_decision_df$decision))

p_2.2 <- ggplot(data=sign_decision_df, aes(x=forcats::fct_rev(decision), y=n, fill=regulation)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.85), width = 0.8) +
  theme_minimal() +
  scale_fill_manual(values=c('#A55A5B','#91B57D', '#979797')) +
  xlab("") +
  ylab("number of interactions") +
  theme(text = element_text(size = 10)) + coord_flip()

pdf("figures/manuscript/p2.2.pdf", width = 5, height = 5)
p_2.2
dev.off()

### 2.3 TFs mode of regulation
mor_TFs_df <- collecTRI %>%
  group_by(source) %>%
  summarize(n_pos = sum(weight > 0),
            n_neg = sum(weight < 0)) %>%
  mutate(mor = case_when(
    n_pos > 0 & n_neg == 0 ~ "activators",
    n_pos > 0 & n_neg > 0 ~ "dual",
    n_pos == 0 & n_neg > 0 ~ "repressors"
  )) %>%
  group_by(mor) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

mor_TFs_df$mor <- factor(mor_TFs_df$mor, levels = unique(mor_TFs_df$mor))

p_2.3 <- ggplot(mor_TFs_df, aes(x="", y=n, fill=mor)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = c("#C0B372", "#91B57D", "#A55A5B"),
                    name = NULL) +
  ggtitle("TFs mode of regulation (MoR)") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 12),
        legend.position = "bottom",
        text = element_text(size = 10))

pdf("figures/manuscript/p2.3.pdf", width = 5, height = 5)
p_2.3
dev.off()
