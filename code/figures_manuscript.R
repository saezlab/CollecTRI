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

collecTRI <- read_csv("output/040722/02_signed_networks/strict_signed_CollecTRI.csv")

networks <- list(doro = doro,
                 regnet = regnet,
                 chea_arch = chea_arch,
                 chea_encode = chea_encode,
                 chea_enrichr = chea_enrichr,
                 chea_GTEx = chea_GTEx,
                 chea_lit = chea_lit,
                 chea_remap = chea_remap,
                 collecTRI = collecTRI)

#networks <- list(doro = doro,
#                 regnet = regnet,
#                 chea = chea,
#                 collecTRI = collecTRI)

## Figure 1 coverage ---------------------------
### 1.2 coverage compared to other networks
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
  mutate(count = log(count))

tmp <- TF_df %>% column_to_rownames("source") %>%
  add_column(total =  rowSums(TF_df %>% column_to_rownames("source"))) %>%
  mutate(unique = case_when(
    total == 1 & doro == 1 ~ "doro",
    total == 1 & regnet == 1 ~ "regnet",
    total == 1 & chea_arch == 1 ~ "chea_arch",
    total == 1 & chea_encode == 1 ~ "chea_encode",
    total == 1 & chea_enrichr == 1 ~ "chea_enrichr",
    total == 1 & chea_GTEx == 1 ~ "chea_GTEx",
    total == 1 & chea_lit == 1 ~ "chea_lit",
    total == 1 & chea_remap == 1 ~ "chea_remap",
    total == 1 & collecTRI == 1 ~ "collecTRI",
    total > 1 ~ "shared"
  ))

#tmp <- TF_df %>% column_to_rownames("source") %>%
#  add_column(total =  rowSums(TF_df %>% column_to_rownames("source"))) %>%
#  mutate(unique = case_when(
#    total == 1 & doro == 1 ~ "doro",
#    total == 1 & regnet == 1 ~ "regnet",
#    total == 1 & chea == 1 ~ "chea",
#    total == 1 & collecTRI == 1 ~ "collecTRI",
#    total > 1 ~ "shared"
#  ))

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
  add_column(type = "TFs")

plot_df_final <- rbind(plot_df, plot_df_1)

plot_df_final$network <- factor(plot_df_final$network, levels = unique(plot_df_final$network))
plot_df_final$unique <- factor(plot_df_final$unique, levels = unique(plot_df_final$unique))
plot_df_final$type <- factor(plot_df_final$type, levels = unique(plot_df_final$type))

p_1.2 <- ggplot(plot_df_final, aes(fill=unique, y=count, x=network)) +
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

pdf("figures/manuscript/p1.2.pdf", width = 2.7, height = 5.5)
p_1.2
dev.off()

### 1.3 Benchmark agnostic
bench_agnositc_res <- read_csv("output/040722/benchmark/agnositc_res.csv")
bench_agnositc_res <- bench_agnositc_res %>%
  filter(!net == "collecTRI_signed") %>%
  mutate(net = recode(net, ABC = "Dorothea ABC", rand = "random network"))
order_net <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_agnositc_res$net <- factor(bench_agnositc_res$net, levels = rev(order_net$net))
p_1.3.1 <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUROC") +
  xlab("")

p_1.3.2 <- bench_agnositc_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUPRC") +
  xlab("")

pdf("figures/manuscript/p1.3.1.pdf", width = 3.5, height = 2)
p_1.3.1
dev.off()

pdf("figures/manuscript/p1.3.2.pdf", width = 3.5, height = 2)
p_1.3.2
dev.off()
## Figure 2 effect of signs ---------------------------
### 2.1 Overview sign decision


### 2.2 Number positive and negative interactions after each step
sign_collecTRI <- read_csv("output/040722/02_signed_networks/strict_signed_CollecTRI.csv")

sign_decision_df <- map_dfr(1:length(unique(sign_collecTRI$decision)), function(i){
  typ <- c("PMID", "keywords", "regulon", "unknown")[1:i]
  typ_net <- sign_collecTRI %>%
    filter(decision %in% typ)

  data.frame(decision = rep(paste(typ, collapse = "\n"), times = 3),
             regulation = c("activation", "inhibition", "unknown"),
             n = c(sum(typ_net$weight > 0),
                   sum(typ_net$weight < 0),
                   nrow(sign_collecTRI) - nrow(typ_net)))



})
sign_decision_df <- sign_decision_df %>% mutate(decision = recode(decision,
                                              "PMID\nkeywords" = "PMID\nTF role",
                                              "PMID\nkeywords\nregulon" = "PMID\nTF role\nRegulon",
                                              "PMID\nkeywords\nregulon\nunknown" = "PMID\nTF role\nRegulon\nUnknown"))
sign_decision_df$decision <- factor(sign_decision_df$decision, levels = unique(sign_decision_df$decision))
sign_decision_df$regulation <- factor(sign_decision_df$regulation, levels = c("inhibition", "activation", "unknown"))

p_2.2 <- ggplot(data=sign_decision_df, aes(x=forcats::fct_rev(decision), y=n, fill=regulation)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.85), width = 0.8) +
  theme_minimal() +
  scale_fill_manual(values=c('#A55A5B','#91B57D', '#979797')) +
  xlab("") +
  ylab("Number of interactions") +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) + coord_flip() +
  guides(fill=guide_legend(title=""))

pdf("figures/manuscript/p2.2.pdf", width = 2.6, height = 2.3)
p_2.2 + theme(legend.position = "none")
dev.off()

pdf("figures/manuscript/p2.2_legend.pdf", width = 3, height = 2.3)
p_2.2  + theme(legend.position = "bottom")
dev.off()

### 2.3 TFs mode of regulation
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

mor_TFs_df$mor <- factor(mor_TFs_df$mor, levels = unique(mor_TFs_df$mor))

p_2.3 <- ggplot(mor_TFs_df, aes(x="", y=n, fill=mor)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = c("#C0B372", "#91B57D", "#A55A5B"),
                    name = NULL) +
  ggtitle("TFs role in regulation") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 9),
        legend.position = "left",
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"))

pdf("figures/manuscript/p2.3.pdf", width = 2.4, height = 1.5)
p_2.3
dev.off()

### 2.4 Sign benchmark
bench_signed_res <- read_csv("output/040722/benchmark/signed_res.csv")
bench_signed_res <- bench_signed_res %>%
  mutate(net = recode(net, ABC = "Dorothea ABC", rand = "random network"))
order_net <- bench_signed_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_signed_res$net <- factor(bench_signed_res$net, levels = rev(order_net$net))
p_2.4.1 <- bench_signed_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUROC") +
  xlab("")

p_2.4.2 <- bench_signed_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUPRC") +
  xlab("")

pdf("figures/manuscript/p2.4.1.pdf", width = 2, height = 2.3)
p_2.4.1
dev.off()

pdf("figures/manuscript/p2.4.2.pdf", width = 2, height = 2.3)
p_2.4.2
dev.off()

### 2.5 Sign benchmark per source
bench_signed_source_res <- read_csv("output/040722/benchmark/signed_source_res.csv")
bench_signed_source_res <- bench_signed_source_res %>%
  mutate(net = recode(net, ABC = "Dorothea ABC", rand = "random network"))
signed_source_summarized <- bench_signed_source_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric %in% c("mcauprc", "mcauroc")) %>%
  group_by(net, source, metric) %>%
  summarise(mean = mean(score)) %>%
  pivot_wider(names_from = metric, values_from = mean) %>%
  rename("Network" = "net")

signed_source_summarized$Network <- factor(signed_source_summarized$Network, levels = rev(order_net$net))

p_2.5 <- ggplot(signed_source_summarized, aes(color=Network, y=mcauprc, x=mcauroc)) +
  geom_point(position=position_jitter(h=0.01,w=0.01), size = 0.6) +
  geom_vline(xintercept=0.5, linetype="dashed", color = "black", alpha = 0.4) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.4) +
  theme_bw() +
  facet_wrap(~source, strip.position = "top") +
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("AUPRC") +
  xlab("AUROC")  +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.spacing=unit(0.1,'cm'))

pdf("figures/manuscript/p_2.5.pdf", width = 4.2, height = 3.5)
p_2.5 + theme(legend.position = "none")
dev.off()

pdf("figures/manuscript/p_2.5_legend.pdf", width = 5.2, height = 3.5)
p_2.5 + theme(legend.position = "bottom")
dev.off()
