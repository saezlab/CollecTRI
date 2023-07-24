# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will construct the figures for the final manuscript

library(tidyverse)
library(ggsignif)
library(rstatix)
library(ggpubr)
library(caret)
library(RColorBrewer)
library(pheatmap)

## Load data---------------------------
### networks
doro <- read_csv("data/networks/filtered_dorothea_ABC.csv")
doro_ABCD <- read_csv("data/networks/filtered_dorothea_ABCD.csv")
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

collecTRI <- read_csv("output/CollecTRI/CollecTRI_GRN.csv")

networks <- list(doro = doro,
                 doro_ABCD = doro_ABCD,
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

# Percentage of activation/repression of all TF-gene links
100/sum(sign_comp_df$Freq) * sign_comp_df$Freq[1] #repression
100/sum(sign_comp_df$Freq) * sign_comp_df$Freq[2] #activation

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

# Percentage of mode of regulation of TFs (activators, repressors, dual)
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

collecTRI$sign_decision %>% table()
collecTRI$weight %>% table()

nrow(collecTRI)
length(unique(collecTRI$source))

### 1.3 coverage compared to other networks
TFs <- map(networks, function(x){table(x$source) %>% as.data.frame() %>% mutate(Freq = 1)})
TF_df <- TFs %>% reduce(full_join, by = "Var1")
TF_df[is.na(TF_df)] <- 0
colnames(TF_df) <- c("source", names(networks))
100/nrow(TF_df) * sum(rowSums(TF_df %>% column_to_rownames("source")) > 1)

interactions <- map(networks, function(x){table(paste0(x$source, "_", x$target)) %>% as.data.frame()  %>% mutate(Freq = 1)})

interactions_df <- interactions %>% reduce(full_join, by = "Var1")
interactions_df[is.na(interactions_df)] <- 0

colnames(interactions_df) <- c("interactions", names(networks))

tmp_1 <- interactions_df %>% column_to_rownames("interactions") %>%
  add_column(total =  rowSums(interactions_df %>% column_to_rownames("interactions"))) %>%
  mutate(unique = case_when(
    total == 1 & doro == 1 ~ "doro",
    total == 1 & doro_ABCD == 1 ~ "doro_ABCD",
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
                          doro_ABCD = "DoRothEA ABCD",
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

plot_df_1 %>% group_by(network) %>%
  reframe(total = sum(count),
          individual_counts = count,
          labels = unique) %>%
  mutate(perc = 100/total * individual_counts) %>%
  filter(labels == "unique") %>%
  pull(perc) %>%
  mean()

tmp <- TF_df %>% column_to_rownames("source") %>%
  add_column(total =  rowSums(TF_df %>% column_to_rownames("source"))) %>%
  mutate(unique = case_when(
    total == 1 & doro == 1 ~ "doro",
    total == 1 & doro_ABCD == 1 ~ "doro_ABCD",
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
                          doro_ABCD = "DoRothEA ABCD",
                          collecTRI = "CollecTRI",
                          chea_remap = "ChEA3 ReMap",
                          chea_lit = "ChEA3 Literature",
                          chea_encode = "ChEA3 ENCODE",
                          pathComp = "Pathway Commons"))

plot_df %>% filter(unique == "unique")

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

plot_df_final %>%
  group_by(network, type) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total))


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


pdf("figures/manuscript/p1.3.pdf", width = 2.8, height = 3.4)
p_1.3
dev.off()

## Figure 2 systematic comparison ---------------------------
### 2.1 Overview benchmark

### 2.2 Benchmark
bench_agnositc_res <- read_csv("output/benchmark/benchmark_res.csv")
bench_agnositc_res <- bench_agnositc_res %>%
  mutate(net = recode(net,
                      chea3_archs4 = "ChEA3 ARCHS4",
                      chea3_GTEx = "ChEA3 GTEx",
                      chea3_enrich = "ChEA3 Enrichr",
                      regnet = "RegNetwork",
                      ABC = "DoRothEA ABC",
                      ABCD = "DoRothEA ABCD",
                      collecTRI = "CollecTRI",
                      chea3_remap = "ChEA3 ReMap",
                      chea3_lit = "ChEA3 Literature",
                      chea3_encode = "ChEA3 ENCODE",
                      rand = "shuffled CollecTRI",
                      pathComp = "Pathway Commons"))
order_net <- bench_agnositc_res %>%
  filter(method == "ulm_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_agnositc_res$net <- factor(bench_agnositc_res$net, levels = rev(order_net$net))

p_2.2.1 <- bench_agnositc_res %>%
  filter(method == "ulm_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  ylim(0.35, 0.79) +
  geom_signif(y_position = c(0.79), xmin = c(11), xmax = c(12),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black", alpha = 0.7) +
  ylab("AUROC") +
  xlab("") + coord_flip()

p_2.2.2 <- bench_agnositc_res %>%
  filter(method == "ulm_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  ylim(0.4, 0.82) +
  geom_signif(y_position = c(0.82), xmin = c(11), xmax = c(12),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7 ,
              #significance level taken from statistics.R
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

## Supp 1 Resources CollecTRI source ---------------------------
raw.file <- "data/CollecTRI_source.tsv"

collecTRI.raw <- read.table(raw.file,
                            sep = "\t",
                            header = TRUE) # load raw resource list
resources <- colnames(collecTRI.raw)[str_detect(colnames(collecTRI.raw), "PMID")] %>%
  str_remove("..PMID") %>%
  str_remove("X.")

egdes_resource <- map(resources, function(resource){
  res_df <- collecTRI.raw[,paste0("X.", resource, "..present")]
  collecTRI.raw$TF.TG[res_df != ""]
})

unique_int <- lapply(1:length(egdes_resource), function(n) length(setdiff(egdes_resource[[n]], unlist(egdes_resource[-n])))) %>% unlist()

df_TRI <- data.frame(resource = rep(resources, times = 2),
           count = c(unique_int, map_dbl(egdes_resource, length) - unique_int),
           type = rep(c("unique", "shared"), each = length(resources)),
           class = "TRI"
           )

tfs_resource <- map(resources, function(resource){
  res_df <- collecTRI.raw[,paste0("X.", resource, "..present")]
  edges <- collecTRI.raw$TF.TG[res_df != ""]
  map_chr(str_split(edges, ":"), 1) %>% unique()
})

unique_tfs <- lapply(1:length(tfs_resource), function(n) length(setdiff(tfs_resource[[n]], unlist(tfs_resource[-n])))) %>% unlist()

df_TFs <- data.frame(resource = rep(resources, times = 2),
                     count = c(unique_tfs, map_dbl(tfs_resource, length) - unique_tfs),
                     type = rep(c("unique", "shared"), each = length(resources)),
                     class = "TF"
)

tgs_resource <- map(resources, function(resource){
  res_df <- collecTRI.raw[,paste0("X.", resource, "..present")]
  edges <- collecTRI.raw$TF.TG[res_df != ""]
  map_chr(str_split(edges, ":"), 2) %>% unique()
})

unique_tgs <- lapply(1:length(tgs_resource), function(n) length(setdiff(tgs_resource[[n]], unlist(tgs_resource[-n])))) %>% unlist()

df_TGs <- data.frame(resource = rep(resources, times = 2),
                     count = c(unique_tgs, map_dbl(tgs_resource, length) - unique_tgs),
                     type = rep(c("unique", "shared"), each = length(resources)),
                     class = "TG"
)

df_resource <- rbind(df_TRI, df_TFs, df_TGs)

df_resource$resource <- factor(df_resource$resource, levels = df_TFs %>% arrange(desc(count)) %>% pull(resource) %>% unique())

p_S1 <- ggplot(df_resource, aes(fill=type, y=count, x=resource)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#b2c6e8", "#496bac")) +
  facet_grid(class ~ ., scales='free_y') +
  #coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        legend.position = "right",
        panel.spacing = unit(1, "lines")) +
  ylab("Total Size") +
  xlab("Resource")

pdf("figures/manuscript/pS1.pdf", width = 6, height = 5)
p_S1
dev.off()

## Supp new -------------
TFs <- map(networks, function(net) unique(net$source)) %>%
  unlist() %>%
  unique()

jaccard_df <- map_dfr(TFs, function(TF){
  regulons <- map(networks, function(net){
      net %>%
        filter(source == TF) %>%
        mutate(edge = paste(source, target, sep = ":")) %>%
        pull(edge)
    })

  dist <- unlist(lapply(combn(regulons, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))

  jaccard_df_full <- cbind(t(combn(names(regulons),2)), dist) %>% as.data.frame()

  jaccard_df_full$dist <- as.numeric(jaccard_df_full$dist)
  jaccard_df_full$TF <- TF

  jaccard_df_full
})

mean_jaccard <- jaccard_df %>%
  group_by(V1, V2) %>%
  summarise(mean_jaccard = mean(dist, na.rm = T))

ord <- table(mean_jaccard$V1) %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  pull(Var1)

ord2 <- table(mean_jaccard$V2) %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  pull(Var1)
mean_jaccard$V1 <- factor( mean_jaccard$V1 , levels = ord)
mean_jaccard$V2 <- factor( mean_jaccard$V2 , levels = ord2)

mean_jaccard <- mean_jaccard %>%
  mutate(V1 = recode(V1,
                      chea_arch = "ChEA3 ARCHS4",
                      chea_GTEx = "ChEA3 GTEx",
                      chea_enrichr = "ChEA3 Enrichr",
                      regnet = "RegNetwork",
                      doro = "DoRothEA ABC",
                      doro_ABCD = "DoRothEA ABCD",
                      collecTRI = "CollecTRI",
                      chea_remap = "ChEA3 ReMap",
                      chea_lit = "ChEA3 Literature",
                      chea_encode = "ChEA3 ENCODE",
                      pathComp = "Pathway Commons")) %>%
  mutate(V2 = recode(V2,
                     chea_arch = "ChEA3 ARCHS4",
                     chea_GTEx = "ChEA3 GTEx",
                     chea_enrichr = "ChEA3 Enrichr",
                     regnet = "RegNetwork",
                     doro = "DoRothEA ABC",
                     doro_ABCD = "DoRothEA ABCD",
                     collecTRI = "CollecTRI",
                     chea_remap = "ChEA3 ReMap",
                     chea_lit = "ChEA3 Literature",
                     chea_encode = "ChEA3 ENCODE",
                     pathComp = "Pathway Commons"))

pS2_new <- ggplot(mean_jaccard, aes(V1, V2, fill= mean_jaccard)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", name = "Mean Jaccard") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("") +
  xlab("")

pdf("figures/manuscript/pS2_jaccard.pdf", width = 6, height = 5)
pS2_new
dev.off()


## Supp 2 Sign ---------------------------
# Add information about sign decision
bench_sign_collecTRI <- read.csv("output/benchmark/benchmark_sign_res.csv")
bench_sign_collecTRI_prior <- bench_sign_collecTRI %>%
  filter(net %in% c("collecTRI_agnostic", "collecTRI_PMID", "collecTRI_TF", "collecTRI_PMID_TF"))
bench_sign_collecTRI_regulon <- bench_sign_collecTRI %>%
  filter(net %in% c("collecTRI_PMID", "collecTRI", "collecTRI_PMID_regulon_TF"))
bench_sign_collecTRI_repression <- bench_sign_collecTRI %>%
  filter(net %in% c("collecTRI_repression", "collecTRI"))

# 2.1 compare prior from PMIDs and TF classificications
bench_sign_collecTRI_prior <- bench_sign_collecTRI_prior %>%
  mutate(net = recode(net,
                      collecTRI_agnostic = "default\nactivation",
                      collecTRI_PMID = "PMID\ndefault\nactivation",
                      collecTRI_TF = "TF role\ndefault\nactivation",
                      collecTRI_PMID_TF = "PMID\nTF role\ndefault\nactivation"))
bench_sign_collecTRI_prior$net <- factor(bench_sign_collecTRI_prior$net, levels = c("default\nactivation",
                                                                                    "PMID\ndefault\nactivation",
                                                                                    "TF role\ndefault\nactivation",
                                                                                    "PMID\nTF role\ndefault\nactivation"))

p_S2.1.1 <- bench_sign_collecTRI_prior %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  geom_signif(y_position = c(0.775), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  ylim(0.6, 0.785) +
  theme(legend.position="none",
        text = element_text(size = 9),
        axis.text.x = element_text(hjust=0.5)) +
  ylab("AUROC") +
  xlab("")

p_S2.1.2 <- bench_sign_collecTRI_prior %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  geom_signif(y_position = c(0.81), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  ylim(0.65, 0.815) +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  ylab("AUPRC") +
  xlab("")


pdf("figures/manuscript/pS2.1.1.pdf", width = 3, height = 2)
p_S2.1.1
dev.off()

pdf("figures/manuscript/pS2.1.2.pdf", width = 3, height = 2)
p_S2.1.2
dev.off()

# 2.2 compare effect of regulon
bench_sign_collecTRI_regulon <- bench_sign_collecTRI_regulon %>%
  mutate(net = recode(net,
                      collecTRI_PMID = "PMID\ndefault\nactivation",
                      collecTRI = "PMID\nRegulon\ndefault\nactivation",
                      collecTRI_PMID_regulon_TF = "PMID\nRegulon\nTF role\ndefault\nactivation"))
bench_sign_collecTRI_regulon$net <- factor(bench_sign_collecTRI_regulon$net, levels = c("PMID\ndefault\nactivation",
                                                                                        "PMID\nRegulon\ndefault\nactivation",
                                                                                        "PMID\nRegulon\nTF role\ndefault\nactivation"))

p_S2.2.1 <- bench_sign_collecTRI_regulon %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  geom_signif(y_position = c(0.775), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  ylim(0.67, 0.785) +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  ylab("AUROC") +
  xlab("")

p_S2.2.2 <- bench_sign_collecTRI_regulon %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  geom_signif(y_position = c(0.825), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  ylim(0.7, 0.835) +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  ylab("AUPRC") +
  xlab("")


pdf("figures/manuscript/pS2.2.1.pdf", width = 2, height = 1.9)
p_S2.2.1
dev.off()

pdf("figures/manuscript/pS2.2.2.pdf", width = 2, height = 1.9)
p_S2.2.2
dev.off()

# 2.3 default activation versus default repression
bench_sign_collecTRI_repression <- bench_sign_collecTRI_repression %>%
  mutate(net = recode(net,
                      collecTRI = "PMID\nRegulon\ndefault\nactivation",
                      collecTRI_repression = "PMID\nRegulon\ndefault\nrepression"))
bench_sign_collecTRI_repression$net <- factor(bench_sign_collecTRI_repression$net, levels = c("PMID\nRegulon\ndefault\nactivation",
                                                                                           "PMID\nRegulon\ndefault\nrepression"))

p_S2.3.1 <- bench_sign_collecTRI_repression %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  geom_signif(y_position = c(0.79), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  ylim(0.49, 0.805) +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  ylab("AUROC") +
  xlab("")

p_S2.3.2 <- bench_sign_collecTRI_repression %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  geom_signif(y_position = c(0.85), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              #significance level taken from statistics.R
              margin_top = 0) +
  theme_minimal() +
  ylim(0.54, 0.87) +
  theme(legend.position="none",
        text = element_text(size = 9)) +
  ylab("AUPRC") +
  xlab("")


pdf("figures/manuscript/pS2.3.1.pdf", width = 1.8, height = 1.8)
p_S2.3.1
dev.off()

pdf("figures/manuscript/pS2.3.2.pdf", width = 1.8, height = 1.8)
p_S2.3.2
dev.off()

# 2.4 overview final sign assignment
sign_collecTRI <- read.csv("output/CollecTRI/CollecTRI_GRN.csv")
sign_collecTRI <- sign_collecTRI %>%
  mutate(sign_decision = recode(sign_decision,
                                "default activation" = "default\nactivation"))

decision_df <- sign_collecTRI %>%
  group_by(sign_decision, weight) %>%
  summarise(total = n()) %>%
  mutate(weight = recode(weight,
                         "1" = "activation",
                         "-1" = "repression")) %>%
  arrange(desc(total))
decision_df <- rbind(decision_df, data.frame(sign_decision = "default\nactivation", weight = "repression", total = 0))
decision_df$sign.decision <- factor(decision_df$sign_decision, levels = unique(decision_df$sign_decision))

p_S2.4 <- ggplot(data=decision_df, aes(x=sign.decision, y=total, fill=weight)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.85), width = 0.8) +
  theme_minimal() +
  scale_fill_manual(values=c('#91B57D','#A55A5B')) +
  xlab("Decision source") +
  ylab("Number of interactions") +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "bottom") +
  guides(fill=guide_legend(title=""))

decision_df %>%
  group_by(sign.decision) %>%
  summarise(all = sum(total))

pdf("figures/manuscript/pS2.4.pdf", width = 2.3, height = 2.9)
p_S2.4
dev.off()

## Supp 3 TF with at least five targets ---------------------------
nTF_df <- map_dfr(names(networks), function(net_i){
  net <- networks[[net_i]]

  data.frame(net = net_i,
             nTFs = c(table(net$source) %>% length(), sum(table(net$source) >= 5)),
             TF = c("all", "at least 5 targets"))
})

nTF_df <- nTF_df %>%
  mutate(net = recode(net,
                      chea_archs4 = "ChEA3 ARCHS4",
                      chea_GTEx = "ChEA3 GTEx",
                      chea_enrich = "ChEA3 Enrichr",
                      regnet = "RegNetwork",
                      doro_ABCD = "DoRothEA ABCD",
                      doro = "DoRothEA ABC",
                      collecTRI = "CollecTRI",
                      chea_remap = "ChEA3 ReMap",
                      chea_lit = "ChEA3 Literature",
                      chea_encode = "ChEA3 ENCODE",
                      rand = "shuffled CollecTRI",
                      pathComp = "Pathway Commons"))

order_net <- nTF_df %>%
  filter(TF == "all") %>%
  arrange(desc(nTFs)) %>%
  pull(net)

nTF_df$net <- factor(nTF_df$net, levels=order_net)
nTF_df$TF <- factor(nTF_df$TF, levels=c("all", "at least 5 targets"))


p_S3 <- ggplot(data=nTF_df, aes(x=net, y=nTFs, fill=TF)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() +
  xlab("") +
  ylab("Number of TFs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  guides(fill=guide_legend(title=""))

pdf("figures/manuscript/pS3.pdf", width = 6, height = 3.4)
p_S3
dev.off()

## Supp 4 Bias ---------------------------
### S4.1 Size difference between TFs in benchmark and background
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
                          regnet = "RegNetwork"))

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

p_S4.1 <- ggboxplot(TF_df, x = "network", y = "nTargets", fill = "benchmark", outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif",
                     tip.length = 0.01, bracket.nudge.y = -2) +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm")) +
  xlab("TFs") +
  ylab("Number of targets")


pdf("figures/manuscript/pS4.1.pdf", width = 6, height = 2.8)
p_S4.1
dev.off()

### S4.2 Correlation per exp
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

labels_plot

p_S4.2 <- ggplot(cor_perExp ,
                 aes(x = pearson.cor, fill = network) ) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "right") +
  xlab("pearson correlation (r)") +
  ylab("Number of experiments") +
  facet_grid(. ~ network, labeller = as_labeller(labels_plot))


pdf("figures/manuscript/pS4.2.pdf", width = 6, height = 2.7)
p_S4.2
dev.off()

## Supp 5 Benchmark per source ---------------------------
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


signed_source_summarized$Network <- factor(signed_source_summarized$Network, levels = unique(signed_source_summarized$Network))

p_S5.1 <- ggplot(signed_source_summarized, aes(color=Network, y=mcauprc, x=mcauroc)) +
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

p_S5.2.1 <- top_performing_method %>%
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

p_S5.2.2 <- top_performing_method %>%
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


pdf("figures/manuscript/pS5.1.pdf", width = 4.4, height = 4.2)
p_S5.1
dev.off()

pdf("figures/manuscript/pS5.2.1.pdf", width = 1.5, height = 2.5)
p_S5.2.1
dev.off()

pdf("figures/manuscript/pS5.2.2.pdf", width = 1.5, height = 2.5)
p_S5.2.2
dev.off()


## Supp 6 weights ---------------------------
### S6.1 Correlation between weighting strategy
corr_matrix <- read_csv("output/weighted_networks/correlation_matrix.csv")%>%
  as.data.frame()
rownames(corr_matrix) <- colnames(corr_matrix)
corr_matrix <- corr_matrix[rownames(corr_matrix) != "unweighted", colnames(corr_matrix) != "unweighted"]

annotation_df <- data.frame(method = map_chr(str_split(colnames(corr_matrix), "_"), 1),
           `promoter length` = map_chr(str_split(colnames(corr_matrix), "_"), 2),
           normalisation = map_chr(str_split(colnames(corr_matrix), "_"), 3)
           ) %>%
  mutate(normalisation = recode(normalisation,
                                gene = "per gene",
                                tf = "per TF",
                                raw = "none"))

rownames(annotation_df) <- colnames(corr_matrix)

cor_heat <- pheatmap::pheatmap(corr_matrix, cluster_rows = T,
                               cluster_cols = T, silent=T)
idxRows <- cor_heat$tree_row$order
idxCols <- cor_heat$tree_col$order
corr_matrix <- corr_matrix[idxRows,idxCols]

# Plot
color <- colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = 'Blues')))(100)

mat_heat <- pheatmap::pheatmap(corr_matrix, color = color,
                               display_numbers=F, number_color='black', border_color=NA,
                               na_col=NA, cellwidth = 30, cellheight = 30,
                               legend=T, fontsize = 10,
                               show_rownames = T, show_colnames = F,treeheight_col = 0,
                               silent=T, annotation_col = annotation_df)
p_S6.1 <- ggplotify::as.ggplot(mat_heat) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

pdf("figures/manuscript/pS6.1.pdf", width = 10, height = 10)
p_S6.1
dev.off()

range(corr_matrix)
### S6.2 Weighting strategy
bench_weights_res <- read_csv("output/benchmark/benchmark_weights_res.csv")
bench_weights_res <- bench_weights_res %>%
  filter(net %in% c("collecTRI", "matRid1raw")) %>%
  mutate(net = recode(net,
                      collecTRI = 'CollecTRI',
                      matRid1raw = 'weighted CollecTRI'))

order_net <- bench_weights_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_weights_res$net <- factor(bench_weights_res$net, levels = rev(order_net$net))

t_auroc <- t.test(bench_weights_res  %>%
         filter(method == "consensus_estimate") %>%
         filter(metric == "mcauroc") %>%
         filter(net == "CollecTRI") %>%
         pull(score), bench_weights_res  %>%
         filter(method == "consensus_estimate") %>%
         filter(metric == "mcauroc") %>%
         filter(net == "weighted CollecTRI") %>%
         pull(score))$p.value

t_auprc <- t.test(bench_weights_res  %>%
         filter(method == "consensus_estimate") %>%
         filter(metric == "mcauprc") %>%
         filter(net == "CollecTRI") %>%
         pull(score), bench_weights_res  %>%
         filter(method == "consensus_estimate") %>%
         filter(metric == "mcauprc") %>%
         filter(net == "weighted CollecTRI") %>%
         pull(score))$p.value

p.adjust(c(t_auroc, t_auprc), method = "BH")

p_S6.2.1 <- bench_weights_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUROC") +
  xlab("")

p_S6.2.2 <- bench_weights_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=net, y=score, fill=net)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm")) +
  ylab("AUPRC") +
  xlab("")

pdf("figures/manuscript/pS6.2.1.pdf", width = 2.5, height = 2.5)
p_S6.2.1
dev.off()

pdf("figures/manuscript/pS6.2.2.pdf", width = 2.5, height = 2.5)
p_S6.2.2
dev.off()

## Use weights as filtering
bench_filtered_res <- read_csv("output/benchmark/benchmark_weights_filtered_res.csv")
bench_filtered_res <- bench_filtered_res %>%
  filter(str_detect(net, "matRid1_") | str_detect(net, "matRid1r") | str_detect(net, "collecTRI")) %>%
  mutate(quantile = case_when(
    str_detect(net, "raw") ~ "0%",
    net == "collecTRI" ~ "0%",
    str_detect(net, "_10") ~ "10%",
    str_detect(net, "_20") ~ "20%",
    str_detect(net, "_30") ~ "30%",
  )) %>%
  mutate(weightingMethod = case_when(
    str_detect(net, "collecTRI") ~ "random",
    str_detect(net, "matRid") ~ "weights"
  ))


order_net <- bench_filtered_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  group_by(net) %>%
  summarize(mean = median(score)) %>%
  arrange(desc(mean))
bench_filtered_res$net <- factor(bench_filtered_res$net, levels = rev(order_net$net))
bench_filtered_res$quantile <- factor(bench_filtered_res$quantile, levels = c("0%", "10%", "20%", "30%"))

p_S6.3.1 <- bench_filtered_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauroc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=quantile)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  guides(fill=guide_legend(title="Quantile of\nremoved edges")) +
  ylab("AUROC") +
  xlab("Basis for the removal of edges from CollecTRI")

p_S6.3.2 <- bench_filtered_res %>%
  filter(method == "consensus_estimate") %>%
  filter(metric == "mcauprc") %>%
  ggplot(aes(x=weightingMethod, y=score, fill=quantile)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 9),
        legend.key.size = unit(0.3, "cm")) +
  guides(fill=guide_legend(title="Quantile of\nremoved edges")) +
  ylab("AUPRC") +
  xlab("Basis for the removal of edges from CollecTRI")


pdf("figures/manuscript/pS6.3.1.pdf", width = 4, height = 2)
p_S6.3.1
dev.off()

pdf("figures/manuscript/pS6.3.2.pdf", width = 4, height = 2)
p_S6.3.2
dev.off()



## Methods TFs covered in benchmark ---------------------------
TFs_in_bench <- map_dfr(names(networks), function(name_net){
  net <- networks[[name_net]]
  data.frame(network = name_net,
             n_TFs = sum(unique(net$source) %in% obs_filtered$TF))
})

write_csv(TFs_in_bench %>% arrange(network), "figures/manuscript/table_TFs_in_benchmark.csv")
