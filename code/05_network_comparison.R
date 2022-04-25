library(tidyverse)
library(UpSetR)
library(ggplotify)
library(patchwork)

# Functions ------------------------------------------------
get_auc_df <- function(df, .type){
  .type <- enquo(.type)
  df %>%
    select(set_name, set_name, !!.type) %>%
    mutate(!!.type := map(!!.type, function(df){
      df %>%
        group_by(run) %>%
        summarize(raw_auc = unique(raw_auc)) %>%
        pull(raw_auc)
    })) %>%
    unnest(cols = c(!!.type)) %>%
    filter(!(set_name=='weighted' &
               (set_name %in% c('aucell','ora','norm_fgsea','fgsea','gsva')))) %>%
    filter(!(set_name=='unweighted' &
               (!set_name %in% c('aucell','ora','norm_fgsea','fgsea','gsva'))))
}

get_auc_boxplot <- function(df, .type, ylabel='AUROC'){
  .type <- enquo(.type)

  order <- df %>%
    group_by(set_name) %>%
    summarize(median=median(!!.type), .groups='drop') %>%
    arrange(desc(median)) %>%
    distinct(set_name) %>%
    pull(set_name)
  df$set_name <- factor(df$set_name, levels = rev(order))

  df
}

# Load data ------------------------------------------------
# Load final networks
input_path_dorothea <- c(file.path('data', "dorothea", "dorothea_A_new.rds"),
                         file.path('data', "dorothea", "dorothea_ABC_new.rds"))
input_path_v1 <- file.path('data', 'networks_v1', 'network_collection_v1.rds')
input_path_v2 <- file.path('data', 'networks_v2', 'network_collection_v2.rds')
input_path_v2_dbTF <- file.path('data', 'networks_dbTF_v2', 'network_collection_v2_dbTF.rds')

network_collection_v1 <- readRDS(input_path_v1)$ExTRI_comp_unrestricted
network_collection_v2 <- readRDS(input_path_v2)$ExTRI_comp_unrestricted
network_collection_v2_dbTF <- readRDS(input_path_v2_dbTF)$ExTRI_comp_unrestricted

networks <- rbind(network_collection_v1, network_collection_v2, network_collection_v2_dbTF)
path_networks <- c(input_path_dorothea, networks$path)

final_networks <- map(path_networks, readRDS)
names(final_networks) <- c("Dorothea A", "Dorothea ABC",
                           "NTNU.1", "NTNU.1 s", "NTNU.1 w", "NTNU.1 s+w",
                           "NTNU.2", "NTNU.2 s", "NTNU.2 w", "NTNU.2 s+w",
                           "NTNU.2 dbTF", "NTNU.2 dbTF s", "NTNU.2 dbTF w", "NTNU.2 dbTF s+w")

# Load benchmark
bmeta_knockTF <- readRDS(file.path('data',"bench", "knockTF_meta.rds"))
bexpr_knockTF <-  readRDS(file.path('data',"bench", "knockTF_exp.rds"))
bmeta_rna <- readRDS(file.path('data',"bench", "rna_meta.rds"))
bexpr_rna <-  readRDS(file.path('data',"bench", "rna_expr.rds"))
perturbedTFs <- unique(c(bmeta_knockTF$target, bmeta_rna$target))

# Load TF activities
decoupler_path <- list.files("output/new", pattern = "decoupler", full.names = T)
decoupler_path_100 <- decoupler_path[!(str_detect(decoupler_path, "2000") | str_detect(decoupler_path, "1000"))]
decoupler_path_1000 <- decoupler_path[str_detect(decoupler_path, "1000")]
decoupler_path_2000 <- decoupler_path[str_detect(decoupler_path, "2000")]

decoupler_res_100 <- map(decoupler_path_100, readRDS)
names(decoupler_res_100) <- names(final_networks)

decoupler_res_1000 <- map(decoupler_path_1000, readRDS)
names(decoupler_res_1000) <- names(final_networks)

decoupler_res_2000 <- map(decoupler_path_2000, readRDS)
names(decoupler_res_2000) <- names(final_networks)[c(1,2,11,12,13,14)]


# Compare unrestricted and restricted network performance ------------------------------------------------
rna_result <- readRDS(file.path('output', 'new', 'estimate_rna_signed_test_networks.rds')) %>% filter(statistic == "consensus")
knockTF_result <- readRDS(file.path('output', 'new', 'estimate_knockTF_signed_test_networks.rds')) %>% filter(statistic == "consensus")
rna_result_unres <- readRDS(file.path('output', 'new', 'estimate_rna_signed_test_networks_unres.rds')) %>% filter(statistic == "consensus")
knockTF_result_unres <- readRDS(file.path('output', 'new', 'estimate_knockTF_signed_test_networks_unres.rds')) %>% filter(statistic == "consensus")

rna_result$set_name <- c("Dorothea ABC", "Dorothea A",
                         "NTNU.2 dbTF w", "NTNU.2 dbTF", "NTNU.2 dbTF w+s", "NTNU.2 dbTF s",
                         "NTNU.1 w", "NTNU.1", "NTNU.1 w+s", "NTNU.1 s",
                         "NTNU.2 w", "NTNU.2", "NTNU.2 w+s", "NTNU.2 s")
rna_result_unres$set_name <- rna_result$set_name

knockTF_result$set_name <- c("Dorothea ABC", "Dorothea A",
                             "NTNU.2 dbTF w", "NTNU.2 dbTF", "NTNU.2 dbTF w+s", "NTNU.2 dbTF s",
                             "NTNU.1 w", "NTNU.1", "NTNU.1 w+s", "NTNU.1 s",
                             "NTNU.2 w", "NTNU.2", "NTNU.2 w+s", "NTNU.2 s")
knockTF_result_unres$set_name <- knockTF_result$set_name


# Generate data-frames
rna_roc_df <- get_auc_df(rna_result, roc)
rna_prc_df <- get_auc_df(rna_result, prc)
rna_roc_unres_df <- get_auc_df(rna_result_unres, roc)
rna_prc_unres_df <- get_auc_df(rna_result_unres, prc)
knockTF_roc_df <- get_auc_df(knockTF_result, roc)
knockTF_prc_df <- get_auc_df(knockTF_result, prc)
knockTF_roc_unres_df <- get_auc_df(knockTF_result_unres, roc)
knockTF_prc_unres_df <- get_auc_df(knockTF_result_unres, prc)

# Generate plots
rna_roc_boxp <- get_auc_boxplot(rna_roc_df, roc, 'AUROC')
rna_prc_boxp <- get_auc_boxplot(rna_prc_df, prc, 'AUPRC')
rna_roc_unres_boxp <- get_auc_boxplot(rna_roc_unres_df, roc, 'AUROC')
rna_prc_unres_boxp <- get_auc_boxplot(rna_prc_unres_df, prc, 'AUPRC')
knockTF_roc_boxp <- get_auc_boxplot(knockTF_roc_df, roc, 'AUROC')
knockTF_prc_boxp <- get_auc_boxplot(knockTF_prc_df, prc, 'AUPRC')
knockTF_roc_unres_boxp <- get_auc_boxplot(knockTF_roc_unres_df, roc, 'AUROC')
knockTF_prc_unres_boxp <- get_auc_boxplot(knockTF_prc_unres_df, prc, 'AUPRC')

rna_roc_wilc <- map_dfr(as.character(unique(rna_roc_boxp$set_name)), function(net_name){
  wilc_results <- wilcox.test(rna_roc_unres_boxp %>% filter(set_name == net_name) %>% pull(roc),
                              rna_roc_boxp %>% filter(set_name == net_name) %>% pull(roc),
                              alternative = "greater")
  c(network = net_name, wilcoxon.p = wilc_results$p.value)
})
rna_roc_wilc$wilcoxon.p <- p.adjust(rna_roc_wilc$wilcoxon.p)

rna_prc_wilc <- map_dfr(as.character(unique(rna_prc_boxp$set_name)), function(net_name){
  wilc_results <- wilcox.test(rna_prc_unres_boxp %>% filter(set_name == net_name) %>% pull(prc),
                              rna_prc_boxp %>% filter(set_name == net_name) %>% pull(prc),
                              alternative = "greater")
  c(network = net_name, wilcoxon.p = wilc_results$p.value)
})
rna_prc_wilc$wilcoxon.p <- p.adjust(rna_prc_wilc$wilcoxon.p)

knockTF_roc_wilc <- map_dfr(as.character(unique(knockTF_roc_boxp$set_name)), function(net_name){
  wilc_results <- wilcox.test(knockTF_roc_unres_boxp %>% filter(set_name == net_name) %>% pull(roc),
                              knockTF_roc_boxp %>% filter(set_name == net_name) %>% pull(roc),
                              alternative = "greater")
  c(network = net_name, wilcoxon.p = wilc_results$p.value)
})
knockTF_roc_wilc$wilcoxon.p <- p.adjust(knockTF_roc_wilc$wilcoxon.p)

knockTF_prc_wilc <- map_dfr(as.character(unique(knockTF_prc_boxp$set_name)), function(net_name){
  wilc_results <- wilcox.test(knockTF_prc_unres_boxp %>% filter(set_name == net_name) %>% pull(prc),
                              knockTF_prc_boxp %>% filter(set_name == net_name) %>% pull(prc),
                              alternative = "greater")
  c(network = net_name, wilcoxon.p = wilc_results$p.value)
})
knockTF_prc_wilc$wilcoxon.p <- p.adjust(knockTF_prc_wilc$wilcoxon.p)

wilc_results <- list(rna_roc_wilc, rna_prc_wilc, knockTF_roc_wilc, knockTF_prc_wilc) %>% reduce(full_join, by = "network")
colnames(wilc_results) <- c("network",  "dorothea_roc.p.adj","dorothea_prc.p.adj", "knockTF_roc.p.adj","knockTF_prc.p.adj")





# Network size ------------------------------------------------
# general size
map_df(final_networks, function(x){
  x_rna <- decoupleR::intersect_regulons(mat = bexpr_rna, network = x, .source = "source", .target = "target",minsize = 5)
  x_knockTF <- decoupleR::intersect_regulons(mat = bexpr_knockTF, network = x, .source = "source", .target = "target",minsize = 5)
  c(nTFs = x$source %>% unique() %>% length(),
    nTFs5 = sum(table(x$source) >= 5),
    nTFsBench5 =  x_rna$source %>% unique() %>% length(),
    nTFsKnock5 =  x_knockTF$source %>% unique() %>% length(),
    nTFsBench =  sum(unique(x_rna$source) %in% unique(bmeta_rna$target)),
    nTFsKnock =  sum(unique(x_knockTF$source) %in% unique(bmeta_knockTF$target)))
})

final_networks <- map(final_networks, function(x){
  x %>% filter(source %in% names(table(x$source))[table(x$source) >= 5])
})


# overlap
# Dataset + TF in benchmark (dbTF, coTF)
overlap_networks <- final_networks[c("Dorothea A", "Dorothea ABC", "NTNU.1", "NTNU.2", "NTNU.2 dbTF")]
TFs <- map(overlap_networks, function(x){table(x$source) %>% as.data.frame() %>% filter(Freq >= 5) %>% mutate(Freq = 1)})

TF_df <- TFs %>% reduce(full_join, by = "Var1")
TF_df[is.na(TF_df)] <- 0

colnames(TF_df) <- c("source", names(overlap_networks))

TF_upset_p <- upset(TF_df, sets = c(names(overlap_networks)), order.by = "freq", text.scale = 1.5) %>% as.ggplot() +
  ggtitle("TFs") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14)
  )


# edges
edges <- map(overlap_networks, function(x){
  tf_5targets <- table(x$source) %>% as.data.frame() %>% filter(Freq >= 5) %>% pull(Var1)
  x <- x %>% filter(source %in% tf_5targets)
  data.frame(edge = unique(paste0(x$source, "_", x$target)),
             count = 1)})

edge_df <- edges %>% reduce(full_join, by = "edge")
edge_df[is.na(edge_df)] <- 0

colnames(edge_df) <- c("edge", names(overlap_networks))

edges_upset_p <- upset(edge_df, sets = c(names(overlap_networks)), order.by = "freq",text.scale = 1.5) %>% as.ggplot() +
  ggtitle("Edges") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 14))

pdf(file = file.path('figures', 'final_comp', 'new', 'overlap_networks.pdf'),
    width = 10, # The width of the plot in inches
    height = 5)
(TF_upset_p | edges_upset_p)+
  plot_layout(guides = 'collect', widths = c(1, 1))
dev.off()



# Mode of regulation ------------------------------------------------
mor_distribution_p <- map(1:length(final_networks), function(i){
  network <- final_networks[[i]]
  if (str_detect(names(final_networks)[[i]], "w") | str_detect(names(final_networks)[[i]], "Dorothea")) {
    ggplot(network, aes(x=mor)) +
      geom_histogram(aes(y=..density..),
                     binwidth = 0.005,
                     colour="black", fill="white") +
      stat_function(
        fun = dnorm,
        args = list(mean = mean(network$mor), sd = sd(network$mor)),
        col = 'black'
      ) + ggtitle(names(final_networks)[i]) +
      theme_bw() + theme(text = element_text(size = 14))
  } else {
    NULL
  }


})

mor_distribution_p_new <- mor_distribution_p[!map_lgl(mor_distribution_p, is.null)]
pdf(file = "figures/final_comp/new/mor_distribution.pdf",
    width = 10,
    height = 15)

cowplot::plot_grid(plotlist = mor_distribution_p_new, ncol = 2)

dev.off()


# Regulon composition ------------------------------------------------
networks_TF <- final_networks[c("Dorothea A", "Dorothea ABC", "NTNU.1 s", "NTNU.2 s", "NTNU.2 dbTF s")]

dist_df <- map_df(networks_TF, function(net){
  TF_split <- net %>% group_by(source) %>% group_split()
  dist <- map_df(TF_split, function(x){
    pos <- x %>% filter(mor > 0) %>% nrow()
    neg <- x %>% filter(mor < 0) %>% nrow()

    c(dist = neg/(neg+pos),  neg = (pos == 0), pos = (neg == 0),mixed = ((pos != 0 & neg != 0)))

  })

  data.frame(counts = colSums(dist)[2:4],
             `regulon type` = c("all negative", "all positiv", "mixed"))

})

dist_df$network <- rep(names(networks_TF), each = 3)
dist_df$`regulon type` <- factor(dist_df$regulon.type, levels = unique(dist_df$regulon.type))

p <- ggplot(data=dist_df, aes(x=network, y=counts, fill=`regulon type`)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_bw() + scale_fill_grey()  + theme(text = element_text(size = 14))

pdf("figures/final_comp/new/regulon_type.pdf", width = 10, height = 10)
plot(p) + ylab("number of TFs")
dev.off()

mean_mor_df <- map(networks_TF, function(net){
  TF_split <- net %>% group_by(source) %>% group_split()
  dist <- map_df(TF_split, function(x){
    c(mean = mean(abs(x$mor)), ntargets = nrow(x), max = max(abs(x$mor)), min = min(abs(x$mor)))
  })

})


dist_plot <- map(names(networks_TF), function(name){
  net <- networks_TF[[name]]
  TF_split <- net %>% group_by(source) %>% group_split()
  dist <- map_df(TF_split, function(x){
    pos <- x %>% filter(mor > 0) %>% nrow()
    neg <- x %>% filter(mor < 0) %>% nrow()

    c(dist = (pos/(neg+pos))*100,  neg = (pos == 0), pos = (neg == 0),mixed = ((pos != 0 & neg != 0)))

  })

  ggplot(dist, aes(x=dist)) + geom_density() + theme_bw() + ggtitle(name) + xlab("% positive edges") + theme(text = element_text(size = 14))

})


pdf("figures/final_comp/new/regulon_distribution.pdf", width = 10, height = 10)
cowplot::plot_grid(plotlist = dist_plot, ncol = 2 )
dev.off()

decoupler_mixed <- map_dfr(unique(decoupler_consensus$network), function(net_name){
  net <- final_networks[[net_name]]
  TF_split <- net %>% group_by(source) %>% group_split()
  dist <- map_df(TF_split, function(x){
    pos <- x %>% filter(mor > 0) %>% nrow()
    neg <- x %>% filter(mor < 0) %>% nrow()

    c(dist = neg/(neg+pos),  neg = (pos == 0), pos = (neg == 0),mixed = ((pos != 0 & neg != 0)), TF = unique(x$source))

  })
  decoupler_consensus %>% filter(network == net_name) %>% filter(!(source %in% dist$TF[as.logical(dist$mixed)]))
})

bmeta_id <- c(paste(bmeta_knockTF$id, bmeta_knockTF$target, sep = "_"))

bmeta_id_rna <- paste(bmeta_rna$id, bmeta_rna$target, sep = "_")
decoupler_mixed_filtered <- decoupler_mixed %>%
  mutate(bench_id = paste(decoupler_mixed$condition,
                          decoupler_mixed$source,
                          sep = "_")) %>%
  filter(bench_id %in% bmeta_id_rna)
decoupler_mixed_filtered$network <- factor(decoupler_mixed_filtered$network,  levels = unique(decoupler_mixed_filtered$network))

p <- ggplot(decoupler_mixed_filtered, aes(x=network, y=abs(score))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf("figures/final_comp/new/test_signeffect_dorothea.pdf", width = 10, height = 7)
plot(p)
dev.off()

