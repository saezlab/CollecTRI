library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggrepel)
library(patchwork)
library(ggplotify)

# Create dir
path_figs <- file.path('figures', 'final_comp', 'knockTF')
dir.create(path_figs, showWarnings = F, recursive = T)

# Plot functions
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

  ggplot(df, aes(x=set_name, y=!!.type)) +
    geom_boxplot() +
    xlab('Network') +
    ylab(ylabel) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

get_auc_scatter <- function(df){
  min_lim <- floor(min(c(df$roc, df$prc)) * 100)/100
  max_lim <- ceiling(max(c(df$roc, df$prc)) * 100)/100
  ggplot(df, aes(x=roc, y=prc, label=set_name)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, max.time=5, max.iter=1000000) +
    theme(text = element_text(size=14)) +
    xlab('AUROC') +
    ylab('AUPRC') +
    xlim(min_lim,max_lim) +
    ylim(min_lim,max_lim) +
    theme_bw()
}

# Read
rna_result <- readRDS(file.path('output', 'estimate_dorothea_sign.rds'))@bench_res
knockTF_result <- readRDS(file.path('output','estimate_knockTF_sign_new.rds'))@bench_res

rna_result <- rna_result %>% filter(statistic == "consensus")
rna_result$set_name <- c("Dorothea A", "Dorothea ABC",
                         "NTNU.2 dbTF w", "NTNU.2 dbTF", "NTNU.2 dbTF w+s", "NTNU.2 dbTF s",
                         "NTNU.1 w", "NTNU.1", "NTNU.1 w+s", "NTNU.1 s",
                         "NTNU.2 w", "NTNU.2", "NTNU.2 w+s", "NTNU.2 s")

knockTF_result <- knockTF_result %>% filter(statistic == "consensus")
knockTF_result$set_name <- c("Dorothea A", "Dorothea ABC",
                             "NTNU.2 dbTF w", "NTNU.2 dbTF", "NTNU.2 dbTF w+s", "NTNU.2 dbTF s",
                             "NTNU.1 w", "NTNU.1", "NTNU.1 w+s", "NTNU.1 s",
                             "NTNU.2 w", "NTNU.2", "NTNU.2 w+s", "NTNU.2 s")

#knockTF_result$set_name <- c("Dorothea A", "Dorothea ABC",
#                             "NTNU.2 0", "NTNU.2 0.1", "NTNU.2 0.2", "NTNU.2 0.3",
#                             "NTNU.2 0.4", "NTNU.2 0.5", "NTNU.2 0.6", "NTNU.2 0.7",
#                             "NTNU.2 0.8", "NTNU.2 0.9", "NTNU.2 1")

# Generate data-frames
rna_roc_df <- get_auc_df(rna_result, roc)
rna_prc_df <- get_auc_df(rna_result, prc)
knockTF_roc_df <- get_auc_df(knockTF_result, roc)
knockTF_prc_df <- get_auc_df(knockTF_result, prc)

rna_auc_df <- rna_roc_df %>%
  left_join(rna_prc_df) %>%
  group_by(set_name) %>%
  summarise(roc = median(roc), prc = median(prc), .groups='drop')
knockTF_auc_df <- knockTF_roc_df %>%
  left_join(knockTF_prc_df) %>%
  group_by(set_name) %>%
  summarise(roc = median(roc), prc = median(prc), .groups='drop')

# Generate plots
rna_roc_boxp <- get_auc_boxplot(rna_roc_df, roc, 'AUROC')
rna_prc_boxp <- get_auc_boxplot(rna_prc_df, prc, 'AUPRC')
knockTF_roc_boxp <- get_auc_boxplot(knockTF_roc_df, roc, 'AUROC')
knockTF_prc_boxp <- get_auc_boxplot(knockTF_prc_df, prc, 'AUPRC')

rna_auc_scatt <- get_auc_scatter(rna_auc_df)
knockTF_auc_scatt <- get_auc_scatter(knockTF_auc_df)

# Merge together and save
pdf(file = file.path(path_figs, 'dorothea_bench_sign.pdf'),
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
((rna_roc_boxp / rna_prc_boxp) | rna_auc_scatt) +
  plot_layout(guides = 'collect', widths = c(1, 2))  +
  plot_annotation(tag_levels = 'A', title = "Dorothea bench")
dev.off()

pdf(file = file.path(path_figs, 'knockTF_bench_sign_new.pdf'),
    width = 7.5, # The width of the plot in inches
    height = 7) # The height of the plot in inches
((knockTF_roc_boxp / knockTF_prc_boxp) | knockTF_auc_scatt)  +
  plot_layout(guides = 'collect', widths = c(1, 2))  +
  plot_annotation(tag_levels = 'A', title = "knockTF bench")
dev.off()

# Test bias of methods
path_networks <- c(dorothea_A = "data/dorothea/dorothea_A_new.rds",
                   dorothea_ABC = "data/dorothea/dorothea_ABC_new.rds",
                   v1_weighted = "data/networks_v1/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v1_weighted_signed = "data/networks_v1/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v1 = "data/networks_v1/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                   v1_signed = "data/networks_v1/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                   v2_weighted = "data/networks_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v2_weighted_signed = "data/networks_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v2 = "data/networks_v2/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                   v2_signed = "data/networks_v2/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                   v2_dbTF_weighted = "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v2_dbTF_weighted_signed = "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v2_dbTF = "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                   v2_dbTF_signed = "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds")

final_networks <- map(path_networks, readRDS)
names(final_networks) <- names(decoupler_act) <- c("Dorothea A", "Dorothea ABC",
                                                   "NTNU.1 w", "NTNU.1 w+s", "NTNU.1", "NTNU.1 s",
                                                   "NTNU.2 w", "NTNU.2 w+s", "NTNU.2", "NTNU.2 s",
                                                   "NTNU.2 dbTF w", "NTNU.2 dbTF w+s", "NTNU.2 dbTF", "NTNU.2 dbTF s")

decoupler_act <- readRDS("output/decoupler_act.rds")
decoupler_act <- decoupler_act[c(5,6,3,4,9,10,7,8,13,14,11,12,1,2)]

corr_plot_list <- map(names(decoupler_act), function(name_act){
  act_net <- decoupler_act[[name_act]]
  ggplot(act_net, aes(abs(score), Freq)) +
    geom_point() +
    ggtitle(paste0(name_act, " Pearson correlation = ", round(cor(abs(act_net$score), act_net$Freq), digits = 3))) +
    theme(text = element_text(size=14)) +
    xlab('Absolute TF activity score') +
    ylab('Regulon size') +
    theme_bw()
})

pdf(file = file.path(path_figs, 'corr_regulon_size.pdf'),
        width = 7.5, # The width of the plot in inches
        height = 9)
cowplot::plot_grid(plotlist = corr_plot_list, ncol = 4)
dev.off()

weighted_results <- decoupler_act[c(3,7,11,14)]
mor_corr_plot_list <- map(names(weighted_results), function(net_name){
  act_res <- weighted_results[[net_name]]
  net <- final_networks[[net_name]]
  mean_mor <- net %>% group_by(source) %>% summarise(mean_mor = mean(abs(mor)))

  act_res_mor <- left_join(act_res, mean_mor)

  ggplot(act_res_mor, aes(abs(score), mean_mor)) +
    geom_point() +
    ggtitle(paste0(net_name, " Pearson correlation = ", round(cor(abs(act_res_mor$score), act_res_mor$mean_mor), digits = 3))) +
    theme(text = element_text(size=14)) +
    xlab('Absolute TF activity score') +
    ylab('Mean mode of regulation') +
    theme_bw()
})

pdf(file = file.path(path_figs, 'corr_mean_mor.pdf'),
    width = 7.5, # The width of the plot in inches
    height = 9)
cowplot::plot_grid(plotlist = mor_corr_plot_list, ncol = 2)
dev.off()

