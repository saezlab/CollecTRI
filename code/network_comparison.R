library(tidyverse)
library(UpSetR)

# final networks
path_networks <- c(dorothea_A = "data/dorothea/dorothea_A.rds",
                   dorothea_ABC = "data/dorothea/dorothea_ABC.rds",
                   v1_weighted = "data/networks_v1/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v1_weighted_signed = "data/networks_v1/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v2_weighted = "data/networks_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v2_weighted_signed = "data/networks_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                   v2_dbTF_weighted = "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                   v2_dbTF_weighted_signed = "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds")

final_networks <- map(path_networks, readRDS)
bmeta_knockTF <- readRDS(file.path('data',"bench", "knockTF_meta.rds"))
bexpr_knockTF <-  readRDS(file.path('data',"bench", "knockTF_expr.rds"))
bmeta_rna <- readRDS(file.path('data',"bench", "rna_meta.rds"))
bexpr_rna <-  readRDS(file.path('data',"bench", "rna_expr.rds"))
bench_TF_df <- data.frame(Var1 = unique(bmeta$target),
                          x = 1)

# Size of networks
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


# plot overlapp
# Dataset + TF in benchmark (dbTF, coTF)
overlap_networks <- final_networks[c("dorothea_A", "dorothea_ABC", "v1_weighted_signed", "v2_weighted_signed", "v2_dbTF_weighted_signed")]
names(overlap_networks) <- c("Dorothea A", "Dorothea ABC", "NTNU.1", "NTNU.2", "NTNU.2 dbTF")
TFs <- map(overlap_networks, function(x){table(x$source) %>% as.data.frame() %>% filter(Freq >= 5) %>% mutate(Freq = 1)})

TF_df <- TFs %>% reduce(full_join, by = "Var1")
TF_df[is.na(TF_df)] <- 0

colnames(TF_df) <- c("source", names(overlap_networks))

TF_upset_p <- upset(TF_df, sets = c(names(overlap_networks)), order.by = "freq", text.scale = 1.5) %>% as.ggplot() + ggtitle("TFs") + theme(
  plot.title = element_text(hjust = 0.5)
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

edges_upset_p <- upset(edge_df, sets = c(names(overlap_networks)), order.by = "freq",text.scale = 1.5) %>% as.ggplot() + ggtitle("Edges") + theme(
  plot.title = element_text(hjust = 0.5)
)

pdf(file = file.path('figures', 'final_comp', 'overlap_networks.pdf'),
    width = 10, # The width of the plot in inches
    height = 5)
(TF_upset_p | edges_upset_p)+
  plot_layout(guides = 'collect', widths = c(1, 1))
dev.off()


# plot mor distributaions
mor_distribution_p <- map(3:length(final_networks), function(i){
  network <- final_networks[[i]]

  ggplot(network, aes(x=mor)) +
    geom_histogram(aes(y=..density..),
                   binwidth = 0.005,
                   colour="black", fill="white") +
    stat_function(
      fun = dnorm,
      args = list(mean = mean(network$mor), sd = sd(network$mor)),
      col = 'red'
    ) + ggtitle(names(path_networks)[i])

})


# Save figures
pdf(file = "figures/final_comp/TF_overlap_knockTF.pdf",
    width = 10,
    height = 7)

TF_upset_p

dev.off()

png(file = "figures/final_comp/TF_overlap_knockTF.png",
    width = 10,
    height = 7, units = "in", res = 800)

TF_upset_p

dev.off()

pdf(file = "figures/final_comp/edge_overlap_knockTF.pdf",
    width = 10,
    height = 7)

edges_upset_p

dev.off()

png(file = "figures/final_comp/edge_overlap_knockTF.png",
    width = 10,
    height = 7, units = "in", res = 800)

edges_upset_p

dev.off()



pdf(file = "figures/final_comp/mor_distribution_knockTF.pdf",
    width = 10,
    height = 15)

cowplot::plot_grid(plotlist = mor_distribution_p, ncol = 2)

dev.off()

png(file = "figures/final_comp/mor_distribution_knockTF.png",
    width = 10,
    height = 15, units = "in", res = 800)

cowplot::plot_grid(plotlist = mor_distribution_p, ncol = 2)

dev.off()

