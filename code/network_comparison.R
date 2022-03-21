library(tidyverse)
library(UpSetR)
library(ggplotify)
library(patchwork)

# final networks
path_networks <- c(dorothea_A = "data/dorothea/dorothea_A_new.rds",
                   dorothea_ABC = "data/dorothea/dorothea_ABC_new.rds",
                   v1_weighted = "data/networks_v1/ExTRI_comp_scaled_FALSE_1_evidence_FALSE.rds",
                   v1_weighted_signed = "data/networks_v1/ExTRI_comp_scaled_TRUE_1_evidence_FALSE.rds",
                   v2_weighted = "data/networks_v2/ExTRI_comp_scaled_FALSE_0.8_evidence_FALSE.rds",
                   v2_weighted_signed = "data/networks_v2/ExTRI_comp_scaled_TRUE_0.8_evidence_FALSE.rds",
                   v2_dbTF_weighted = "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_0.8_evidence_FALSE.rds",
                   v2_dbTF_weighted_signed = "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_0.8_evidence_FALSE.rds")

final_networks <- map(path_networks, readRDS)
bmeta_knockTF <- readRDS(file.path('data',"bench", "knockTF_meta.rds"))
bexpr_knockTF <-  readRDS(file.path('data',"bench", "knockTF_exp.rds"))
bmeta_rna <- readRDS(file.path('data',"bench", "rna_meta.rds"))
bexpr_rna <-  readRDS(file.path('data',"bench", "rna_expr.rds"))


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

pdf(file = file.path('figures', 'final_comp', 'overlap_networks.pdf'),
    width = 10, # The width of the plot in inches
    height = 5)
(TF_upset_p | edges_upset_p)+
  plot_layout(guides = 'collect', widths = c(1, 1))
dev.off()


# plot mor distributaions
names(final_networks) <- c("Dorothea A", "Dorothea ABC", "NTNU.1 w", "NTNU.1 w+s", "NTNU.2 w", "NTNU.2 w+s", "NTNU.2 dbTF w", "NTNU.2 dbTF w+s")
mor_distribution_p <- map(1:length(final_networks), function(i){
  network <- final_networks[[i]]

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

})


# Save figures
pdf(file = "figures/final_comp/TF_overlap_knockTF.pdf",
    width = 10,
    height = 7)

TF_upset_p

dev.off()

pdf(file = "figures/final_comp/edge_overlap_knockTF.pdf",
    width = 10,
    height = 7)

edges_upset_p

dev.off()

pdf(file = "figures/final_comp/mor_distribution_knockTF.pdf",
    width = 10,
    height = 15)

cowplot::plot_grid(plotlist = mor_distribution_p, ncol = 2)

dev.off()


# Porpotion of negative signs per TF
# How many TFs have Mixed signs, only pos, only neg

networks_TF <- final_networks[c("Dorothea A", "Dorothea ABC", "NTNU.1 w+s", "NTNU.2 w+s", "NTNU.2 dbTF w+s")]

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

pdf("figures/regulon_type.pdf", width = 10, height = 10)
plot(p) + ylab("number of TFs")
dev.off()

mean_mor_df <- map(networks_TF, function(net){
  TF_split <- net %>% group_by(source) %>% group_split()
  dist <- map_df(TF_split, function(x){
    c(mean = mean(abs(x$mor)), ntargets = nrow(x), max = max(abs(x$mor)), min = min(abs(x$mor)))
  })

})

map(mean_mor_df, function(x){
  ggplot(x, aes(ntargets, mean)) + geom_point() + theme_bw()
  cor(x$ntargets, x$mean)
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


pdf("figures/regulon_distribution.pdf", width = 10, height = 10)
cowplot::plot_grid(plotlist = dist_plot, ncol = 2 )
dev.off()

