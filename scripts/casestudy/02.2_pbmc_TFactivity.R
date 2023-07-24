## We load the required packages
library(Seurat)
library(cluster)
library(decoupleR)
library(tidyverse)
library(ggsignif)
library(ggplot2)
library(ggpubr)


# Load processed pbmc data from decoupler
data <- readRDS("output/case_study/single_cell/pbmc3k_processed.rds")

collecTRI <- read_csv("output/CollecTRI/CollecTRI_GRN.csv")
collecTRI$mor <- collecTRI$weight
mat <- as.matrix(data@assays$RNA@data)


# TF activity estimation ---------------------------------------------------------------------------------
if(!file.exists("output/case_study/single_cell/pbmc3k_activity.rds")){
  activity_scores <- run_ulm(mat = mat, network = collecTRI) # takes a few hours
  system("say activity estimation finished!")
  activity_scores_wide <- activity_scores %>%
    dplyr::filter(statistic == "ulm") %>%
    dplyr::select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()

  activity_scores_assay <- CreateAssayObject(counts = activity_scores_wide)

  data[["TFact"]] <- activity_scores_assay


  # TF expression ---------------------------------------------------------------------------------
  mat_TF <- as.matrix(data@assays$RNA@data)[rownames(as.matrix(data@assays$RNA@data)) %in% rownames(activity_scores_wide),]
  TF_expr <- CreateAssayObject(counts = mat_TF)
  data[["TFexpr"]] <- TF_expr

  saveRDS(data, "output/case_study/single_cell/pbmc3k_activity.rds")
} else {
  data <- readRDS("output/case_study/single_cell/pbmc3k_activity.rds")
}

# UMAP ---------------------------------------------------------------------------------
data <- Seurat::ScaleData(data, assay = "TFact")
data <- Seurat::ScaleData(data, assay = "TFexpr")

data <- RunPCA(data, assay = "TFact", verbose = FALSE, reduction.name = "act.pca", reduction.key = "actPCA_", features = GetAssayData(data, "TFact", slot = "data") %>% rownames())
data <- RunPCA(data, assay = "TFexpr", verbose = FALSE, reduction.name = "expr.pca", reduction.key = "exprPCA_", features = GetAssayData(data, "TFexpr", slot = "data") %>% rownames())

data <- RunUMAP(data, assay = "TFact", reduction = 'act.pca', dims = 1:30, verbose = FALSE, reduction.name = "act.umap", reduction.key = "actUMAP_")
data <- RunUMAP(data, assay = "TFexpr", reduction = 'expr.pca', dims = 1:30, verbose = FALSE, reduction.name = "expr.umap", reduction.key = "exprUMAP_")


# Activity versus expression ---------------------------------------------------------------------------------
old_idents <- Idents(data)
new_idents <- recode_factor(old_idents,
                            "Memory CD4 T" = "Memory CD4+ T cells",
                            "B" = "B cells",
                            "Naive CD4 T" = "Naive CD4+ T cells",
                            "CD14+ Mono" = "CD14+ Monocytes",
                            "CD8 T" = "CD8+ T cells",
                            "FCGR3A+ Mono" = "FCGR3A+ Monocytes",
                            "NK" = "Natural Killer Cells",
                            "DC" = "Dendritic cells")

Idents(data) <- new_idents
DefaultAssay(data = data) <- "TFact"
markersAct <- FindAllMarkers(data,
                             logfc.threshold = 0.5,
                             test.use = "wilcox",
                             min.pct = 0.15,
                             only.pos = T)

DefaultAssay(data = data) <- "TFexpr"
markersExpr <- FindAllMarkers(data,
                              logfc.threshold = 0.5,
                              test.use = "wilcox",
                              min.pct = 0.15,
                              only.pos = T)

celltypes <- unique(Idents(data))

write_csv(markersExpr, "output/case_study/single_cell/SuppFile3_TF_expression_markers.csv")
write_csv(markersAct, "output/case_study/single_cell/SuppFile4_TF_activity_markers.csv")

markEx <- markersExpr %>% mutate(id = paste(gene, cluster, sep = "_")) %>% pull(id)
markAc <- markersAct %>% mutate(id = paste(gene, cluster, sep = "_")) %>% pull(id)

markEx[!markEx %in% markAc]
diffTFs <- map_dfr(celltypes, function(celltype){
  act_gene <- markersAct %>%
    filter(cluster == celltype) %>%
    filter(p_val_adj <= 0.05) %>%
    pull(gene)

  expr_gene <- markersExpr %>%
    filter(cluster == celltype) %>%
    filter(p_val_adj <= 0.05) %>%
    pull(gene)

  data.frame(celltype = celltype,
             n = c(intersect(act_gene, expr_gene) %>% length(),
                   sum(!act_gene %in% expr_gene),
                   sum(!expr_gene %in% act_gene)),
             type = c("shared", "activity", "expression"))
})

order_ct <- diffTFs %>%
  group_by(celltype) %>%
  summarise(total = sum(n)) %>%
  arrange(total) %>%
  pull(celltype)

diffTFs$celltype <- factor(diffTFs$celltype, levels = order_ct)
pDiffTFs <- ggplot(diffTFs, aes(fill=type, y=n, x=celltype)) +
  geom_bar(position="stack", stat="identity")  +
  scale_fill_manual(values=c("#99b9a3", "#E8D4B2", "#b2c6e8")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines")) +
  ylab("Total number of differentially\nexpressed or active TFs") +
  xlab("Cell type")

diffTFs %>% group_by(type) %>% summarize(total = sum(n)) %>% pull(total) %>% sum()
diffTFs %>% group_by(type) %>% summarize(total = sum(n))
# Features B cells, NK cells ---------------------------------------------------------------------------------
DefaultAssay(data = data) <- "TFact"
pPAX5act <- (FeaturePlot(data, features = c("PAX5"), reduction = "umap",pt.size = 0.01) &
               scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red')) +
  ggtitle('PAX5 activity') +
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', linewidth = 0.3),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("UMAP 1") +
  ylab ("UMAP 2")

DefaultAssay(object = data) <- "RNA"
pPAX5exp <- FeaturePlot(data, features = c("PAX5"), pt.size = 0.01) + ggtitle('PAX5 expression')+
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', linewidth = 0.3),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("UMAP 1") +
  ylab ("UMAP 2")

DefaultAssay(data = data) <- "TFact"
pEOMESact <- (FeaturePlot(data, features = c("EOMES"), reduction = "umap", pt.size = 0.01) & # or PRDM1 https://www.frontiersin.org/articles/10.3389/fimmu.2020.01945/full#:~:text=Multiple%20transcriptional%20factors%20mediate%20regulation,critical%20for%20NK%20cell%20maturation.
                scale_colour_gradient2(low = 'blue', mid = 'lightgrey', high = 'red')) +
  ggtitle('EOMES activity')+
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', linewidth = 0.3),
        axis.ticks = element_line(colour = "black", linewidth = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("UMAP 1") +
  ylab ("UMAP 2")

DefaultAssay(object = data) <- "RNA"
pEOMESexp <- FeaturePlot(data, features = c("EOMES"), pt.size = 0.01) + ggtitle('EOMES expression')+
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("UMAP 1") +
  ylab ("UMAP 2")

# Violin Plots ---------------------------------------------------------------------------------
vln_df <- rbind(data.frame(PAX5 = c(data[["TFact"]]@data["PAX5", names(Idents(data))[Idents(data) == "B cells"]],
                                    data[["TFexpr"]]@data["PAX5",names(Idents(data))[Idents(data) == "B cells"]]),
                           TF = rep(c("TF activity", "TF expression"), each = sum(Idents(data) == "B cells")),
                           celltype = "B cells"),
                data.frame(PAX5 = c(data[["TFact"]]@data["PAX5", names(Idents(data))[!Idents(data) == "B cells"]],
                                    data[["TFexpr"]]@data["PAX5",names(Idents(data))[!Idents(data) == "B cells"]]),
                           TF = rep(c("TF activity", "TF expression"), each = sum(!Idents(data) == "B cells")),
                           celltype = "others"))

sum(vln_df %>% filter(TF == "TF expression" & celltype == "B cells") %>% pull(PAX5) != 0)
vln_df %>% filter(TF == "TF expression" & celltype == "B cells") %>% pull(PAX5) %>% length()

# Add noise as done for Seurat violin plot
PAX5_vln_df <- vln_df
noise <- rnorm(n = length(x = PAX5_vln_df[, "PAX5"])) / 100000
PAX5_vln_df$PAX5 <- PAX5_vln_df$PAX5  + noise

markersAct %>% filter(cluster == "B cells") %>% filter(gene == "PAX5")
markersExpr %>% filter(cluster == "B cells") %>% filter(gene == "PAX5")


PAX5_vln_df
PAX5_vln <- ggplot(PAX5_vln_df %>%
                     filter(TF == "TF activity"), aes(x = celltype, y = PAX5, fill = celltype)) +
  geom_violin(adjust = 1,trim=TRUE, scale = "width", linewidth = 0.2) +
  geom_jitter(size = 0, stroke = 0, shape = 16)  +
  theme_minimal() +
  geom_signif(y_position = c(5.5), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              margin_top = 0)+
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(c(-2, 8.1)) +
  xlab("") +
  ylab("Activity")

PAX5_vln_2 <- ggplot(PAX5_vln_df %>%
                     filter(TF == "TF expression"), aes(x = celltype, y = PAX5, fill = celltype)) +
  geom_violin(adjust = 1,trim=TRUE, scale = "width", linewidth = 0.2) +
  geom_jitter(size = 0.2, stroke = 0, shape = 16)  +
  theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("") +
  ylab("Expression")

# EOMES
vln_df <- rbind(data.frame(EOMES = c(data[["TFact"]]@data["EOMES", names(Idents(data))[Idents(data) == "Natural Killer Cells"]],
                                    data[["TFexpr"]]@data["EOMES",names(Idents(data))[Idents(data) == "Natural Killer Cells"]]),
                           TF = rep(c("TF activity", "TF expression"), each = sum(Idents(data) == "Natural Killer Cells")),
                           celltype = "NK Cells"),
                data.frame(EOMES = c(data[["TFact"]]@data["EOMES", names(Idents(data))[!Idents(data) == "Natural Killer Cells"]],
                                    data[["TFexpr"]]@data["EOMES",names(Idents(data))[!Idents(data) == "Natural Killer Cells"]]),
                           TF = rep(c("TF activity", "TF expression"), each = sum(!Idents(data) == "Natural Killer Cells")),
                           celltype = "others"))

sum(vln_df %>% filter(TF == "TF expression" & celltype == "NK Cells") %>% pull(EOMES) != 0)
vln_df %>% filter(TF == "TF expression" & celltype == "NK Cells") %>% pull(EOMES) %>% length()

# Add noise as done for Seurat violin plot
EOMES_vln_df <- vln_df
noise <- rnorm(n = length(x = EOMES_vln_df[, "EOMES"])) / 100000
EOMES_vln_df$EOMES <- EOMES_vln_df$EOMES  + noise

markersAct %>% filter(cluster == "Natural Killer Cells") %>% filter(gene == "EOMES")
markersExpr %>% filter(cluster == "Natural Killer Cells") %>% filter(gene == "EOMES")


EOMES_vln_df
EOMES_vln <- ggplot(EOMES_vln_df %>%
                     filter(TF == "TF activity"), aes(x = celltype, y = EOMES, fill = celltype)) +
  geom_violin(adjust = 1,trim=TRUE, scale = "width", linewidth = 0.2) +
  geom_jitter(size = 0, stroke = 0, shape = 16)  +
  theme_minimal() +
  geom_signif(y_position = c(7), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              margin_top = 0)+
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylim(c(-2, 8.5)) +
  xlab("") +
  ylab("Activity")

EOMES_vln_2 <- ggplot(EOMES_vln_df %>%
                       filter(TF == "TF expression"), aes(x = celltype, y = EOMES, fill = celltype)) +
  geom_violin(adjust = 1,trim=TRUE, scale = "width", linewidth = 0.2) +
  geom_jitter(size = 0.2, stroke = 0, shape = 16)  +
  theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("") +
  ylab("Expression")

# Calculate Silhouette score ---------------------------------------------------------------------------------
# silhouette metric
reductions <- c("expr.pca", "act.pca")

silhouette_df <- map_dfr(reductions, function(red){
  reduction <- red
  dims <- 1:30
  dist.matrix <- dist(x = Embeddings(object = data[[reduction]])[, dims])
  clusters <- Idents(data)
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)

  data.frame(sil.score = sil[, 3],
             reduction = red)

})
silhouette_df <- silhouette_df %>%
  mutate(reduction = recode(reduction,
                            act.pca = "TF activity",
                            expr.pca = "TF expression"))

t.test(silhouette_df$sil.score[silhouette_df$reduction == "TF activity"],
       silhouette_df$sil.score[silhouette_df$reduction == "TF expression"] )$p.value
t.test(silhouette_df$sil.score[silhouette_df$reduction == "TF activity"],
       silhouette_df$sil.score[silhouette_df$reduction == "TF expression"] )$statistic

p <- ggplot(silhouette_df, aes(x=reduction, y=sil.score, fill = reduction)) +
  geom_boxplot(outlier.size=0.2, lwd=0.2) +
  scale_colour_manual(values=c("#99b9a3", "#E8D4B2")) +
  scale_fill_manual(values=c("#99b9a3", "#E8D4B2")) +
  geom_signif(y_position = c(0.5), xmin = c(1), xmax = c(2),
              annotation = c("***"), tip_length = 0.02, size = 0.25, textsize = 2.7,
              margin_top = 0) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),) +
  ylab("Average silhouette width") +
  xlab("") +
  ylim(c(-0.3, 0.56))
p
# Save plots ---------------------------------------------------------------------------------
pdf("figures/manuscript/p4.1.pdf", width = 3, height = 2.8)
DimPlot(data, reduction = "umap", pt.size = 0.01) + theme_minimal() +
  theme(legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.position = "none",
        axis.line = element_line(colour = 'black', size = 0.3),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("UMAP 1") +
  ylab ("UMAP 2")
dev.off()

pdf("figures/manuscript/p4.1_legend.pdf", width = 8, height = 3)
DimPlot(data, reduction = "umap", pt.size = 0.01) +
  guides(color = guide_legend(override.aes = list(size=2), ncol=3) ) +
  theme(legend.key.size = unit(0.1, "cm"),
        legend.text=element_text(size=9),
        legend.position = "bottom")
dev.off()


pdf("figures/manuscript/p4.2.pdf", width = 2.7, height = 3.2)
pDiffTFs
dev.off()

pdf("figures/manuscript/p4.3.1.pdf", width = 1.7, height = 4)
ggarrange(pPAX5act + xlab ("") + theme(plot.margin = unit(c(0,0,0,0), "lines")), pEOMESact + theme(plot.margin = unit(c(0,0,0,0), "lines")), ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

pdf("figures/manuscript/p4.3.2.pdf", width = 1.7, height = 4)
ggarrange(pPAX5exp + xlab ("") + ylab("") + theme(plot.margin = unit(c(0,0,0,0), "lines")), pEOMESexp + ylab("") + theme(plot.margin = unit(c(0,0,0,0), "lines")), ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

pdf("figures/manuscript/p4.3.3.pdf", width = 1.4, height = 4)
ggarrange(PAX5_vln, PAX5_vln_2, EOMES_vln, EOMES_vln_2, common.legend = TRUE, heights = c(2, 1.7, 2, 1.7), nrow = 4, ncol = 1, legend = "none")
dev.off()

pdf("figures/manuscript/p4.4.pdf", width = 1.2, height = 3.4)
p
dev.off()
