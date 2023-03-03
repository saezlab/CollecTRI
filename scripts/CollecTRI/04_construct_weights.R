# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will construct weighted networks from
#' 1) matrixRider
#' 2) FIMO
#' 3) (RcisTarget)


library(tidyverse)

## Load data ---------------------------
# define version and load signed network
file.version <- "040722"
collecTRI_homogenized <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv")) %>%
  mutate(edge = paste(source, target, sep = "."))

matrixRider.res.1 <-  readRDS(file.path("output", file.version, "03_weighting_strategies", "matrixRider_res_1000bp.rds")) %>%
  mutate(edge = paste(source, target, sep = ".")) %>%
  arrange(edge)
matrixRider.res.10 <-  readRDS(file.path("output", file.version, "03_weighting_strategies", "matrixRider_res_10000bp.rds") )%>%
  mutate(edge = paste(source, target, sep = ".")) %>%
  arrange(edge)
FIMO.res.1 <- readRDS(file.path("output", file.version, "03_weighting_strategies", "FIMO_res_1000bp.rds")) %>%
  mutate(edge = paste(source, target, sep = ".")) %>%
  arrange(edge)
FIMO.res.1 <- FIMO.res.1[!duplicated(FIMO.res.1$edge),]
FIMO.res.10 <- readRDS(file.path("output", file.version, "03_weighting_strategies", "FIMO_res_10000bp.rds")) %>%
  mutate(edge = paste(source, target, sep = ".")) %>%
  arrange(edge)
FIMO.res.10 <- FIMO.res.10[!duplicated(FIMO.res.10$edge),]

weights_networks <- list(matrixRider.res.1, matrixRider.res.10, FIMO.res.1, FIMO.res.10)
names(weights_networks) <- c("matrixRider.res.1", "matrixRider.res.10", "FIMO.res.1", "FIMO.res.10")

dir.create(file.path("output", file.version, "04_weighted_networks"), showWarnings = FALSE)

set.seed(123)

## weight networks ---------------------------
norm_net <- map_dfc(names(weights_networks), function(idx_net){
  net <- weights_networks[[idx_net]]
  net <- left_join(net, collecTRI_homogenized %>% select(edge, weight)) %>%
    rename("mor" = weight)

  if(str_detect(idx_net, "matrixRider")){
    net_filtered <- net %>%
      mutate(bind_aff = bind_aff - min(bind_aff) + 1)

    write_csv(net_filtered %>%
                mutate(weight = bind_aff * mor) %>%
                select(source, target, weight),
              file.path("output", file.version, "04_weighted_networks", paste0(idx_net, ".raw.csv")))

    net_gene <- map_dfr(unique(net_filtered$target), function(tg){
      GRN_gene <- net_filtered %>% dplyr::filter(target == tg)

      GRN_gene %>% mutate(bind_aff_gene = (GRN_gene$bind_aff) / (max(GRN_gene$bind_aff, na.rm = T)))
    })

    write_csv(net_gene %>%
                mutate(weight = bind_aff_gene * mor) %>%
                select(source, target, weight),
              file.path("output", file.version, "04_weighted_networks", paste0(idx_net, ".geneNorm.csv")))

    net_TF <- map_dfr(unique(net_filtered$source), function(tf){
      GRN_tf <- net_filtered %>% dplyr::filter(source == tf)

      GRN_tf %>% mutate(bind_aff_tf = (GRN_tf$bind_aff) / (max(GRN_tf$bind_aff, na.rm = T)))
    })

    write_csv(net_TF %>%
                mutate(weight = bind_aff_tf * mor) %>%
                select(source, target, weight),
              file.path("output", file.version, "04_weighted_networks", paste0(idx_net, ".tfNorm.csv")))

    weights_df <- data.frame(edge = net_filtered$edge,
               raw = net_filtered$bind_aff *net_filtered$mor,
               gene = net_gene$bind_aff_gene * net_gene$mor,
               tf =  net_TF$bind_aff_tf * net_TF$mor)
    colnames(weights_df) <- paste0(idx_net, "_", colnames(weights_df))
    weights_df

  } else if (str_detect(idx_net, "FIMO")){
    net_filtered <- net %>%
      ungroup() %>%
      mutate(score = score - min(score) + 1)

    write_csv(net_filtered %>%
                mutate(weight = score * mor) %>%
                select(source, target, weight),
              file.path("output", file.version, "04_weighted_networks", paste0(idx_net, ".raw.csv")))

    net_gene <- map_dfr(unique(net_filtered$target), function(tg){
      GRN_gene <- net_filtered %>% dplyr::filter(target == tg)

      GRN_gene %>% mutate(bind_aff_gene = (GRN_gene$score) / (max(GRN_gene$score, na.rm = T)))
    })

    write_csv(net_gene %>%
                mutate(weight = bind_aff_gene * mor) %>%
                select(source, target, weight),
              file.path("output", file.version, "04_weighted_networks", paste0(idx_net, ".geneNorm.csv")))

    net_TF <- map_dfr(unique(net_filtered$source), function(tf){
      GRN_tf <- net_filtered %>% dplyr::filter(source == tf)

      GRN_tf %>% mutate(bind_aff_tf = (GRN_tf$score) / (max(GRN_tf$score, na.rm = T)))
    })

    write_csv(net_TF %>%
                mutate(weight = bind_aff_tf * mor) %>%
                select(source, target, weight),
              file.path("output", file.version, "04_weighted_networks", paste0(idx_net, ".tfNorm.csv")))

    weights_df <- data.frame(edge = net_filtered$edge,
                             raw = net_filtered$score *net_filtered$mor,
                             gene = net_gene$bind_aff_gene * net_gene$mor,
                             tf =  net_TF$bind_aff_tf * net_TF$mor)
    colnames(weights_df) <- paste0(idx_net, "_", colnames(weights_df))
    weights_df
  }
})
rownames(norm_net) <- norm_net$matrixRider.res.1_edge
norm_net <- norm_net[!str_detect(colnames(norm_net), "edge")]
norm_net <- cbind(norm_net, unweighted = sign(norm_net$matrixRider.res.1_raw))
cor(norm_net)

## filter networks ---------------------------
net_filtered <- map(names(weights_networks), function(idx_net){
  print(idx_net)
  net <- weights_networks[[idx_net]]
  net <- left_join(net, collecTRI_homogenized %>% select(edge, weight))

  if (str_detect(idx_net, "matrixRider")){
    net <- net %>%
      rename("score" = bind_aff)
  }

  net_filtered <- map(c(10,20,30), function(q_n){
    n_to_keep <- nrow(net) - round(nrow(net)/100 * q_n)
    net_filtered <- net %>%
      ungroup() %>%
      arrange(desc(score)) %>%
      slice(1:n_to_keep) %>%
      select(source, target, weight)

    write_csv(net_filtered, file.path("output", file.version, "04_weighted_networks", "filtered", paste0(idx_net, ".filtered", q_n,".csv")))
    net_filtered
  })
})

net_unweighted <- weights_networks$matrixRider.res.1 %>%
  left_join(collecTRI_homogenized %>% select(edge, weight)) %>%
  select(source, target, weight)

net_random10 <- net_unweighted[sample(nrow(net_unweighted), nrow(net) - round(nrow(net)/100 * 10)), ]
net_random20 <- net_unweighted[sample(nrow(net_unweighted), nrow(net) - round(nrow(net)/100 * 20)), ]
net_random30 <- net_unweighted[sample(nrow(net_unweighted), nrow(net) - round(nrow(net)/100 * 30)), ]

write_csv(net_random10, file.path("output", file.version, "04_weighted_networks", "filtered", "random10.csv"))
write_csv(net_random20, file.path("output", file.version, "04_weighted_networks", "filtered", "random20.csv"))
write_csv(net_random30, file.path("output", file.version, "04_weighted_networks", "filtered", "random30.csv"))
