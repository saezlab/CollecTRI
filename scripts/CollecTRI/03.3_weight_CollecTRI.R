#' In this script we will assign the weights to CollecTRI and try different
#' normalisation strategies. To run this script you need to first calculate
#' weights using matrixRider and FIMO as descriped in 03.1_weights_matrixRider.R
#' and 03.2_weights_FIMO.R.

library(tidyverse)
set.seed(713)

## Load data ---------------------------
output.folder <- "output"
collecTRI <- read_csv( "output/CollecTRI/CollecTRI.csv") %>%
  dplyr::mutate(edge = paste(source, target, sep = "."))

weigths.files <- list.files(file.path("output", "weighted_networks", "raw"), full.names = T)

weights_networks <- map(weigths.files, read_csv)
weights_networks <- map(weights_networks, function(net){
  net %>%
    dplyr::mutate(edge = paste(source, target, sep = ".")) %>%
    dplyr::arrange(edge)
})
names(weights_networks) <- str_remove(list.files(file.path("output", "weighted_networks", "raw")), ".csv")

dir.create(file.path(output.folder, "weighted_networks"), showWarnings = FALSE)

## weight networks ---------------------------
norm_net <- map_dfc(names(weights_networks), function(idx_net){
  net <- weights_networks[[idx_net]]
  net <- left_join(net, collecTRI %>% select(edge, weight), by = "edge") %>%
    rename("mor" = weight)

  if(str_detect(idx_net, "matrixRider")){
    net_filtered <- net %>%
      mutate(bind_aff = bind_aff - min(bind_aff) + 1)

    write_csv(net_filtered %>%
                mutate(weight = bind_aff * mor) %>%
                select(source, target, weight),
              file.path(output.folder,"weighted_networks", paste0(idx_net, "_raw.csv")))

    net_gene <- map_dfr(unique(net_filtered$target), function(tg){
      GRN_gene <- net_filtered %>% dplyr::filter(target == tg)

      GRN_gene %>% mutate(bind_aff_gene = (GRN_gene$bind_aff) / (max(GRN_gene$bind_aff, na.rm = T)))
    })
    net_gene <- net_gene %>%
      dplyr::arrange(edge)

    write_csv(net_gene %>%
                mutate(weight = bind_aff_gene * mor) %>%
                select(source, target, weight),
              file.path(output.folder,"weighted_networks", paste0(idx_net, "_geneNorm.csv")))

    net_TF <- map_dfr(unique(net_filtered$source), function(tf){
      GRN_tf <- net_filtered %>% dplyr::filter(source == tf)

      GRN_tf %>% mutate(bind_aff_tf = (GRN_tf$bind_aff) / (max(GRN_tf$bind_aff, na.rm = T)))
    })
    net_TF <- net_TF %>%
      dplyr::arrange(edge)

    write_csv(net_TF %>%
                mutate(weight = bind_aff_tf * mor) %>%
                select(source, target, weight),
              file.path(output.folder,"weighted_networks", paste0(idx_net, "_tfNorm.csv")))

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
              file.path(output.folder,"weighted_networks", paste0(idx_net,"_raw.csv")))

    net_gene <- map_dfr(unique(net_filtered$target), function(tg){
      GRN_gene <- net_filtered %>% dplyr::filter(target == tg)

      GRN_gene %>% mutate(bind_aff_gene = (GRN_gene$score) / (max(GRN_gene$score, na.rm = T)))
    })

    net_gene <- net_gene %>%
      dplyr::arrange(edge)

    write_csv(net_gene %>%
                mutate(weight = bind_aff_gene * mor) %>%
                select(source, target, weight),
              file.path(output.folder,"weighted_networks", paste0(idx_net, "_geneNorm.csv")))

    net_TF <- map_dfr(unique(net_filtered$source), function(tf){
      GRN_tf <- net_filtered %>% dplyr::filter(source == tf)

      GRN_tf %>% mutate(bind_aff_tf = (GRN_tf$score) / (max(GRN_tf$score, na.rm = T)))
    })

    net_TF <- net_TF %>%
      dplyr::arrange(edge)

    write_csv(net_TF %>%
                mutate(weight = bind_aff_tf * mor) %>%
                select(source, target, weight),
              file.path(output.folder,"weighted_networks", paste0(idx_net, "_tfNorm.csv")))

    weights_df <- data.frame(edge = net_filtered$edge,
                             raw = net_filtered$score *net_filtered$mor,
                             gene = net_gene$bind_aff_gene * net_gene$mor,
                             tf =  net_TF$bind_aff_tf * net_TF$mor)
    colnames(weights_df) <- paste0(idx_net, "_", colnames(weights_df))
    weights_df
  }
})


## compare correlation of weights ---------------------------
rownames(norm_net) <- norm_net$matrixRider_10000bp_edge
norm_net <- norm_net[!str_detect(colnames(norm_net), "edge")]
norm_net <- cbind(norm_net, unweighted = sign(norm_net$matrixRider_1000bp_raw))

write_csv(cor(norm_net) %>% as.data.frame(),
          file.path(output.folder,"weighted_networks", "correlation_matrix.csv"))


## filter networks ---------------------------
net_filtered <- map(names(weights_networks), function(idx_net){
  net <- weights_networks[[idx_net]]
  net <- left_join(net, collecTRI %>% select(edge, weight), by = "edge")

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

    write_csv(net_filtered, file.path(output.folder, "weighted_networks", paste0(idx_net, "_filtered", q_n,".csv")))
    net_filtered
  })
})

net_unweighted <- weights_networks$matrixRider_1000bp %>%
  left_join(collecTRI %>% select(edge, weight), by = "edge") %>%
  select(source, target, weight)

net_random10 <- net_unweighted[sample(nrow(net_unweighted), nrow(net_unweighted) - round(nrow(net_unweighted)/100 * 10)), ]
net_random20 <- net_unweighted[sample(nrow(net_unweighted), nrow(net_unweighted) - round(nrow(net_unweighted)/100 * 20)), ]
net_random30 <- net_unweighted[sample(nrow(net_unweighted), nrow(net_unweighted) - round(nrow(net_unweighted)/100 * 30)), ]

write_csv(net_random10, file.path(output.folder, "weighted_networks", "random10.csv"))
write_csv(net_random20, file.path(output.folder, "weighted_networks", "random20.csv"))
write_csv(net_random30, file.path(output.folder, "weighted_networks", "random30.csv"))

