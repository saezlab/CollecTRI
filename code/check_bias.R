library(tidyverse)

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


# Test bias of methods ------------------------------------------------
all_bias_plots <- function(decoupler_res, n_permutations, method_filter){
  bias_ntargets_p <- map(names(decoupler_res), function(network){
    net <- final_networks[[network]]
    decoupler_res_df <- decoupler_res[[network]]
    nTargets <- table(net$source) %>%
      as.data.frame() %>%
      rename("source" = Var1)

    decoup <- left_join(decoupler_res_df, nTargets, by = "source")
    decoup <- decoup %>% filter(!(is.na(Freq)))

    decoup_per_exp <- decoup %>% filter(statistic == method_filter) %>% group_by(condition) %>% group_split()

    correlations_df <- map_dfr(decoup_per_exp, function(exp){
      all_correlation <- c()
      for (i in c(nrow(exp), 100, 75, 50, 25)){
        exp_filtered <- exp %>% arrange(desc(abs(score))) %>% slice(1:i)
        corr_pearson <- cor(abs(exp_filtered$score), exp_filtered$Freq, method = "pearson")
        corr_spearman <- cor(abs(exp_filtered$score), exp_filtered$Freq, method = "spearman")

        all_correlation <- append(all_correlation, c(corr_pearson, corr_spearman))
      }
      names(all_correlation) <- c("corr_pearson_all", "corr_spearman_all",
                                  "corr_pearson_100", "corr_spearman_100",
                                  "corr_pearson_75", "corr_spearman_75",
                                  "corr_pearson_50", "corr_spearman_50",
                                  "corr_pearson_25", "corr_spearman_25")
      all_correlation
    })

    correlations_df <- correlations_df %>%
      add_column(experiment = unique(decoup$condition)) %>%
      column_to_rownames("experiment")


    correltion_plots <- map(colnames(correlations_df), function(corr_type){
      ggplot(correlations_df, aes_string(x=corr_type)) +
        geom_histogram(color="black", fill="grey", binwidth = 0.02) +
        xlim(-1,1) +
        ggtitle(paste0(network, " ", corr_type, " = ", round(mean(correlations_df[,corr_type]), digits = 2))) +
        theme_bw()
    })

    cowplot::plot_grid(plotlist = correltion_plots, ncol = 2)

  })

  pdf(paste0("figures/final_comp/new/bias/bias_regulonsize_", method_filter, "_corr_", n_permutations, ".pdf"), width = 10, height = 20)
  print(bias_ntargets_p)
  dev.off()


  # Check individual experiments
  bias_scatterPerExp_p <- map(names(decoupler_res), function(network){
    net <- final_networks[[network]]
    decoupler_res_df <- decoupler_res[[network]]
    nTargets <- table(net$source) %>%
      as.data.frame() %>%
      rename("source" = Var1)

    decoup <- left_join(decoupler_res_df, nTargets, by = "source")
    decoup <- decoup %>% filter(!(is.na(Freq)))

    decoup_per_exp <- decoup %>% filter(statistic == method_filter) %>% group_by(condition) %>% group_split()

    i <- sample(1:length(decoup_per_exp), 6)

    scatter_plots <- map(i, function(x){
      ggplot(decoup_per_exp[[x]], aes(x=Freq, y = abs(score))) +
        geom_point() +
        theme_bw() +
        ggtitle(paste0(network, " ", x, " , corr =",
                       round(cor(decoup_per_exp[[x]]$Freq, abs(decoup_per_exp[[x]]$score)), digits = 2)
        )
        )
    })

    cowplot::plot_grid(plotlist = scatter_plots, ncol = 3)

  })

  pdf(paste0("figures/final_comp/new/bias/bias_scatterPerExp_", method_filter, "_", n_permutations, ".pdf"), width = 15, height = 10)
  print(bias_scatterPerExp_p)
  dev.off()

  # Check mean activity per TF
  bias_meanAct_p <- map(names(decoupler_res), function(network){
    net <- final_networks[[network]]
    decoupler_res_df <- decoupler_res[[network]]
    nTargets <- table(net$source) %>%
      as.data.frame() %>%
      rename("source" = Var1)

    decoup <- left_join(decoupler_res_df, nTargets, by = "source")
    decoup <- decoup %>% filter(!(is.na(Freq)))

    decoup_per_tf <- decoup %>% filter(statistic == method_filter) %>% group_by(source) %>% summarise(mean_act = mean(abs(score)),
                                                                                                      targets = unique(Freq))

    ggplot(decoup_per_tf, aes(x=targets, y = mean_act)) +
      geom_point() +
      ggtitle(paste0(network, ", pearson_corr = ",  round(cor(decoup_per_tf$mean_act, decoup_per_tf$targets), digits = 2))) +
      theme_bw()

  })

  pdf(paste0("figures/final_comp/new/bias/bias_meanAct_", method_filter, "_", n_permutations, ".pdf"), width = 10, height = 20)
  print(cowplot::plot_grid(plotlist = bias_meanAct_p, ncol = 2))
  dev.off()

  # Check activity for specific TFs

  bias_perTF_p <- map(names(decoupler_res), function(network){
    net <- final_networks[[network]]
    decoupler_res_df <- decoupler_res[[network]]
    nTargets <- table(net$source) %>%
      as.data.frame() %>%
      rename("source" = Var1)

    decoup <- left_join(decoupler_res_df, nTargets, by = "source")
    decoup <- decoup %>% filter(!(is.na(Freq)))

    decoup <- decoup %>% mutate(Freq_range = case_when(Freq < 10 ~ "< 10",
                                                       Freq >= 10 & Freq < 20 ~ "10-19",
                                                       Freq >= 20 & Freq < 30 ~ "20-29",
                                                       Freq >= 30 & Freq < 50 ~ "30-49",
                                                       Freq >= 50 & Freq < 100 ~ "50-99",
                                                       Freq >= 100 & Freq < 200 ~ "100-199",
                                                       Freq >= 200  ~ "> 200"))
    TFs_freq_range <- decoup %>% filter(statistic == method_filter) %>% group_split(Freq_range)
    set.seed(123)
    TFs_random <- map(TFs_freq_range, function(x){
      sample(unique(x$source), 4)})  %>% unlist()

    boxplot_df <- decoup %>% filter(statistic == method_filter) %>% filter(source %in% TFs_random) %>% arrange(Freq)
    boxplot_df$Freq_range <- factor(boxplot_df$Freq_range, levels = unique(boxplot_df$Freq_range))

    ggplot(boxplot_df, aes(x=Freq_range, y=abs(score), fill=source)) +
      geom_boxplot() + ggtitle(network) +
      theme_bw()
  })

  pdf(paste0("figures/final_comp/new/bias/bias_perTF_", method_filter, "_", n_permutations, ".pdf"), width = 20, height = 40)
  print(cowplot::plot_grid(plotlist = bias_perTF_p, ncol = 2))
  dev.off()
}


all_bias_plots(decoupler_res = decoupler_res_100, n_permutations = "100", method_filter = "ulm")

all_bias_plots(decoupler_res = decoupler_res_100, n_permutations = "100", method_filter = "mlm")

all_bias_plots(decoupler_res = decoupler_res_100, n_permutations = "100", method_filter = "norm_wsum")
all_bias_plots(decoupler_res = decoupler_res_1000, n_permutations = "1000", method_filter = "norm_wsum")
all_bias_plots(decoupler_res = decoupler_res_2000, n_permutations = "2000", method_filter = "norm_wsum")

all_bias_plots(decoupler_res = decoupler_res_100, n_permutations = "100", method_filter = "consensus")
all_bias_plots(decoupler_res = decoupler_res_1000, n_permutations = "1000", method_filter = "consensus")
all_bias_plots(decoupler_res = decoupler_res_2000, n_permutations = "2000", method_filter = "consensus")
