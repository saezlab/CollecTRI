library(decoupleRBench)
library(tidyverse)

#### NTNU networks v2.0 ####
setwd("/net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR")

# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
bexample_url <- file.path('data',"bench", "knockTF_expr_dorotheaA.rds")
bmeta_url <- file.path('data',"bench",  "knockTF_meta_dorotheaA.rds")
source_url <- file.path('data', "dorothea", "dorothea_A.rds")
input_path <- file.path('data', 'networks_v2', 'network_collection_v2.rds')
output_path <- file.path('figures', 'final_comp', 'knockTF', 'abs')
manual_collection <- TRUE

# Design contains statistical methods that take weights into account
design_row <-
  tibble(
    set_name = c("Dorothea_A", "Dorothea_ABC"), # name of the set resource
    bench_name = "dbd", # name of the benchmark data
    stats_list = list( # a list of the stats to call
      c(
        "wsum",
        "ulm",
        "mlm"
      )
    ),
    opts_list = list(list( # list of options for each stat method
      wsum = list(times = 100,
                  .mor = "mor"),
      ulm = list(.mor = "mor"),
      mlm = list(.mor = "mor")
    )),
    bexpr_loc = bexample_url, # benchmark data location
    bmeta_loc = bmeta_url, # metadata location
    source_loc = c(source_url, "data/dorothea/dorothea_ABC.rds"), # set source location
    source_col = "source", # source name of the gene set source
    target_col = "target", # target name of the set source'
    filter_col = "confidence", # column by which we wish to filter
    filter_crit = list(c("A"), c("A", "B", "C")) # criteria by which we wish to filter
  )

network_collection <- readRDS(input_path)

if(manual_collection){
  network_collection <- list(ExTRI_comparison = data.frame(path = c("data/networks_v1/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                                                                    "data/networks_v1/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                                                                    "data/networks_v1/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                                                                    "data/networks_v1/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                                                                    "data/networks_v2/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                                                                    "data/networks_v2/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                                                                    "data/networks_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                                                                    "data/networks_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds",
                                                                    "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_none_TRUE.rds",
                                                                    "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_none_TRUE.rds",
                                                                    "data/networks_dbTF_v2/ExTRI_comp_scaled_FALSE_1_evidence_TRUE.rds",
                                                                    "data/networks_dbTF_v2/ExTRI_comp_scaled_TRUE_1_evidence_TRUE.rds")))
}

map(names(network_collection), function(name){
  networks <- network_collection[[name]]
  input_tibble <- design_row
  # input tibble for one run with the weighted network
  # and one where the weight is removed
  for (i in 1:nrow(networks)){
    input_tibble <- input_tibble %>% add_row(design_row[1,] %>%
                                               mutate(set_name = str_remove(str_remove(networks$path[i], ".rds"), "data/")) %>%
                                               mutate(source_loc = networks$path[i]))
  }
  input_tibble <- input_tibble %>%
    mutate(bexpr_loc = c("data/bench/knockTF_expr_dorotheaA.rds", "data/bench/knockTF_expr_dorotheaABC.rds",
                        "data/bench/knockTF_expr_NTNUv1.rds", "data/bench/knockTF_expr_NTNUv1.rds", "data/bench/knockTF_expr_NTNUv1.rds", "data/bench/knockTF_expr_NTNUv1.rds",
                        "data/bench/knockTF_expr_NTNUv2.rds", "data/bench/knockTF_expr_NTNUv2.rds", "data/bench/knockTF_expr_NTNUv2.rds", "data/bench/knockTF_expr_NTNUv2.rds",
                        "data/bench/knockTF_expr_NTNUv2_dbTF.rds", "data/bench/knockTF_expr_NTNUv2_dbTF.rds", "data/bench/knockTF_expr_NTNUv2_dbTF.rds", "data/bench/knockTF_expr_NTNUv2_dbTF.rds")) %>%
    mutate(bmeta_loc = c("data/bench/knockTF_meta_dorotheaA.rds", "data/bench/knockTF_meta_dorotheaABC.rds",
                         "data/bench/knockTF_meta_NTNUv1.rds", "data/bench/knockTF_meta_NTNUv1.rds", "data/bench/knockTF_meta_NTNUv1.rds", "data/bench/knockTF_meta_NTNUv1.rds",
                         "data/bench/knockTF_meta_NTNUv2.rds", "data/bench/knockTF_meta_NTNUv2.rds", "data/bench/knockTF_meta_NTNUv2.rds", "data/bench/knockTF_meta_NTNUv2.rds",
                         "data/bench/knockTF_meta_NTNUv2_dbTF.rds", "data/bench/knockTF_meta_NTNUv2_dbTF.rds",  "data/bench/knockTF_meta_NTNUv2_dbTF.rds",  "data/bench/knockTF_meta_NTNUv2_dbTF.rds"))


  # run decoupleRBenchmark
  estimate <- run_benchmark(
    .design = input_tibble, # provide input tibble
    .minsize = 5, # filter gene sets with size < 5
    .form = TRUE, # format the benchmark results
    .perform = TRUE, # evaluate benchmarking performance
    .silent = FALSE, # silently run the pipeline
    .downsample_pr = TRUE, # downsample TNs for precision-recall curve
    .downsample_roc = TRUE, # downsample TNs for ROC
    .downsample_times = 100, # downsampling iterations
    .url_bool = FALSE # whether to load from url
  )
saveRDS(estimate, paste0(output_path, "/estimate_knockTF_abs.rds"))
  # extract each auc per permutation run
  auc_downsampling <- lapply(estimate@bench_res$roc,
                             function(x) x %>% group_by(run) %>% summarize(raw_auc = unique(raw_auc)) %>% pull(raw_auc))

  estimate@bench_res <- add_column(estimate@bench_res, auc_downsampling)

  prc_downsampling <- lapply(estimate@bench_res$prc,
                             function(x) x %>% group_by(run) %>% summarize(raw_auc = unique(raw_auc)) %>% pull(raw_auc))

  estimate@bench_res <- add_column(estimate@bench_res, prc_downsampling)

  # boxplot sorted by methods with lowest p-value (t-test)
  boxplot_tibble <- bind_cols(auc = unlist(auc_downsampling),
                              prc = unlist(prc_downsampling),
                              statistic = rep(estimate@bench_res$statistic, each = 100),
                              network = rep(estimate@bench_res$set_name, each = 100))

  boxplot_tibble$network <- factor(boxplot_tibble$network, levels = boxplot_tibble %>%
                                     group_by(network) %>%
                                     summarise(mean = mean(auc)) %>% arrange(mean) %>% pull(network))


  auc <- ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = auc)) +
    geom_boxplot() + theme_grey(base_size = 14)

  prc <- ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = prc)) +
    geom_boxplot() + theme_grey(base_size = 14)

saveRDS(boxplot_tibble, paste0(output_path, "/bench_res_knockTF_abs.rds"))
  ggsave(paste0(output_path, "/auc/", name, ".pdf"), device = "pdf", width = 20, height = 8, plot = auc)
  ggsave(paste0(output_path, "/prc/", name, ".pdf"), device = "pdf", width = 20, height = 8, plot = prc)
  ggsave(paste0(output_path, "/auc/", name, ".png"), device = "png", width = 20, height = 8, plot = auc)
  ggsave(paste0(output_path, "/prc/", name, ".png"), device = "png", width = 20, height = 8, plot = prc)

})

system("say done")
