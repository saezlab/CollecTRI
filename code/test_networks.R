library(decoupleRBench)
library(dplyr)
library(purrr)
library(tidyverse)

#### NTNU networks v2.0 ####

# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
bexample_url <- file.path('data', "rna_expr.rds")
bmeta_url <- file.path('data', "rna_meta.rds")
source_url <- file.path('data', "dorothea_filtered.rds")

# Design contains statistical methods that take weights into account
design_row <-
  tibble(
    set_name = "Dorothea_A", # name of the set resource
    bench_name = "dbd", # name of the benchmark data
    stats_list = list( # a list of the stats to call
      c(
        "udt",
        "mdt",
        "wsum",
        "wmean",
        "ulm",
        "mlm",
        "viper"
      )
    ),
    opts_list = list(list( # list of options for each stat method
      udt = list(.mor = "mor",
                 .likelihood = "likelihood"),
      mdt = list(.mor = "mor",
                 .likelihood = "likelihood"),
      wsum = list(times = 100,
                  .mor = "mor",
                  .likelihood = "likelihood"),
      wmean = list(times = 100,
                   .mor = "mor",
                   .likelihood = "likelihood"),
      ulm = list(.mor = "mor",
                 .likelihood = "likelihood"),
      mlm = list(.mor = "mor",
                 .likelihood = "likelihood"),
      viper = list(verbose = FALSE,
                   minsize = 0,
                   .mor = "mor",
                   .likelihood = "likelihood",
                   pleiotropy = TRUE,
                   eset.filter = FALSE)
    )),
    bexpr_loc = bexample_url, # benchmark data location
    bmeta_loc = bmeta_url, # metadata location
    source_loc = source_url, # set source location
    source_col = "source", # source name of the gene set source
    target_col = "target", # target name of the set source'
    filter_col = "confidence", # column by which we wish to filter
    filter_crit = list(c("A")) # criteria by which we wish to filter
  )

network_collection <- readRDS("data/network_collection.rds")

map(names(network_collection), function(name){
  networks <- network_collection[[name]]
  input_tibble <- design_row
  # input tibble for one run with the weighted network
  # and one where the weight is removed
  for (i in 1:nrow(networks)){
    input_tibble <- input_tibble %>% add_row(design_row %>%
                                               mutate(set_name = str_remove(str_remove(networks$path[i], ".rds"), "data/")) %>%
                                               mutate(source_loc = networks$path[i]))
  }

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

  ggsave(paste0("figures/auc/", name, ".pdf"), device = "pdf", width = 20, height = 8, plot = auc)
  ggsave(paste0("figures/prc/", name, ".pdf"), device = "pdf", width = 20, height = 8, plot = prc)

})
