library(decoupleRBench)
library(tidyverse)

#### NTNU networks v2.0 ####
setwd("/net/data.isilon/ag-saez/bq_smueller/NTNUdecoupleR")

# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
bexample_url <- file.path('data',"bench", "rna_expr_dorotheaA_new.rds")
bmeta_url <- file.path('data',"bench",  "rna_meta_dorotheaA_new.rds")
source_url <- file.path('data', "dorothea", "dorothea_A_new.rds")

input_path_v1 <- file.path('data', 'networks_v1', 'network_collection_v1.rds')
input_path_v2 <- file.path('data', 'networks_v2', 'network_collection_v2.rds')
input_path_v2_dbTF <- file.path('data', 'networks_dbTF_v2', 'network_collection_v2_dbTF.rds')

output_path <- file.path('output')

# Design contains statistical methods that take weights into account
design_row <-
tibble(
  set_name = c("Dorothea_A_new", "Dorothea_ABC_new"), # name of the set resource
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
  bexpr_loc = c(bexample_url, "data/bench/rna_expr_dorotheaABC_new.rds"), # benchmark data location
  bmeta_loc = c(bmeta_url, "data/bench/rna_meta_dorotheaABC_new.rds"), # metadata location
  source_loc = c(source_url, "data/dorothea/dorothea_ABC_new.rds"), # set source location
  source_col = "source", # source name of the gene set source
  target_col = "target", # target name of the set source'
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c("A"), c("A", "B", "C")) # criteria by which we wish to filter
)

network_collection_v1 <- readRDS(input_path_v1)$ExTRI_comp_unrestricted
network_collection_v2 <- readRDS(input_path_v2)$ExTRI_comp_unrestricted
network_collection_v2_dbTF <- readRDS(input_path_v2_dbTF)$ExTRI_comp_unrestricted

networks <- rbind(network_collection_v1, network_collection_v2, network_collection_v2_dbTF)


input_tibble <- design_row
# input tibble for one run with the weighted network
# and one where the weight is removed
for (i in 1:nrow(networks)){
  input_tibble <- input_tibble %>% add_row(design_row[1,] %>%
                                             mutate(set_name = str_remove(str_remove(networks$path[i], ".rds"), "data/")) %>%
                                             mutate(source_loc = networks$path[i]))
}

rna_exp_net <- rep(c("data/bench/rna_expr_NTNU_v1_unres.rds",
                         "data/bench/rna_expr_NTNU_v2_unres.rds",
                         "data/bench/rna_expr_NTNU_v2_dbTF_unres.rds"), each = 4)
rna_meta_net <- rep(c("data/bench/rna_meta_NTNU_v1_unres.rds",
                          "data/bench/rna_meta_NTNU_v2_unres.rds",
                          "data/bench/rna_meta_NTNU_v2_dbTF_unres.rds"), each = 4)

input_tibble <- input_tibble %>%
  mutate(bexpr_loc = c(input_tibble$bexpr_loc[1:2], rna_exp_net)) %>%
  mutate(bmeta_loc = c(input_tibble$bmeta_loc[1:2], rna_meta_net))


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
saveRDS(estimate@bench_res, paste0(output_path, "/estimate_rna_signed_test_networks_unres.rds"))

