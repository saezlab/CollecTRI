library(decoupleRBench)
library(dplyr)
library(purrr)
library(tidyverse)

#### NTNU networks v1.0 ####
# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
bexample_url <- file.path('data', "dorothea_bench_expr.rds")
bmeta_url <- file.path('data', "dorothea_bench_meta.rds")
source_url <- file.path('data', "GRN_agnostic.rds")

# Design contains statistical methods that take weights into account
design_row <-
  tibble(
    set_name = "GRN_agnostic", # name of the set resource
    bench_name = "dbd", # name of the benchmark data
    stats_list = list( # a list of the stats to call
      c(
        "scira",
        "pscira",
        "mean",
        "viper"
      )
    ),
    opts_list = list(list( # list of options for each stat method
      scira = list(.mor = "mor",
                   .likelihood = "likelihood"),
      pscira = list(times = 100,
                    .mor = "mor",
                    .likelihood = "likelihood"),
      mean = list(times = 100,
                  .mor = "mor",
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

# input tibble for one run with the weighted network
# and one where the weight is removed
input_tibble <- bind_rows(
  design_row,
  design_row %>%
    mutate(set_name ="GRN_weighted_evidence") %>%
    mutate(source_loc = file.path('data', "GRN_weighted_evidence.rds")),
  design_row %>%
    mutate(set_name ="GRN_signed") %>%
    mutate(source_loc = file.path('data', "GRN_signed.rds")),
  #design_row %>%
  #  mutate(set_name ="GRN_signed_0.8") %>%
  #  mutate(source_loc = file.path('data', "GRN_signed_0.8.rds")),
  design_row %>%
    mutate(set_name ="GRN_signed_weighted_evidence") %>%
    mutate(source_loc = file.path('data', "GRN_signed_weighted_evidence.rds")),
  #design_row %>%
  #  mutate(set_name ="GRN_signed_weighted_evidence_0.8") %>%
  #  mutate(source_loc = file.path('data', "GRN_signed_weighted_evidence_0.8.rds")),
  design_row %>%
    mutate(set_name ="GRN_signed_weighted_evidence_sign") %>%
    mutate(source_loc = file.path('data', "GRN_signed_weighted_evidence_sign.rds")),
  #design_row %>%
  #  mutate(set_name ="GRN_signed_weighted_evidence_sign_0.8") %>%
  #  mutate(source_loc = file.path('data', "GRN_signed_weighted_evidence_sign_0.8.rds")),
  design_row %>%
    mutate(set_name ="GRN_signed_weighted_discrepancy_sign") %>%
    mutate(source_loc = file.path('data', "GRN_signed_weighted_discrepancy_sign.rds")),
  #design_row %>%
  #  mutate(set_name ="GRN_signed_weighted_discrepancy_sign_0.8") %>%
  #  mutate(source_loc = file.path('data', "GRN_signed_weighted_discrepancy_sign_0.8.rds")),
  design_row %>%
    mutate(set_name ="Dorothea_A") %>%
    mutate(source_loc = file.path('data', "dorothea_filtered.rds"))
)


# run decoupleRBenchmark
estimate <- run_benchmark(
  .design = input_tibble, # provide input tibble
  .minsize = 5, # filter gene sets with size < 5
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = TRUE, # silently run the pipeline
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
boxplot_tibble$network <- factor(boxplot_tibble$network, levels = unique(estimate@bench_res$set_name))

ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = auc)) +
  geom_boxplot()

ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = prc)) +
  geom_boxplot()



#### NTNU networks v2.0 ####

# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
bexample_url <- file.path('data', "dorothea_bench_expr.rds")
bmeta_url <- file.path('data', "dorothea_bench_meta.rds")
source_url <- file.path('data', "dorothea_filtered.rds")

# Design contains statistical methods that take weights into account
design_row <-
  tibble(
    set_name = "Dorothea_A", # name of the set resource
    bench_name = "dbd", # name of the benchmark data
    stats_list = list( # a list of the stats to call
      c(
        "scira",
        "pscira",
        "mean",
        "viper"
      )
    ),
    opts_list = list(list( # list of options for each stat method
      scira = list(.mor = "mor",
                   .likelihood = "likelihood"),
      pscira = list(times = 100,
                    .mor = "mor",
                    .likelihood = "likelihood"),
      mean = list(times = 100,
                  .mor = "mor",
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
    .silent = TRUE, # silently run the pipeline
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
    geom_boxplot()

  roc <- ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = prc)) +
    geom_boxplot()

  ggsave(paste0("figures/auc/", name, ".png"), plot = auc, width = 1750, height = 800, units = "px")
  ggsave(paste0("figures/roc/", name, ".png"), plot = roc, width = 1750, height = 800, units = "px")

})


#### old ####
networks <- network_collection[["final_comp"]]

input_tibble <- design_row
# input tibble for one run with the weighted network
# and one where the weight is removed
for (i in 1:nrow(networks)){
  input_tibble <- input_tibble %>% add_row(design_row %>%
                        mutate(set_name =paste(networks$GRNcuration[i], networks$curated_counts[i], networks$mor[i], networks$weight[i], sep = "_")) %>%
                        mutate(source_loc = file.path('data', paste0(networks$GRNcuration[i], "_",  networks$curated_counts[i], "_", networks$mor[i], "_", networks$weight[i], ".rds"))))
}

input_tibble <- input_tibble %>% add_row(design_row %>%
                                           mutate(set_name = "GRNcuration_ExTRI_scaled_FALSE_evidence_unfiltered") %>%
                                           mutate(source_loc = file.path('data',"GRNcuration_ExTRI_scaled_FALSE_scaledunfiltered.rds"))) %>%
  add_row(design_row %>%
            mutate(set_name = "GRNcuration_ExTRI_scaled_FALSE_none_unfiltered") %>%
            mutate(source_loc = file.path('data',"GRNcuration_ExTRI_scaled_FALSE_none_unfiltered.rds")))


# run decoupleRBenchmark
estimate <- run_benchmark(
  .design = input_tibble, # provide input tibble
  .minsize = 5, # filter gene sets with size < 5
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = TRUE, # silently run the pipeline
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

boxplot_tibble <- boxplot_tibble %>% filter(!network == "GRNcuration_ExTRI_scaled_FALSE_none_unfiltered") %>% filter(!network == "GRNcuration_new_scaled_TRUE_evidence")

boxplot_tibble <- boxplot_tibble %>% mutate(network = replace(network, network == "GRNcuration_ExTRI_scaled_FALSE_none" , "GRN_agnostic"),
                                            network = replace(network, network == "GRNcuration_ExTRI_scaled_FALSE_evidence" , "GRN_weighted"),
                                            network = replace(network, network == "GRNcuration_ExTRI_scaled_TRUE_evidence" , "GRN_signed_weighted"),
                                            network = replace(network, network == "GRNcuration_ExTRI_tot_standard_TRUE_evidence_tot" , "GRN_total_signed_weighted"),
                                            network = replace(network, network == "GRNcuration_ExTRI_scaled_FALSE_evidence_unfiltered" , "GRN_weighted_unfiltered"))


boxplot_tibble$network <- factor(boxplot_tibble$network, levels = boxplot_tibble %>%
                                   group_by(network) %>%
                                   summarise(mean = mean(auc)) %>% arrange(mean) %>% pull(network))


ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = auc)) +
  geom_boxplot()

ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = prc)) +
  geom_boxplot()


## Compare networks with same number of TFs
input_tibble_size <- design_row %>% add_row(design_row %>%
                                              mutate(set_name = "ExTRI_signed_weighted_resticted") %>%
                                              mutate(source_loc = "data/ExTRI_signed_weighted_restricted.rds")) %>%
  add_row(design_row %>%
            mutate(set_name = "ExTRI_weighted_unfiltered_resticted") %>%
            mutate(source_loc = "data/ExTRI_weighted_unfiltered_restricted.rds"))

input_tibble_size <- design_row %>% add_row(design_row %>%
                                              mutate(set_name = "ExTRI_signed_weighted") %>%
                                              mutate(source_loc = "data/GRNcuration_ExTRI_scaled_TRUE_evidence.rds")) %>%
  add_row(design_row %>%
            mutate(set_name = "ExTRI_signed_weighted_unrestricted") %>%
            mutate(source_loc = "data/GRN_tmp.rds"))


# run decoupleRBenchmark
estimate <- run_benchmark(
  .design = input_tibble_size, # provide input tibble
  .minsize = 5, # filter gene sets with size < 5
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = TRUE, # silently run the pipeline
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

boxplot_tibble <- boxplot_tibble %>% filter(!network == "GRNcuration_ExTRI_scaled_FALSE_none_unfiltered") %>% filter(!network == "GRNcuration_new_scaled_TRUE_evidence")

boxplot_tibble <- boxplot_tibble %>% mutate(network = replace(network, network == "GRNcuration_ExTRI_scaled_FALSE_none" , "GRN_agnostic"),
                                            network = replace(network, network == "GRNcuration_ExTRI_scaled_FALSE_evidence" , "GRN_weighted"),
                                            network = replace(network, network == "GRNcuration_ExTRI_scaled_TRUE_evidence" , "GRN_signed_weighted"),
                                            network = replace(network, network == "GRNcuration_ExTRI_tot_standard_TRUE_evidence_tot" , "GRN_total_signed_weighted"),
                                            network = replace(network, network == "GRNcuration_ExTRI_scaled_FALSE_evidence_unfiltered" , "GRN_weighted_unfiltered"))


boxplot_tibble$network <- factor(boxplot_tibble$network, levels = boxplot_tibble %>%
                                   group_by(network) %>%
                                   summarise(mean = mean(auc)) %>% arrange(mean) %>% pull(network))


ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = auc)) +
  geom_boxplot()

ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = prc)) +
  geom_boxplot()
