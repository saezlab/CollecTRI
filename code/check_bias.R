# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we will check if TFs with more
#' targets by default get a higher score in the
#' prediction and thus automatically perform better overall

library(tidyverse)
library(ggsignif)
library(ggpubr)
library(rstatix)

## Load networks and TF activities ---------------------------
file.version <- "040722"

# networks
collecTRI <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))
dorothea <- read_csv("data/raw/dorothea_ABC.csv")[2:5]
regnet <- read_csv("data/raw/regnetwork.csv")

obs <- read_csv('data/knockTF_meta.csv')
msk <- obs$logFC < -1
obs_filtered <- obs[msk,]

TFs_bench <- obs_filtered$TF %>% unique()


# TF activities
decoupler_path <- list.files(file.path("output", file.version, "decoupler"), full.names = T)
act <-  map(decoupler_path, read_csv)
names(act) <-  list.files(file.path("output", file.version, "decoupler"))

## Check difference of TFs in benchmark ---------------------------
TFs_collecTRI <- table(collecTRI$source) %>% as.data.frame() %>%
  add_column(network = "collecTRI")
TFs_doro <- table(dorothea$source) %>% as.data.frame() %>%
  add_column(network = "dorothea")
TFs_regnet <- table(regnet$source) %>% as.data.frame() %>%
  add_column(network = "regnet")

TF_df <- rbind(TFs_collecTRI, TFs_doro, TFs_regnet) %>%
  mutate(benchmark = case_when(
    Var1 %in% TFs_bench ~ "TF in benchmark",
    !Var1 %in% TFs_bench ~ "TF not in benchmark"
  )) %>%
  rename("source" = "Var1") %>%
  rename("nTargets" = "Freq") %>%
  filter(nTargets >= 5)

TF_df$network <- as.factor(TF_df$network)

stat.test <- TF_df %>%
  group_by(network) %>%
  t_test(nTargets ~ benchmark) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

bxp <- ggplot(TF_df, aes(x=network, y=nTargets, color=benchmark)) +
  geom_boxplot()

bxp <- ggboxplot(
  TF_df, x = "network", y = "nTargets", fill = "benchmark"
) +
  theme_minimal()

stat.test <- stat.test %>%
  add_xy_position(x = "network", dodge = 0.8)
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0.01,
  bracket.nudge.y = -2
)





TF_df %>%
  group_by(network, benchmark) %>%
  summarise(mean_targets = mean(nTargets)) %>%
  pivot_wider(names_from = benchmark, values_from = mean_targets) %>%
  mutate(ratio =`TF in benchmark`/`TF not in benchmark`)


tmp <- ggplot(TF_df, aes(x=network, y=nTargets, fill=benchmark)) +
 geom_boxplot() +
  ylab("Number of targets")

stats.test <- TF_df %>%
  group_by(network) %>%
  t_test(nTargets ~ benchmark) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


stat.test <- stats.test %>%
  add_xy_position(x = "network", dodge = 0.8)
tmp +
  stat_pvalue_manual(
    stat.test, label = "p.adj", tip.length = 0.01,
    bracket.nudge.y = -2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


ggplot(TF_df %>% filter(benchmark == "TF not in benchmark"), aes(x=network, y=nTargets, fill=benchmark)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c( "collecTRI", "regnet")),
              map_signif_level=TRUE)

ggplot(TF_df %>% filter(network == "dorothea"), aes(x=benchmark, y=Freq, fill=network)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("TF in benchmark", "TF not in benchmark")),
              map_signif_level=TRUE)


## Test bias ------------------------------------------------
# correlation per experiment
split_func <- function(x, by) {
  r <- diff(range(x))
  out <- seq(0, r - by - 1, by = by)
  c(round(min(x) + c(0, out - 0.51 + (max(x) - max(out)) / 2), 0), max(x))
}


test2 <- act_df %>% filter(network == "collecTRI") %>% arrange(desc(Freq))

test2$source <- factor(test2$source, levels = unique(test2$source))

ggplot(test2, aes(x = source, y = abs(act))) +
  geom_boxplot()
test2 %>% group_by(source) %>% summarise(targets = mean(Freq))

test <- map_dfr(names(act), function(act_i){
  act_df <- act[[act_i]] %>%
    pivot_longer(!...1,
                 names_to = "source",
                 values_to = "act"
    )

  if(str_detect(string = act_i, pattern = "collecTRI")){
    act_df <- act_df %>% left_join(TF_df %>% filter(network == "collecTRI"), by = "source")
  } else if (str_detect(string = act_i, pattern = "doro")){
    act_df <- act_df %>% left_join(TF_df %>% filter(network == "dorothea"), by = "source")
  }

  map_dfr(unique(act_df$...1), function(exp){
    df_exp <- act_df %>% filter(...1 == exp)
    data.frame(experiment = exp,
               pearson.cor = cor(abs(df_exp$act), df_exp$Freq, method = "pearson"),
               network = act_i)
  })
})

ggplot(test, aes(x = network, y = pearson.cor)) +
  geom_boxplot()

map(names(act)[3:4], function(act_i){
  act_df <- act[[act_i]] %>% column_to_rownames("...1")
  act_summarized <- map_dfr(1:ncol(act_df), function(col_i){
    data.frame(source = colnames(act_df)[col_i],
               median = mean(abs(act_df[[col_i]])),
               std = sd(abs(act_df[[col_i]])))
  }) %>%
    arrange(desc(median))

  if(str_detect(string = act_i, pattern = "collecTRI")){
    act_summarized <- left_join(act_summarized, TF_df %>% filter(network == "collecTRI"))
  } else if (str_detect(string = act_i, pattern = "doro")){
    act_summarized <- left_join(act_summarized, TF_df %>% filter(network == "dorothea"))
  }

  ggplot(act_summarized, aes(x = log(Freq), y= std)) +
    geom_point() +
    ggtitle(paste(act_i,
                  "spearman",
                  cor(log(act_summarized$Freq), act_summarized$std, method = "spearman") %>% round(digits = 3),
                  "pearson",
                  cor(log(act_summarized$Freq), act_summarized$std, method = "pearson") %>% round(digits = 3),
                  sep = " "))


})
