library(tidyverse)
library(ggbreak)

file.version <- "040722"

collecTRI <- read_csv(file.path("output", file.version, "02_signed_networks", "strict_signed_collecTRI.csv"))
dorothea <- read_csv("data/raw/dorothea_ABC.csv")[2:5]
#mat <- read_csv('data/knockTF_expr.csv')
obs <- read_csv('data/knockTF_meta.csv')
msk <- obs$logFC < -1
obs <- obs[msk,]

TFs_bench <- obs$TF %>% unique()


collecTRI <- collecTRI %>%
  mutate(edge = paste(source, target, sep = ":"))

dorothea <- dorothea %>%
  mutate(edge = paste(source, target, sep = ":")) %>%
  mutate(weight = case_when(
    weight > 0 ~ 1,
    weight < 0 ~ -1
  ))

merged <- full_join(collecTRI, dorothea, by = c("edge", "source", "target"))
merged <- merged %>%
  mutate(weight = case_when(
    !is.na(weight.x) & decision != "unknown" ~ weight.x,
    !is.na(weight.x) & is.na(weight.y) & decision == "unknown" ~ weight.x,
    is.na(weight.x) ~  weight.y,
    decision == "unknown" & !is.na(weight.y) ~ weight.y

  ))

shared_TFs <- unique(collecTRI$source)[unique(collecTRI$source) %in% unique(dorothea$source)]

n_shared <- sum(rowSums(is.na(merged)) == 0)
n_doro <- nrow(dorothea)-n_shared
n_collecTRI <- nrow(collecTRI)-n_shared

n_shared_TF <- sum(unique(collecTRI$source) %in% unique(dorothea$source))
n_doro_TF <- length(unique(dorothea$source))-n_shared_TF
n_collecTRI_TF <- length(unique(collecTRI$source))-n_shared_TF
n_shared_TF_bench <- sum(shared_TFs %in% TFs_bench)
n_doro_TF_bench <- sum(unique(dorothea$source)[!unique(dorothea$source) %in% shared_TFs] %in% TFs_bench)
n_collecTRI_TF_bench <- sum(unique(collecTRI$source)[!unique(collecTRI$source) %in% shared_TFs] %in% TFs_bench)



plot.df <- data.frame(x = c(n_shared,
                            n_doro,
                            n_collecTRI,
                            n_shared_TF,
                            n_doro_TF,
                            n_collecTRI_TF,
                            n_shared_TF_bench,
                            n_doro_TF_bench,
                            n_collecTRI_TF_bench),
              network = rep(c("shared", "dorothea", "collecTRI"), times = 3),
              label = rep(c("edge", "TF", "TF in bench"), each = 3))

ggplot(plot.df %>% filter(!label == "edge")) +
  aes(fill=network, y=x, x=label) +
  geom_bar(stat="identity", color="black", position="stack", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  theme(text = element_text(size = 12)) + theme(legend.title = element_blank())

ggplot(plot.df %>% filter(label == "edge")) +
  aes(fill=network, y=x, x=label) +
  geom_bar(stat="identity", color="black", position="stack", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  theme(text = element_text(size = 12)) + theme(legend.title = element_blank())

dir.create(file.path("output", file.version, "05_merged_network"), showWarnings = FALSE)
write_csv(merged %>% select(source, target, weight, TF.category, PMID), file.path("output", file.version, "05_merged_network", "doro_collecTRI.csv"))
