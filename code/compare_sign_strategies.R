# Copyright (c) Sophia Müller-Dott [2022]
# sophia.mueller-dott@uni-heidelberg.de

#' In this script we compare the results of the different signing startegies
#'

library(tidyverse)
## ---------------------------
files_dbTF <- list.files("output/040722/signed_networks", pattern = "dbTF", full.names = T)
files <- list.files("output/040722/signed_networks", full.names = T)[!list.files("output/040722/signed_networks", full.names = T) %in% files_dbTF]

## ---------------------------
networks <- map(files, read.csv)
names(networks) <- str_remove(str_remove(files, "output/040722/signed_networks/"), "_signed_CollecTRI_040722.csv")

networks_sign <- map(names(networks), function(sign_dec){
  net <- networks[[sign_dec]]
   net.res <- net %>% mutate(TF.TG = paste(source, target, sep = ".")) %>%
     select(TF.TG, weight)
   colnames(net.res) <- c("TF.TG", sign_dec)
   net.res
})

sign_collapse <- reduce(networks_sign, full_join)
sign_collapse_long <- sign_collapse %>%
  pivot_longer(!TF.TG, names_to = "sign_strategy", values_to = "sign") %>%
  group_by(sign_strategy) %>%
  summarise(counts_pos = sum(sign == 1),
            counts_neg = sum(sign == -1)) %>%
  pivot_longer(!sign_strategy, names_to = "sign", values_to = "counts")

ggplot(data=sign_collapse_long,
            aes(x=sign_strategy, y=counts, fill=sign)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

cor(sign_collapse %>% column_to_rownames("TF.TG"))

## ---------------------------
networks_dbTF <- map(files_dbTF, read.csv)
names(networks_dbTF) <- str_remove(str_remove(files_dbTF, "output/040722/signed_networks/"), "_signed_CollecTRI_dbTF_040722.csv")

networks_dbTF_sign <- map(names(networks_dbTF), function(sign_dec){
  net <- networks_dbTF[[sign_dec]]
  net.res <- net %>% mutate(TF.TG = paste(source, target, sep = ".")) %>%
    select(TF.TG, weight)
  colnames(net.res) <- c("TF.TG", sign_dec)
  net.res
})

sign_collapse_dbTF <- reduce(networks_dbTF_sign, full_join)
sign_collapse_dbTF_long <- sign_collapse_dbTF %>%
  pivot_longer(!TF.TG, names_to = "sign_strategy", values_to = "sign") %>%
  group_by(sign_strategy) %>%
  summarise(counts_pos = sum(sign == 1),
            counts_neg = sum(sign == -1)) %>%
  pivot_longer(!sign_strategy, names_to = "sign", values_to = "counts")

ggplot(data=sign_collapse_dbTF_long,
       aes(x=sign_strategy, y=counts, fill=sign)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

cor(sign_collapse_dbTF %>% column_to_rownames("TF.TG"))
