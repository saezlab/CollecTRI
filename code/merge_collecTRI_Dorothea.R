library(tidyverse)

collecTRI <- read_csv("output/040722/final_dbTF/signed_collecTRI.csv")
dorothea <- read_csv("data/raw/dorothea_ABC.csv")[2:5]


collecTRI <- collecTRI %>%
  mutate(edge = paste(source, target, sep = ":"))

dorothea <- dorothea %>%
  mutate(edge = paste(source, target, sep = ":"))

sum(!collecTRI$source %in% dorothea$source)
sum(!collecTRI$edge %in% dorothea$edge)

plot.df <- data.frame(x = c(sum(unique(collecTRI$source) %in% unique(dorothea$source)),
                            sum(!unique(collecTRI$source) %in% unique(dorothea$source)),
                            sum(collecTRI$edge %in% dorothea$edge),
                            sum(!collecTRI$edge %in% dorothea$edge)),
              label = c("TF", "TF", "TF-gene", "TF-gene"),
           network = c("covered in ABD", "new",
                       "covered in ABD", "new"))

plot.df$network <- factor(plot.df$network, levels = unique(rev(plot.df$network)))
library(ggbreak)
ggplot(plot.df %>% filter(label == "TF-gene")) +
  aes(fill=network, y=x, x=label) +
  geom_bar(stat="identity", color="black", position="stack", width=0.7)+
  theme_minimal() + xlab("") + ylab("") +
  scale_fill_manual(labels = c("new", "covered in\nDorothea ABC"), values=c('#477AA3','#97CAF3')) +
  theme(text = element_text(size = 12)) + theme(legend.title = element_blank()) + coord_flip()

ggplot(plot.df, aes(fill=network, y=x, x=label)) +
  geom_bar(position="stack", stat="identity") +
  ggtitle("Studying 4 species..") +
  xlab("")
dorothea_filter <- dorothea %>%
  filter(!edge %in% collecTRI$edge)

collecTRI_doro <- rbind(collecTRI, dorothea_filter)
write_csv(collecTRI_doro[1:4], "output/040722/doro_collecTRI.csv")


knockTF_expr <- read_csv("data/knockTF_expr.csv") %>%
  column_to_rownames("...1")
knockTF_meta <- read_csv("data/knockTF_meta.csv")


msk = knockTF_meta$'logFC' < -1
mat = knockTF_expr[msk,]
obs = knockTF_meta[msk,]
TFs_bench <- obs$TF %>% unique()


edge_cTRI <- collecTRI %>% filter(target %in% TFs_bench) %>% pull(edge)
edge_dorothea <- dorothea %>% filter(confidence == "A") %>% filter(target %in% TFs_bench) %>% pull(edge)

sum(edge_cTRI %in% edge_dorothea)
sum(!edge_cTRI %in% edge_dorothea)
sum(!edge_dorothea %in% edge_cTRI)

