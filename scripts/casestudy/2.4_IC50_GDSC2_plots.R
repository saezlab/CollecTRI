library(VennDiagram)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
library(data.table)

CollecTRI <- read.csv("data/Networks/CollecTRI.csv")
Dorothea <- read.csv("data/Networks/dorothea_ABC.csv")
CollecTRI_assoc <- read.delim("IC50_GDSC2_CollecTRI_resAll.tsv", skip = 1)
Dorothea_assoc <- read.delim("IC50_GDSC2_Dorothea_resAll.tsv", skip = 1)
IC50_ulm_CollecTRI <- read.csv("IC50_ulm_CollecTRI.csv")
IC50_ulm_Dorothea <- read.csv("IC50_ulm_Dorothea.csv")

ulm_flt <- IC50_ulm_CollecTRI %>%
  filter(Pval < 0.05)

CollecTRI_flt <- CollecTRI_assoc %>%
  filter(p.value..adj.. < 0.05) %>%
  separate(X.Model, c("Drug","TF"), remove = F, sep = ":")

length(unique(CollecTRI_flt$Drug))
length(unique(CollecTRI_flt$TF))
# Generate barplot with the TF with the most associations
tot_associations <- CollecTRI_flt  %>%
  group_by(TF) %>%
  summarize(tot_assoc = n())

top_20 <- tot_associations %>%
  top_n(20,wt = tot_assoc)

setdiff(top_20$TF,ulm_flt$TF)


tot_associations %>%
  top_n(10,wt = tot_assoc) %>%
  ggplot(aes(x = reorder(TF, tot_assoc, decreasing = T), y = tot_assoc))+
  geom_bar(stat = "identity", fill = "steelblue") +
  ylab("Number of associations") + xlab("Transcription Factor")+
  theme_minimal() +
  theme(text = element_text(size=20))

# Generate figure for correlation of regulon size and number of associations
regulon_size <- CollecTRI %>%
  group_by(source) %>%  summarize(reg = n())

TF_char <- left_join(tot_associations,regulon_size, by = c("TF" = "source"))

ggplot(TF_char, aes(x = tot_assoc, y = reg)) +
  geom_point() +
  labs(x = "Total associations with drug responses", y = "Regulon size") +
  theme_minimal()

# Calculate jaccard similarity between the regulons of TFs associated with the same drug

#Filter drugs that have more than 1 TF associated
more_than1_TF <- CollecTRI_flt %>%
  group_by(Drug) %>%
  summarise(TF_n = n()) %>%
  filter(TF_n > 1)

pairs <- CollecTRI_flt %>%
  select(TF,Drug) %>%
  filter(Drug %in% more_than1_TF$Drug) %>%
  split(~Drug)

#Generate pairwise combinations among the TFs associated with each drug
pairs <- lapply(pairs, function(x) combn(x$TF, 2, FUN = function(x) paste(x[1], x[2], sep = "_")))

combined_df <- do.call(rbind, lapply(pairs, function(df) data.frame(Combination = df)))

# Check that pairs are not repeated, e.g. A-B and B-A.
combined_df$SortedComb <- apply(combined_df, 1, function(row) paste(sort(unlist(strsplit(row, "_"))), collapse = "_"))

# Count the occurrences of each combination
pair_counts <- table(combined_df$SortedComb)

# Print the pair counts
print(pair_counts)
freq <- combined_df %>% group_by(SortedComb) %>% summarise(Tot_pairs = n()) %>% mutate(Freq = Tot_pairs/nrow(CollecTRI_flt))

in_pair <- separate(freq, SortedComb, into = c("TF_A", "TF_B"), sep = "_")

regulons <- subset(CollecTRI, source %in% c(in_pair$TF_A, in_pair$TF_B))
num_sources <- length(unique(regulons$source))
similarity_matrix <- matrix(0, nrow = num_sources, ncol = num_sources, dimnames = list(unique(regulons$source), unique(regulons$source)))

# Calculate Jaccard similarity coefficient for each pair of sources
for (i in 1:num_sources) {
  for (j in i:num_sources) {
    source1 <- unique(regulons$source)[i]
    source2 <- unique(regulons$source)[j]
    targets1 <- unique(regulons$target[regulons$source == source1])
    targets2 <- unique(regulons$target[regulons$source == source2])
    intersection <- length(intersect(targets1, targets2))
    union <- length(union(targets1, targets2))
    similarity <- intersection / union
    similarity_matrix[source1, source2] <- similarity
    similarity_matrix[source2, source1] <- similarity
  }
}

similarity_df <- as.data.frame(as.table(similarity_matrix)) %>%
  unite(SortedComb, c(Var1,Var2)) %>%
  filter(SortedComb %in% freq$SortedComb) %>%
  rename(Similarity = "Freq")

# Print the similarity matrix
print(similarity_matrix)

sim_freq <- merge(freq, similarity_df)
ggplot(sim_freq, aes(x = Similarity, y = Tot_pairs)) +
  geom_point() +
  theme_minimal() +
  xlab("Regulon similarity") + ylab("Total common associations")

# Comparison with DoRothEA ABC results
Dorothea_flt <- Dorothea_assoc %>%
  filter(p.value..adj.. < 0.05)%>%
  separate(X.Model, c("Drug","TF"), remove = F, sep = ":")

length(unique(Dorothea_flt$Drug))
length(unique(Dorothea_flt$TF))

length(intersect(CollecTRI_flt$X.Model, Dorothea_flt$X.Model))

all_activities <- inner_join(IC50_ulm_CollecTRI, IC50_ulm_Dorothea, by = c("Sample","TF"))

ggscatter(all_activities, x = "Activity.x", y = "Activity.y",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,     xlab =  "CollecTRI ulm scores", ylab = "DoRothEA ABC ulm scores",                             # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "pearson", label.x = 3, label.y = 30)

CollecTRI <- select(CollecTRI, source, target)
Dorothea <- select(Dorothea, source, target)
common_tfs <- intersect(CollecTRI_flt$TF, Dorothea_flt$TF)

# Create an empty data frame to store the results
jaccard_similarities <- data.frame(tf = character(0), similarity = numeric(0))

# Iterate over the TFs
for (tf in common_tfs) {
  set1 <- subset(CollecTRI, source == tf)
  set2 <- subset(Dorothea, source == tf)
  targets1 <- set1[,"target"]
  targets2 <- set2[,"target"]
  intersection <- length(intersect(targets1, targets2))
  union <- length(union(targets1, targets2))
  similarity <- intersection / union

  # Add the TF and its similarity to the data frame
  jaccard_similarities <- rbind(jaccard_similarities, data.frame(tf = tf, similarity = similarity))
}

jaccard_similarities %>%
 # filter(similarity > 0) %>%
  ggplot( aes(x=reorder(tf, similarity, decreasing = T), y=similarity)) +
  geom_line() +
  geom_point(shape=21, color="black", fill="#69b3a2", size=3) +
  theme_minimal() +  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size = 10)) +
  xlab("TF")

sum(jaccard_similarities$similarity < 0.1)

ggsave("jaccard_regulons.png", width = 20, height = 5)

# Plot model coefficients with CollecTRI and Dorothea
common <- inner_join(CollecTRI_assoc,Dorothea_assoc, by ="X.Model")
ggplot(common, aes(x=Coefficient.x, y=Coefficient.y)) +
  geom_point() +
  xlab("Model coefficient - CollecTRI") +
  ylab("Model coefficient - Dorothea ABC")


# Plot association for drugs targeting the same pathway

GDSC2_data <- read.delim("data/GDSC/GDSC2_fitted_dose_response.csv")

GDSC2_data_flt <- GDSC2_data %>%
  select(DRUG_ID, DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME) %>%
  rename(ID = "DRUG_ID") %>%
  mutate(ID = as.character(ID))

path_assoc <- CollecTRI_flt %>%
  group_by(Drug) %>%
  summarize(n = n()) %>%
  left_join(GDSC2_data_flt, by = c("Drug" = "ID")) %>%
  distinct()

pdf("Target_pathway_heatmaps.pdf", height = 21)

for (pathway in unique(path_assoc$PATHWAY_NAME)) {
  drugs <- subset(path_assoc, PATHWAY_NAME == pathway)

  associations <- CollecTRI_flt %>%
    inner_join(drugs) %>%
    select(TF, Coefficient, Drug)

  data_wide <- dcast(setDT(associations), Drug ~ TF, value.var="Coefficient")
  data_wide[is.na(data_wide)] <- 0

  data_wide <- data_wide %>% column_to_rownames("Drug")
  h <-Heatmap(t(data_wide), row_title = pathway)
  print(h)
}

dev.off()
