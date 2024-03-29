## Load libraries and data
library(tidyverse)
library(decoupleR)
library(magrittr)

## Load files
dorothea_ABC <- read.csv("data/networks/dorothea_ABC.csv")
CollecTRI <- read.csv("output/CollecTRI/CollecTRI_GRN.csv") %>% rename(mor = weight)

# In this use case we use data from CPTAC and three cancer types:
# UCEC: Uterine Corpus Endometrial Carcinoma
# LUAD: Lung Adenocarcinoma
# CCRCC: Clear Cell Renal Cell Carcinoma
download.file("https://zenodo.org/record/7773985/files/ucec_counts_tvalues.csv?download=1", file.path("data", "CPTAC_DEGs", "ucec_counts_tvalues.csv"))
download.file("https://zenodo.org/record/7773985/files/luad_counts_tvalues.csv?download=1", file.path("data", "CPTAC_DEGs", "luad_counts_tvalues.csv"))
download.file("https://zenodo.org/record/7773985/files/ccrcc_counts_tvalues.csv?download=1", file.path("data", "CPTAC_DEGs", "ccrcc_counts_tvalues.csv"))

# Read data in a list of dataframes
file_names = list.files(path ="data/CPTAC_DEGs", pattern="*.csv", full.names = T)
file_list = lapply(file_names, read.csv)
file_list = setNames(file_list, gsub("data/CPTAC_DEGs/|.csv|_counts_tvalues","",file_names))

# Format input for decoupleR
decoupler_inputs <- lapply(file_list, function(x) as.data.frame(x) %>%
                             set_rownames(.$ID) %>% dplyr::select(NATvsTUM_t) %>% filter(!is.na(NATvsTUM_t)))

res_decoupler_CollecTRI <- lapply(decoupler_inputs, function(x) run_ulm(as.matrix(x),network = CollecTRI,.source='source',
                                                                         .target='target', minsize = 5))

dir.create(file.path("output", "case_study", "cptac"), showWarnings = FALSE)
save(res_decoupler_CollecTRI, file = "output/case_study/cptac/res_decoupler_CollecTRI.RData")
#load("output/case_study/cptac/res_decoupler_CollecTRI.RData")
res_decoupler_CollecTRI_flt <- lapply(res_decoupler_CollecTRI, function(x) x %>% dplyr::filter(statistic == "ulm" & p_value < 0.05 & condition == "NATvsTUM_t"))

res_decoupler_CollecTRI_flt <- Map(cbind, res_decoupler_CollecTRI_flt, Cancer_type = names(res_decoupler_CollecTRI_flt)) %>%
  purrr::reduce(rbind) %>%
  select(source,score,p_value,Cancer_type) %>%
  mutate(Network = "CollecTRI") %>%
  mutate(score = round(score, digits = 2))

## Quantify the number of TFs per cancer type
res_decoupler_CollecTRI_flt %>%
  group_by(Cancer_type) %>%
  summarise(n())

res_decoupler_CollecTRI_flt <- res_decoupler_CollecTRI_flt %>%
  split(~Cancer_type)

res_decoupler_CollecTRI_flt <- lapply(res_decoupler_CollecTRI_flt, function(x) x %>% arrange(score))

res_decoupler_df <- res_decoupler_CollecTRI_flt  %>%
  purrr::reduce(rbind) %>%
  select(-Network)

write.table(res_decoupler_df, "output/case_study/cptac/SuppFile2_TFactivities_cptac.csv",sep = ",", row.names = F, quote = F)

## Define function for plotting TF activities as a dotplot
plot_TFs <- function(df) {
  ggplot(data=df,aes( y = reorder(source,score), x=Network, color = score, size = -log10(p_value))) +
    geom_point() +
    theme_minimal() +
    scale_color_gradientn(colors = c("blue","white","red"), limits = c(-5,5), breaks = c(-5,0,5), guide = "none") +
    ylab("") +
    xlab("")+
    #   facet_wrap(~Cancer_type, scales = "free", nrow = 3) +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1),legend.position = "none")
}

## Individual plots per cancer type
plot_TFs(res_decoupler_CollecTRI_flt[["ucec"]])
plot_TFs(res_decoupler_CollecTRI_flt[["ccrcc"]])
plot_TFs(res_decoupler_CollecTRI_flt[["luad"]])


## Print predicted TFs that are not part of the DoRothEA ABC resource
setdiff(res_decoupler_CollecTRI_flt[["ucec"]]$source, dorothea_ABC$source)
setdiff(res_decoupler_CollecTRI_flt[["ccrcc"]]$source, dorothea_ABC$source)
setdiff(res_decoupler_CollecTRI_flt[["luad"]]$source, dorothea_ABC$source)
