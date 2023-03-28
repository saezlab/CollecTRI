get_Lambert <- function(path){
  # Download Lambert list of TFs
  url <- 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867418301065-mmc2.xlsx'
  fname <- file.path(path, 'lambert.csv')
  download.file(url, fname)
  Sys.sleep(1)

  # Select the ones annotated as TF
  df <- readxl::read_xlsx(fname, sheet=2, skip=1)
  df <- dplyr::select(df, Name, `...4`)
  df <- dplyr::filter(df, `...4`=='Yes')
  df <- dplyr::select(df, Name)

  # Write
  write.csv2(df, file.path(path, 'lambert.csv'), row.names=F)
}

get_Lovering <- function(path){
  # Download Lovering list of TFs
  url <- 'https://ars.els-cdn.com/content/image/1-s2.0-S1874939921000833-mmc1.xlsx'
  fname <- file.path(path, 'lovering.csv')
  download.file(url, fname)
  Sys.sleep(1)

  # Select the ones annotated as TF
  df <- readxl::read_xlsx(fname, sheet=2, col_names = T)
  df <- dplyr::select(df, `HGNC approved gene symbol`)
  df <- dplyr::rename(df, "Name" = `HGNC approved gene symbol`)

  # Write
  write.csv2(df, file.path(path, 'lovering.csv'), row.names=F)
}

get_TFclass <- function(path){
  # Download TFclass list of TFs
  url <- 'http://tfclass.bioinf.med.uni-goettingen.de/suppl/tfclass.ttl.gz'
  fname <- file.path(path, 'raw', 'TFclass.ttl.gz')
  download.file(url, fname)
  Sys.sleep(1)

  # Read the file line by line and store each line in a character vector
  ttf_file <- gzfile(file.path(path, 'raw', 'TFclass.ttl.gz'), "r")
  lines <- readLines(ttf_file)
  close(ttf_file)

  # Extract human TFs
  # Define the regular expression
  regex <- "Homo_sapiens_\\w+" #select human TFs
  # Use grep to find all matching words in the lines vector
  matching_words <- grep(regex, lines, value = TRUE)

  # Extract mouse TFs
  # Define the regular expression
  regex_mouse <- "Mus_musculus_\\w+" #select human TFs
  # Use grep to find all matching words in the lines vector
  matching_words_mouse <- grep(regex_mouse, lines, value = TRUE)

  # Extract rat TFs
  # Define the regular expression
  regex_rat <- "Rattus_norvegicus_\\w+" #select human TFs
  # Use grep to find all matching words in the lines vector
  matching_words_rat <- grep(regex_rat, lines, value = TRUE)

  # select TF symbols
  tfs <- stringr::str_match(matching_words, "#Homo_sapiens_([[:alnum:]_-]+)")[, 2]
  tfs <- unique(tfs[!is.na(tfs) ])

  tfs_mouse <- stringr::str_match(matching_words_mouse, "#Mus_musculus_([[:alnum:]_-]+)")[, 2]
  tfs_mouse <- unique(tfs_mouse[!is.na(tfs_mouse) ])

  tfs_rat <- stringr::str_match(matching_words_rat, "Rattus_norvegicus_([[:alnum:]_-]+)")[, 2]
  tfs_rat <- unique(tfs_rat[!is.na(tfs_rat) ])

  df <- data.frame(Name = union(union(tfs, tfs_mouse), tfs_rat))

  # Write
  write.csv2(df, file.path(path, 'TFclass.csv'), row.names=F)
}

get_dbTFs <- function(path){
  # Load TFs classified as coTFs according to GO (GO:0003700)
  fname <- file.path(path, 'raw', 'QuickGO-annotations-1677666121160-20230301.tsv')

  # Select the ones annotated as TF
  df <- utils::read.table(fname, sep = "\t", header = T)
  tfs <- dplyr::pull(df, SYMBOL)
  df <- data.frame("Name" = unique(tfs))

  # Write
  write.csv2(df, file.path(path, 'dbTF_quickGO.csv'), row.names=F)
}

get_coTFs <- function(path){
  # Load human TFs classified as coTFs according to GO (GO:0003712)
  fname <- file.path(path, 'raw', 'QuickGO-annotations-1676559707182-20230216.tsv')

  # Select the ones annotated as human TF
  df <- utils::read.table(fname, sep = "\t", header = T)
  tfs <- dplyr::pull(df, SYMBOL)

  # mouse coTFs according to GO (GO:0003712)
  fname_mouse <- file.path(path, 'raw', 'QuickGO-annotations-1678184433156-20230307.tsv')
  df_mouse <- utils::read.table(fname_mouse, sep = "\t", header = T)
  tfs_mouse <- dplyr::select(df_mouse, SYMBOL)

  #conversion to human
  human_gene_names <- OmnipathR::homologene_download(target = 9606, source = 10090, id_type = "genesymbol")
  human_gene_names <- dplyr::rename(human_gene_names, "SYMBOL" = genesymbol_source)
  tfs_mouse <- dplyr::left_join(tfs_mouse, human_gene_names, by = "SYMBOL", multiple = "all")
  tfs_mouse <- dplyr::pull(tfs_mouse, genesymbol_target)
  tfs_mouse <- tfs_mouse[!is.na(tfs_mouse)]

  # rat coTFs according to GO (GO:0003712)
  fname_rat <- file.path(path, 'raw', 'QuickGO-annotations-1678184598281-20230307.tsv')
  df_rat <- utils::read.table(fname_rat, sep = "\t", header = T)
  tfs_rat <- dplyr::select(df_rat, SYMBOL)

  #conversion to human
  human_gene_names <- OmnipathR::homologene_download(target = 9606, source = 10116, id_type = "genesymbol")
  human_gene_names <- dplyr::rename(human_gene_names, "SYMBOL" = genesymbol_source)
  tfs_rat <- dplyr::left_join(tfs_rat, human_gene_names, by = "SYMBOL", multiple = "all")
  tfs_rat <- dplyr::pull(tfs_rat, genesymbol_target)
  tfs_rat <- tfs_rat[!is.na(tfs_rat)]

  df <- data.frame("Name" = union(union(tfs, tfs_mouse), tfs_rat))

  # Write
  write.csv2(df, file.path(path, 'coTF_quickGO.csv'), row.names=F)
}

get_GTFs <- function(path){
  # Load TFs classified as GTFs according to GO (GO:0140223)
  fname <- file.path(path, 'raw', 'QuickGO-annotations-1676559912711-20230216.tsv')

  # Select the ones annotated as TF
  df <- utils::read.table(fname, sep = "\t", header = T)
  tfs <- dplyr::pull(df, SYMBOL)
  df <- data.frame("Name" = unique(tfs))

  # Write
  write.csv2(df, file.path(path, 'GTF_quickGO.csv'), row.names=F)
}

get_TFcategory <- function(path){
  # Create dir
  dir.create(path, showWarnings = F, recursive = T)

  # Get networks
  get_Lambert(path)
  get_Lovering(path)
  get_TFclass(path)
  get_dbTFs(path) #manually downloaded
  get_coTFs(path) #manually downloaded
  get_GTFs(path) #manually downloaded
}


# Run
# get_TFcategory('data/TFcategory')

generate_TFlist <- function(path){
  tf_files <- list.files(path, pattern = ".csv")
  tf_lists <- purrr::map(file.path(path, tf_files), read.csv)
  names(tf_lists) <- stringr::str_remove(tf_files, ".csv")

  dbTFs <- unique(c(tf_lists$lambert$Name,
                    tf_lists$lovering$Name,
                    tf_lists$TFclass$Name,
                    tf_lists$dbTF_quickGO$Name))

  TFs <- data.frame(TF = dbTFs,
                    class = "dbTF")

  coTFs <- tf_lists$coTF_quickGO$Name[!tf_lists$coTF_quickGO$Name %in% TFs$TF]

  TFs <- rbind(TFs,
               data.frame(TF = coTFs,
                          class = "coTF"))

  GTFs <- tf_lists$GTF_quickGO$Name[!tf_lists$GTF_quickGO$Name %in% TFs$TF]

  TFs <- rbind(TFs,
               data.frame(TF = GTFs,
                          class = "GTF"))

  # Write
  write.csv2(TFs, file.path(path, 'TF_category.csv'), row.names=F)
}

# Run
# generate_TFlist('data/TFcategory')
