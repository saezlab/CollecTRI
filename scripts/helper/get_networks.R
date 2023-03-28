filter_for_TFs <- function(df){
  # Read TFs
  fname=file.path('data', 'TFcategory', 'TF_category.csv')
  TFcategory <- read.csv2(fname)
  colnames(TFcategory) <- c("source", "TF.category")
  df <- dplyr::filter(df, source %in% TFcategory$source)
  df <- dplyr::left_join(df, TFcategory)
}


get_regnetwork <- function(path, confs=c('High', 'Medium', 'Low'), filterTFs = T){
  temp <- tempfile()
  download.file("https://regnetworkweb.org/download/human.zip",temp, 'curl', extra='--insecure')
  df <- read.table(unz(temp, "human.source"))
  unlink(temp)

  colnames(df) <- c('tf', 'weight', 'target', 'likelihood')
  df['weight'] <- 1
  df <- df[, -which(names(df) %in% c("likelihood"))]

  # Filter by Lambert
  #df <- filter_by_Lambert(df, path)

  # Remove miRNA targets
  df <- dplyr::filter(df, !grepl('hsa-', target))
  df <- dplyr::filter(df, !grepl('MIR', target))

  # Add confidence
  df['confidence'] <- "A"
  colnames(df) <- c("source", "weight", "target", "confidence")
  df <- dplyr::distinct(df, source, target, weight, .keep_all = TRUE)

  filter_idx <- ""

  if (filterTFs) {
    df <- filter_for_TFs(df)
    filter_idx <- "filtered_"
  }

  # Write
  readr::write_csv(df, file.path(path, paste0(filter_idx, 'regnetwork.csv')))
}


unify_chea3_names <- function(tf_list) {
  # Insert '-' for NKX entries
  nkx_last_number <- stringr::str_sub(stringr::str_subset(tf_list, "NKX"), -1)
  nkx_rest <- stringr::str_sub(stringr::str_subset(tf_list, "NKX"), 1,-2)
  nkx_corrected <- paste(nkx_rest, nkx_last_number, sep="-")
  # Correct NKX with '-' and mutate aliases
  tf_list <- stringr::str_subset(tf_list, "NKX.*", negate = TRUE)
  tf_list <- append(tf_list, nkx_corrected)
  tf_list <- stringr::str_replace(tf_list, "ZNF875", "HKR1")
  tf_list <- stringr::str_replace(tf_list, "TBXT", "T")
  tf_list <- stringr::str_replace(tf_list, "CBLL2","ZNF645")
  tf_list <- stringr::str_replace(tf_list, "ZNF788P", "ZNF788")
  tf_list <- stringr::str_replace(tf_list, "ZUP1", "ZUFSP")
  tf_list
}


get_chea3 <- function(path, filterTFs = T){
  # Download chea3
  url <- 'https://maayanlab.cloud/chea3/assets/tflibs/{name}.gmt'
  names <- c(
    'ARCHS4_Coexpression',
    'ENCODE_ChIP-seq',
    'Enrichr_Queries',
    'GTEx_Coexpression',
    'Literature_ChIP-seq',
    'ReMap_ChIP-seq'
  )
  for (name in names) {
    get <- stringr::str_glue(url)
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    download.file(get, fname)
    Sys.sleep(1)
  }

  # Extract and merge
  df <- lapply(names, function(name){
    # Read file
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    f <- file(fname,open="r")
    lines <-readLines(f)
    close(f)

    dfs <- lapply(lines, function(line){
      # Remove tabs
      line <- stringr::str_split(line, '\t')[[1]]

      # Select first element as TF
      tf <- stringr::str_split(line[1], '_')[[1]][1]

      # The rest are targets
      target <- line[2:length(line)]

      # Repeat elements
      tf <- rep(tf, length(target))
      confidence <- rep(name, length(target))

      # Save as df
      data.frame('tf'=tf, 'confidence'=confidence, 'target'=target)
    })
    dfs <- do.call(rbind, dfs)
  })

  # Format
  df <- do.call(rbind, df)
  df['mor'] <- 1
  df['likelihood'] <- 1

  # Change mislabeled TFs
  df[['tf']] <- unify_chea3_names(df[['tf']])

  # Filter by Lambert
  #df <- filter_by_Lambert(df, path)

  # Remove duplicates
  df <- dplyr::distinct(df, tf, target, confidence, .keep_all = TRUE)
  df <- dplyr::mutate(df, weight = mor*likelihood)

  df <- df[, -which(names(df) %in% c("likelihood", "mor"))]

  colnames(df) <- c("source", "confidence", "target", "weight")

  filter_idx <- ""

  if (filterTFs) {
    df <- filter_for_TFs(df)
    filter_idx <- "filtered_"
  }

  # Save
  readr::write_csv(df, file.path(path, paste0(filter_idx, 'chea3.csv')))

  # Remove tmp files
  for (name in names) {
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    file.remove(fname)
  }
}

get_pathwayCommons <- function(path, filterTFs = T){
  url <- 'https://www.pathwaycommons.org/archives/PC2/v11/PathwayCommons11.All.hgnc.sif.gz'
  fname <- file.path(path, 'raw', 'pathwaycommons.sif.gz')
  download.file(url, fname)
  Sys.sleep(1)

  pathwayCommons <- read.table(file.path(path, 'raw', 'pathwaycommons.sif.gz'), header = F)

  pathwayCommons <- dplyr::rename(pathwayCommons, "weight" = V2, "source" = V1, "target" = V3)
  pathwayCommons <- dplyr::filter(pathwayCommons, weight == "controls-expression-of")
  pathwayCommons <- dplyr::mutate(pathwayCommons, weight = 1)
  pathwayCommons <- cbind(pathwayCommons, data.frame(mor = 1))

  filter_idx <- ""

  if (filterTFs) {
    pathwayCommons <- filter_for_TFs(pathwayCommons)
    filter_idx <- "filtered_"
  }

  # Write
  readr::write_csv(pathwayCommons, file.path(path, paste0(filter_idx, 'pathwayCommons.csv')))
}


get_dorothea <- function(path, confs=c('A', 'B', 'C'), filterTFs = T){
  df <- decoupleR::get_dorothea(levels = confs)
  df <- cbind(df, weight = 1)

  filter_idx <- ""

  if (filterTFs) {
    df <- filter_for_TFs(df)
    filter_idx <- "filtered_"
  }

  # Write
  readr::write_csv(df, file.path(path, paste0(filter_idx, 'dorothea_ABC.csv')))
}


get_data <- function(path, filterTFs = T){
  # Create dir
  dir.create(path, showWarnings = F, recursive = T)

  # Get networks
  get_regnetwork(path, filterTFs = filterTFs)
  get_chea3(path, filterTFs = filterTFs)
  get_pathwayCommons(path, filterTFs = filterTFs)
  get_dorothea(path, filterTFs = filterTFs)
}

# Run
# get_data('data/networks')
# get_data('data/networks', filterTFs = F)
