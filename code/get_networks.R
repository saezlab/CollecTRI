# from Pau

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


filter_by_Lambert <- function(df, path){
  # Read TFs
  fname=file.path(path, 'lambert.csv')
  lambert <- read.csv2(fname)$Name
  dplyr::filter(df, tf %in% lambert)
}


get_regnetwork <- function(path, confs=c('High', 'Medium', 'Low')){
  temp <- tempfile()
  download.file("https://regnetworkweb.org/download/human.zip",temp)
  data <- read.table(unz(temp, "human.source"))
  unlink(temp)

  colnames(data) <- c('tf', 'weight', 'target', 'likelihood')
  data['weight'] <- 1
  data <- data[, -which(names(data) %in% c("likelihood"))]

  # Filter by Lambert
  #df <- filter_by_Lambert(data, path)

  # Remove miRNA targets
  df <- dplyr::filter(df, !grepl('hsa-', target))
  df <- dplyr::filter(df, !grepl('MIR', target))

  # Add confidence
  df['confidence'] <- "A"
  colnames(df) <- c("source", "weight", "target", "confidence")
  df <- dplyr::distinct(df, source, target, weight, .keep_all = TRUE)

  # Write
  readr::write_csv(df, file.path(path, 'regnetwork.csv'))
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


get_chea3 <- function(path){
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

  # Save
  readr::write_csv(df, file.path(path, 'chea3.csv'))

  # Remove tmp files
  for (name in names) {
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    file.remove(fname)
  }
}



get_data <- function(path){
  # Create dir
  dir.create(path, showWarnings = F, recursive = T)

  # Get networks
  get_Lambert(path)
  get_regnetwork(path)
  get_chea3(path)
}


# Run
get_data('data/raw')
