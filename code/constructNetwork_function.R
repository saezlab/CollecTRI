#' This function takes the curated table from Eirini and transforms it into a GRN used for decoupler
#' mor and likelihood can be added to the network
#'
#' @param GRNcuration Table with curation information for interaction
#' @param curated_counts whether curated counts should be transformed. Can be "standard" (uses the counts from the table), "binary" (transforms all evidence to 1 indepent of the number of counts) or "scaled" (scales the counts within each resource from 0 to 1)
#' @param mor whether mode of reaction should be added (TRUE or FALSE), if mor should be added rm_bidirectional_rows should be set to TRUE
#' @param mor_filter can be a numeric value between 0-1 so that only for a given ration a negative sign will be assigned
#' @param weight how the interactions should be weighted. Can be "none", "evidence","evidence_sign" or "discrepancy_sign". "Evidence" is supposed to be used if no mode or reaction is FAKSE. It is just based on the amount of evidence ignoring their direction.
#' @param rm_bidirectional_rows removes rows were no sign can be assigned. Based on missing information in TotPositive + TotNeg or same number of information for both. Should be set to TRUE if mor is TRUE or evidence_sign is used
#' @return Constructed network for decoupler.
#'
#' @import tidyverse
#'
#' @export

constructNetwork <- function(GRNcuration, curated_counts = "standard", mor = TRUE, mor_filter = NULL, weight = "none", rm_bidirectional_rows = TRUE){
  # change number of curated information to binaries or scaled values (0-1)
  if (curated_counts == 'binary') {
    GRNcuration <- select(GRNcuration, -totNeg,-totPositive,-totUnknown)
    GRNcuration[GRNcuration > 0] <- 1
    GRNcuration <- GRNcuration %>% add_column(totNeg = rowSums(GRNcuration[,grepl( "Negative", colnames(GRNcuration), fixed = TRUE)])) %>%
      add_column(totPositive = rowSums(GRNcuration[,grepl( "Positive", colnames(GRNcuration), fixed = TRUE)])) %>%
      add_column(totUnknown = rowSums(GRNcuration[,grepl( "Unknown", colnames(GRNcuration), fixed = TRUE)]))
  } else if (curated_counts == "scaled"){
    GRNcuration <- select(GRNcuration, -totNeg,-totPositive,-totUnknown)
    GRNcuration <- apply(GRNcuration,2,function(x){x/max(x)}) %>% replace_na(0) %>% as.data.frame()
    GRNcuration <- GRNcuration %>% add_column(totNeg = rowSums(GRNcuration[,grepl( "Negative", colnames(GRNcuration), fixed = TRUE)])) %>%
      add_column(totPositive = rowSums(GRNcuration[,grepl( "Positive", colnames(GRNcuration), fixed = TRUE)])) %>%
      add_column(totUnknown = rowSums(GRNcuration[,grepl( "Unknown", colnames(GRNcuration), fixed = TRUE)]))
  }

  if (rm_bidirectional_rows){
    # remove rows with no information about direction at all
    GRNcuration <- GRNcuration %>% dplyr::filter(rowSums(GRNcuration[c("totNeg", "totPositive")]) > 0)

    # remove rows with same information for positive and negative
    nrow(GRNcuration %>% dplyr::filter(GRNcuration$totNeg == GRNcuration$totPositive)) # number of connections lost
    GRNcuration <- GRNcuration %>% dplyr::filter(!GRNcuration$totNeg == GRNcuration$totPositive)
  }


  ## GRN construction
  # change to network structure for decoupler
  GRN <- rownames(GRNcuration) %>%
    str_split(':', simplify = TRUE) %>%
    as.data.frame()

  GRN <- GRN %>% rename(source = V1, target = V2)

  # Add generic weight, likelihood and confidence (1,1,A)
  GRN <- GRN %>% add_column(likelihood = 1, mor = 1, confidence = "A")


  ## Change mode of reaction
  if (mor) {
    morEvidence <- GRNcuration %>% select(c("totNeg","totPositive")) %>% max.col()
    if (is.numeric(mor_filter)){
      ratio <- GRNcuration %>% rowwise() %>% mutate(weight = max(c(totNeg, totPositive))/sum(c(totNeg, totPositive, totUnknown))) %>% pull(weight)
      mor <- morEvidence %>% replace(morEvidence == 1, -1) %>% replace(morEvidence == 2, 1)
      morRatio <- data.frame(mor = mor,
                            ratio = ratio)
      mor <- morRatio %>% mutate(mor = replace(mor, morRatio$ratio < mor_filter, 1)) %>% pull(mor)
      GRN <- GRN %>% select(-mor) %>% add_column(mor)
    } else {
      mor <- morEvidence %>% replace(morEvidence == 1, -1) %>% replace(morEvidence == 2, 1)
      GRN <- GRN %>% select(-mor) %>% add_column(mor)
    }
  }


  ## Change weights
  if (weight == "evidence") {
    # weight based on amount of supporting information
    likelihood <- GRNcuration %>% rowwise() %>% mutate(weight = sum(c(totNeg, totPositive, totUnknown))/sum(!grepl("tot", colnames(GRNcuration)))) %>% pull(weight)
    GRN <- GRN %>% select(-likelihood) %>% add_column(likelihood)
  } else if (weight == "evidence_sign"){
    # weight based on amount of information for sign
    likelihood <- GRNcuration %>% rowwise() %>% mutate(weight = max(c(totNeg, totPositive, totUnknown))/sum(c(totNeg, totPositive, totUnknown))) %>% pull(weight)
    GRN <- GRN %>% select(-likelihood) %>% add_column(likelihood)
  } else if (weight == "discrepancy_sign"){
    # weight based on discrepancies in sign
    likelihood <- GRNcuration %>% rowwise() %>% mutate(weight = abs(totNeg*0+totPositive*1+totNeg*-1)/sum(c(totNeg, totPositive, totUnknown))) %>% pull(weight)
    GRN <- GRN %>% select(-likelihood) %>% add_column(likelihood)
  } else if (weight == "evidence_tot"){
    max_evidence <- sum(max(GRNcuration$totNeg) + max(GRNcuration$totPositive) + max(GRNcuration$totUnknown))
    likelihood <- GRNcuration %>% rowwise() %>% mutate(weight = sum(c(totNeg, totPositive, totUnknown))/max_evidence) %>% pull(weight)
    GRN <- GRN %>% select(-likelihood) %>% add_column(likelihood)
  }

  return(GRN)

}
