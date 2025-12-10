library(rentrez)
library(tidyverse)
library(taxize)

# Set API key by logging to NCBI ####
set_entrez_key("5d353264467978b2c10088ca2cd2b09adf08")

# API-key rate limit: 10 requests/sec so delay = 0.1 sec ####
REQUEST_DELAY <- 0.12   # buffer to be safe

# function to turn accession number into title, taxid, organism ####
# ------------------------------------------ #
get_taxid_from_accession <- function(accession) {

  Sys.sleep(REQUEST_DELAY)

  rec <- tryCatch(
    entrez_summary(db = "nucleotide", id = accession),
    error = function(e) return(NULL)
  )
  
  if (is.null(rec)) {
    return(list(
      accession = accession,
      match     = NA,
      taxid     = NA,
      organism  = NA
    ))
  }

  list(
    accession = accession,
    match     = rec$title,
    taxid     = rec$taxid,
    organism  = rec$organism
  )
}

# turn taxid into order, class, family ####
# ------------------------------------------
get_ranks_from_taxid <- function(taxid) {

  Sys.sleep(REQUEST_DELAY)

  cl <- classification(taxid, db = "ncbi")

  cl <- cl[[1]]

  class_val <- cl$name[cl$rank == "class"]
  order_val <- cl$name[cl$rank == "order"]
  family_val <- cl$name[cl$rank == "family"]

    list(
        class = ifelse(length(class_val)==0, NA, class_val),
        order = ifelse(length(order_val)==0, NA, order_val),
        family = ifelse(length(family_val)==0, NA, family_val)
    )
}


# Main function to get taxonomy from accession numbers ####
# ------------------------------------------
get_taxonomy_from_accessions <- function(accessions) {

  message("Step 1/2: Fetching accession summaries...")
  base_info <- map_df(accessions, get_taxid_from_accession)

  message("Step 2/2: Fetching taxonomy lineage...")
  lineage_info <- map_df(base_info$taxid, get_ranks_from_taxid)

  bind_cols(base_info, lineage_info)
}
