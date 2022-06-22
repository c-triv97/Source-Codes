#source codes for downloading sequences, aligning primers and producing some metadata
source("https://raw.githubusercontent.com/c-triv97/Source-Codes/main/ComputR_cDNA2ORF.R")
server <- "http://rest.ensembl.org/" 
con <- "application/json"

find_gene <- function(gene_name, species = "rattus_norvegicus"){
  server = "http://rest.ensembl.org/" 
  con = "application/json"
  ext = paste0("lookup/symbol/", species, "/",
               gene_name, "?expand=1", 
               sep ="")
  
  get_gene = fetch_endpoint(server, ext, con)
  
  return(get_gene)
}

get_sequence <- function(seq_ids, type){ 
  #seq ids can be found in meta_data from find_gene function
  #can only be performed on ids one at a time. id type should match type requested
  server = "http://rest.ensembl.org/" 
  con = "application/json"
  
  ext = if (type == "genomic"){
    paste0("sequence/id/", seq_ids, "?type=genomic")
  }
  else if(type == "cds"){
      paste0("sequence/id/", seq_ids, "?type=cds")
  }
  else if(type == "cdna"){
    paste0("sequence/id/", seq_ids, "?type=cdna")
  }
  else if(type == "protein"){
    paste0("sequence/id/", seq_ids, "?type=protein")
  }
  else {
    print("wrong 'type' inputted")
  }
    
  seq = fetch_endpoint(server, ext, con)
  
  return(seq$seq)
}

library(stringi)

reverseComp <- function(seq){
  return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", seq)))
}


