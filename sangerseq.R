# source code for data import from sanger sequencing .ab1 and .scf files 
library(sangerseqR)
library(tidyverse)

data.import.ab1 <- function(pat= ".ab1", path = getwd(), names = c("none")){
  data.files = list.files(pattern = pat, 
                          recursive = T, 
                          path = path)
  data = lapply(data.files, function(x){
    read.abif(x)
  })
  
  data = lapply(data, function(x){
    sangerseq(x)
  })
  
  if (length(names) > 1) {
    names(data) = names
  } else {
    names(data) = data.files
  }
  
  return(data)
}

data.import.scf <- function(pat = ".scf", path = getwd(), names = c("none")){
  
  data.files = list.files(pattern = pat, 
                          recursive = T, 
                          path = path)
  
  data = lapply(data.files, function(x){
    read.scf(x)
  })
  
  if (length(names) > 1) {
    names(data) = names
  } else {
    names(data) = data.files
  }
  
  return(data)
}

extract.files <- function(sample.pattern = "",
                          data.files, 
                          names = c("none")){
  samples <- grep(sample.pattern, names(data.files))
  
  samples1 <- data.files[samples]
  
  if (length(names) > 1) {
    names(samples1) = names
  } else {
    names(samples1) = names(samples1)
  }
  
  return(samples1)
}


