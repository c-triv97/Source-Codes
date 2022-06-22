library(httr)
library(rentrez)
library(jsonlite)
library(rentrez)
library(XML)
library(stringr)
library(Biostrings)
library(ggplot2)
library(tidyverse)
source("https://raw.githubusercontent.com/c-triv97/Source-Codes/main/SeqFindeR.R")

spp1 <- find_gene("Spp1")

cds <- lapply(spp1$Transcript$id, function(x){
  get_sequence(seq_ids = x, type = "cds")
})

cdna <- lapply(spp1$Transcript$id, function(x){
  get_sequence(seq_ids = x, type = "cdna")
})  

genomic <- get_sequence(spp1$id, type = "genomic")

exons <- lapply(spp1$Transcript$Exon[[1]]$id, function(x){
  get_sequence(seq_ids = x, type = "cdna")
})

exons.dat <- as.data.frame(do.call(rbind, exons))

start.end <- list()
for (n in 1:length(exons)){
  start.end[[n]] = str_locate(genomic, exons[[n]])
}

exons.dat <- cbind(exons.dat, as.data.frame(do.call(rbind, start.end)))%>%
  mutate(id = spp1$Transcript$Exon[[1]]$id)
  
primer.in.seqs <- lapply(cdna, function(x){
  list = findORFsinSeq(x)
  
  plot = plotPotentialStartsAndStops(x)
  print(plot)
  return(list)
})

l = list()
for (n in 1:length(spp1$Transcript$id)){
  l[[n]] = str_locate(cdna[[n]], cds[[n]])
}

primer.fwd1 <- ""

primer.rev1 <- ""

primer.fwd2 <-""

primer.rev2 <- ""

p1 <- ggplot()+
  geom_rect(aes(xmin=0, xmax=nchar(genomic), ymin = 0.3, ymax = 0.35), fill = "grey")+
  geom_rect(data = exons.dat, aes(xmin = start, xmax = end, ymin = 0.5, ymax = 0.55), fill = "red")+
  geom_text(data = exons.dat, aes(x = start, y = 0.55, label = c(1:7)))+
  geom_label(aes(x= 500, y= 0.35, label = spp1$id))+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(labels = unlist(strsplit(genomic, split = "")), breaks = c(1:nchar(genomic)))+
  theme_void()

