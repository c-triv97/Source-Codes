library(VariantAnnotation)
library(ggplot2)
library(stringr)
library(tidyverse)

extract_basic_info <- function(VCF_obj, sample_name) {
    
    genotype <- geno(VCF_obj)$GT

    tbl_dat <- as.data.frame(table(genotype)) %>%
           mutate(genotype = str_replace(genotype,
                  pattern = "/",
                  replacement = "|")) %>%
           group_by(genotype) %>%
           summarise(n = sum(Freq))

    p1 <- ggplot(tbl_dat,aes(x = genotype, y = n, fill = genotype))+
        geom_bar(stat='Identity') +
        scale_y_continuous(trans = "sqrt") + 
        labs(x="", y="Counts", fill="", subtitle = paste0(sample_name, " counts by genotype"),
        caption = paste0("total count: ", sum(tbl_dat$n))) +
        theme_light(base_size = 22)

    depth <- geno(VCF_obj)$DP

    print(summary(as.vector(depth)))

    depth <- as.data.frame(depth)
    colnames(depth) <- c("Sample")

    p2 <- ggplot(depth,aes(x=Sample)) +
        geom_histogram(fill = "blue", alpha = 0.6) +
        labs(x = "", y = "Counts", caption = paste0(sample_name, ": approximate read depth")) +
        scale_x_log10() +
        theme_light(base_size = 22)

    GQ <- geno(VCF_obj)$GQ

    print(summary(as.vector(GQ)))

    GQ <- as.data.frame(GQ)

    colnames(GQ) <- c("Sample")

    p3 <- ggplot(as.data.frame(GQ), aes(x=Sample)) +
        geom_histogram(fill = "blue", alpha = 0.6) +
        labs(x = "", y = "Counts", caption = paste0(sample_name, "phred score (genotype calling quality)")) +
        theme_light(base_size = 22)
    
    return(list(p1, p2, p3))
}

extracting_het <- function(VCF_obj){
    dat <- as.data.frame(genotype) %>%
                  mutate(genotype = str_replace(genotype,
                  pattern = "/",
                  replacement = "|"))

    var_1 <- rownames(dat)[
    dat$genotype=="1|1"]

    variant_info <- rowRanges(VCF_obj) #extract GRanges obj from the VCF object 
# extract info 
    varTab1 <- data.frame(variant=names(variant_info)[names(variant_info) %in% var_1],
                        chr=as.vector(seqnames(variant_info)[names(variant_info) %in% var_1]),
                        start=start(variant_info)[names(variant_info) %in% var_1],
                        end=end(variant_info)[names(variant_info) %in% var_1],
                        stringsAsFactors = FALSE)
# ref allele retrieved from ref(vcf)
    ref_base <- ref(VCF_obj)[rownames(VCF_obj) %in% var_1]
    ref_base[1:2]
    varTab1$refBase <- as.character(ref_base)
# alt alleles are retrieved from alt(vcf)
    alt_base <- alt(VCF_obj)[rownames(VCF_obj) %in% var_1]
    alt_base <- lapply(alt_base, function(x){ as.character(x)} ) #this takes time 
    varTab1$altBase <- unlist(alt_base)
            

}