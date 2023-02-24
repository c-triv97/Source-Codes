library(tidyverse)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(GOSemSim)

SemSimData <- function(organism, ont="BP"){

    SemSim = godata(organism, ont = ont)

    return(SemSim)
}

gse_clusterprofiler <- function(data, 
                                organism,
                                keys = "ENSEMBL",
                                padj = "BH",
                                printplot = TRUE){
    if (is.data.frame(data)){
        genes <- data$logFC
        names(genes) <- data$ensembl_gene_id
    } else {
        print(head(data))
    }

    gse = gseGO(geneList = data,
                    ont = "ALL",
                    keyType = keys,
                    minGSSize = 3,
                    maxGSSize = 800,
                    pvalueCutoff = 0.05,
                    verbose = TRUE,
                    OrgDb = organism,
                    pAdjustMethod = padj, 
                    eps = 0)
    
    if (printplot == TRUE){
        plot = dotplot(gse, showCategory = 10, split = ".sign") +
                        facet_grid(~.sign)+
                        theme_bw(base_size = 18)
        print(plot)
    }
    return(gse)
}

kegg_clusterprofiler <- function(data, # input needs to have entrez ids 
                                 keys = "ncbi-geneid",
                                 organism = "mmu",
                                 padj = "BH",
                                 printplot = TRUE){
    if (is.data.frame(data)){
        dat <- data$logFC
        names(dat) <- data$entrezgene_id
    } else {
        print(head(data))
        dat = data
    }

    gse <- gseKEGG(geneList = dat,
                   organism = organism,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = padj,
                   keyType = keys)
    
    if (printplot == TRUE){
        plot = dotplot(gse, showCategory = 10, split = ".sign", orderBy = "p.adjust") +
                        theme_bw(base_size = 18)
        print(plot)
    }

    return(gse)
}

plotting_common <- function(datalist, comparison_groups, order, save = FALSE, n=10){

    GOlist <- lapply(datalist[c(comparison_groups)], function(x){
        go_ids = x$ID %>%
            droplevels()
        return(go_ids)
    })

    inCommon <- Reduce(intersect, GOlist)

    dat <- rbindlist(datalist[c(comparison_groups)], idcol = "comparison") %>%
            dplyr::mutate(InCommon = ifelse(ID %in% inCommon, TRUE, FALSE))

    ggplot(dat %>%
        group_by(comparison) %>%
        arrange(p.adjust) %>%
        slice_head(n = n)) +
        geom_point(aes(x = -log(p.adjust), y = fct_reorder(Description, p.adjust),
                        color = NES, shape = comparison),
                alpha = 0.6, 
                size = 8) +
        theme_bw(base_size = 20) +
        labs(subtitle = "Top 10 most significantly enriched GO terms in each comparison group", 
             y = "Description")
    
    print(last_plot())

    return(dat)

    if (save == TRUE)
    ggsave(last_plot(), filename = "comparison_plot.pdf", height = 12, widht = 8)
}