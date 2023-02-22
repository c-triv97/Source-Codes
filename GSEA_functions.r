library(tidyverse)
library(fgsea)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(goseq)

goseq_object <- function(data,
                         test_cats = c("GO:BP", "GO:MF", "GO:CC"), 
                         comparison){
    data = data

    sig_genes = data$FDR < 0.01 & !is.na(data$FDR)
    genes = as.integer(sig_genes)
    names(genes) = data$ensembl_gene_id

    pwf = nullp(genes, "mm10", "ensGene", bias.data = data$medianTxLength)
    print(head(pwf))

    goResults = list()

    for (i in 1:length(test_cats)){

        results = goseq(pwf,
                        "mm10",
                        "ensGene",
                        test.cats = test_cats[i])
    
        plot <- results %>%
                top_n(10, wt = -over_represented_pvalue) %>%
                mutate(hitPerc = numDEInCat / numInCat * 100) %>%
                ggplot(aes(x = hitPerc,
                y = term,
                color = over_represented_pvalue,
                size = numDEInCat)) +
                geom_point() +
                expand_limits(x=0) +
                labs(x = "Hits (%)", y = test_cats[i],
                    color = "p value",
                    size = "count", 
                    subtitle = comparison)+
                theme_pubr(base_size = 18, legend = "right") +
                scale_size_continuous(range = c(10, 20))+
                scale_color_viridis_c()

        print(plot)
        ggsave(plot, filename = paste0("plots/",
                                       comparison,
                                       i,
                                       ".pdf", 
                                       sep = ""),
               width = 12)

        goResults[[i]] = results
    }

    return(goResults)
}

fgsea_output <- function(dataframe, pathways, filename){

    gseaDat = dataframe

    ranks = gseaDat$logFC

    names(ranks) = gseaDat$entrezgene_id

    fgseaRes = fgsea::fgsea(pathways,
                            ranks, 
                            minSize = 15,
                            maxSize = 500)

    write.table(fgseaRes, file = paste0(filename, "fgsea.txt"), sep = "\t",
                row.names = FALSE)

    top10up = fgseaRes %>%
                filter(ES > 0) %>%
                top_n(10, wt = -padj)
    
    top10down = fgseaRes %>%
                filter(ES < 0) %>%
                top_n(10, wt = -padj)

    all = bind_rows(top10up, top10down)

    plot = ggplot(all, aes(y = pathway, x = NES, size = size)) +
        geom_point(alpha = 0.6, aes(color = padj)) +
        guides(size = guide_legend(ncol = 1),
            color = "none") +
        theme_pubr(legend = "right") +
        scale_size_continuous(range = c(5, 10)) +
        scale_color_continuous()

    print(plot)

    return(fgseaRes)
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
                   nPerm = 10000,
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
        geom_point(aes(x = -log(p.adjust), y = Description, size = NES,
                        color = comparison),
                alpha = 0.6) +
        theme_bw(base_size = 20) +
        labs(subtitle = "Top 10 most significantly enriched GO terms in each comparison group")
    
    print(last_plot())

    return(dat)

    if (save == TRUE)
    ggsave(last_plot(), filename = "comparison_plot.pdf", height = 12, widht = 8)
}
