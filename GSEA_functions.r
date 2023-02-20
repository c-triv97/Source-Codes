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
                theme_pubr(base_size = 16, legend = "right") +
                scale_size_continuous(range = c(10, 20))+
                scale_color_viridis_c()

        print(plot)
        ggsave(plot, filename = paste0("plots/",
                                       comparison,
                                       i,
                                       ".pdf", 
                                       sep = ""),
               width = 16)

        goResults[[i]] = results
    }

    return(goResults)
}
