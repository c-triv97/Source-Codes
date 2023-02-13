## plots for DESeq2 

myMAplot <- function(x, title){
    dat = as.data.frame(x)%>%
            mutate(sig = ifelse(padj <.05, "TRUE", "FALSE"))
    
    top_n = rownames(dat %>%
                     arrange((abs(padj))) %>%
                     dplyr::slice(1:10))
          
    
    top_genes <- getBM(filters = "ensembl_gene_id",
                       attributes = c("external_gene_name", 
                                      "ensembl_gene_id"),
                       values = top_n, mart = ensembl)%>%
                mutate(external_gene_name = case_when(external_gene_name == ""  ~
                       ensembl_gene_id,
                       TRUE ~ external_gene_name))

    top_n <- merge(top_genes, dat, by.x = 2, by.y = 0)

    colors = c("FALSE" = "#636363", "TRUE" = "#1f78b4")

    plot = ggplot(dat, aes(x = baseMean, y = log2FoldChange, color = sig))+
            geom_point(alpha = 0.6, size = 3)+
            geom_line(aes(y = 0), color = "black", linewidth = 2)+ 
            geom_label_repel(data = top_n, aes(label = external_gene_name), show.legend = FALSE)+
            scale_x_log10()+
            theme_pubr(base_size = 16, legend = "none")+
            scale_color_manual(values = colors)+
            labs(subtitle = title, 
                 y = "Log Fold Change", 
                 x = "Mean of Normalised Counts")
    
    return(plot)
}
