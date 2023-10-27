library(ggplot2)
library(ggpubr)
library(tidyverse)
source("https://raw.githubusercontent.com/c-triv97/Source-Codes/main/ThesisTheme.r")
biomart.packages <- c("DESeq2",
                      "biomaRt",
                      "vsn",
                      "rtracklayer",
                      "GenomicFeatures",
                      "org.Rn.eg.db",
                      "ensembldb",
                      "AnnotationHub",
                      "tximport",
                      "rhdf5"
)

pacman::p_load(biomart.packages, character.only = TRUE)

ah <- AnnotationHub()
#EnsDb.rat <- ah[["AH113790"]] # shrsp
EnsDb.rat <- ah[["AH113793"]]
txs <- transcripts(EnsDb.rat, return.type = "DataFrame")
gxs <- genes(EnsDb.rat, return.type = "DataFrame")

tx2gene <- as.data.frame(txs) %>% 
           relocate("tx_id_version", "gene_id") 


pca_full <- function(dds, trans = "rlog", title = "", info, anno, keep = "loadings"){ # info should be df with sample id and condition, anno is gene id matching dds plus common symbol, keep decides keeping loading plot or the PCA output object (for recreating any plot)
    dds = dds 

    colnames(info) = c("sample", "condition")

    gxs = as.data.frame(anno) # anno is a gene name to symbol conversion 
    colnames(gxs) = c("gene_id", "gene_name")

    if (trans == "rlog"){
        fit_dds = rlog(dds, fitType = "local")
    } else {
        fit_dds = vst(dds, fitType = "local")
    }

    pca_matrix = assay(fit_dds)

    pca = prcomp(t(pca_matrix),
                 center = TRUE,
                 scale = TRUE) #pca requires matrix with cols presrenting variables and rows representing samples 

    print(summary(pca))

    eigenvalues = pca$sdev^2

    eigenvalues = tibble(PC = factor(1:length(eigenvalues)), 
                         variance = eigenvalues) %>% 
                         mutate(pct = variance/sum(variance)*100) %>%
                         mutate(pct_cum = cumsum(pct))

    eigenplot = eigenvalues %>%
                ggplot(aes(x = PC))+
                geom_col(aes(y= pct, fill = PC)) +
                geom_line(aes(y = pct_cum, group=1), color = "red")+
                geom_point(aes(y =  pct_cum), color = "red") +
                scale_y_continuous(expand = c(0,0), guide = "prism_minor") +
                guides(fill = "none")+
                theme_thesis + 
                labs(x = "Principal Component",
                     y = "Percent Variance (%)",
                     caption="Red line represents cumulative variance explained", 
                     subtitle = paste0("PCA Eigenvalues ", title, sep = ""))
    
    print(eigenplot)

    pca_scores = pca$x %>%
                 as_tibble(rownames = "Sample") %>%
                 mutate(condition = info$condition, 
                        Sample = info$sample)

    pca_plot = ggplot(pca_scores, aes(PC1, PC2, color = condition)) +
               geom_point(size = 8, alpha = 0.5) +
               xlab(paste0("PC1: ", round(eigenvalues$pct[1], 2),"% variance")) +
               ylab(paste0("PC2: ", round(eigenvalues$pct[2], 2),"% variance")) +
               geom_text_repel(aes(label = Sample), size = 8, show.legend = FALSE) +
               scale_y_continuous(guide = "prism_minor", limits = c(min(pca_scores$PC2), 
                                                                    max(pca_scores$PC2))) + 
               scale_x_continuous(guide = "prism_minor", limits = c(min(pca_scores$PC1), 
                                                                    max(pca_scores$PC1))) + 
               theme_thesis +
               scale_color_brewer(palette = "Dark2") +
               labs(subtitle = paste0("PCA Analysis", title, sep=""))
    
    print(pca_plot)

    # loadings 
    pc_loadings <- pca$rotation %>%
                   as_tibble(rownames = "gene")

    top_genes <- pc_loadings %>% # select only the PCs we are interested in
                 select(gene, PC1, PC2) %>% # convert to a "long" format
                 pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% # for each PC
                 group_by(PC) %>%  # arrange by descending order of loading
                 arrange(desc(abs(loading))) %>% 
                 slice_head(n=10) %>% # take the 10 top rows
                 pull("gene") %>% # pull the gene column as a vector
                 unique() # ensure only unique genes are retained

    annotation = gxs %>% 
                 filter(gene_id %in% top_genes) %>% 
                 mutate(gene_name = case_when(gene_name == "" ~ gene_id, 
                                    .default = gene_name))

    top_loadings <- pc_loadings %>% 
                    filter(gene %in% top_genes) %>% 
                    merge(annotation, by.x = "gene", by.y = "gene_id")
                    
    loading_plot = ggplot(data = top_loadings) +
                    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                                arrow = arrow(length = unit(0.1, "in")),
                                color = "blue", alpha = 0.8, linewidth = 2) +
                    geom_label_repel(aes(x = PC1, y = PC2, label = gene_name), size = 4, max.overlaps = Inf) +
                    theme_thesis + 
                    scale_x_continuous(guide = "prism_minor",
                                       limits = c(min(top_loadings$PC1)-0.01, 
                                                  max(top_loadings$PC1)+0.01))+
                    scale_y_continuous(guide = "prism_minor",
                                       limits = c(min(top_loadings$PC2)-0.01, 
                                                  max(top_loadings$PC2)+0.01))+
                    labs(x = paste0("PC1 (", round(eigenvalues$pct[1],1), "%)", sep=""), 
                        y = paste0("PC2 (", round(eigenvalues$pct[2],1), "%)", sep=""), 
                        subtitle = paste0("PC loading plot", title, sep = ""))
    
    
    print(loading_plot)

    if (keep == "loadings"){
        return(loading_plot)
    } else {
        return(pca)
    }

}
## MA Plots R
myMAplot <- function(data, title) {
        dat = as.data.frame(data)
        
        dat = dat %>% 
        dplyr::mutate(sig = ifelse(padj <.05, "TRUE", "FALSE"))

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

        plot = ggplot(dat, 
                  aes(x = baseMean, y = log2FoldChange, color = sig)) +
            geom_point(alpha = 0.6, size = 3)+
            geom_line(aes(y = 0), color = "black", linewidth = 2) +
            geom_label_repel(data = top_n, 
                             aes(label = external_gene_name),
                             show.legend = FALSE,
                             max.overlaps = Inf,
                             nudge_x = 1,
                             nudge_y = 2*sign(top_n$log2FoldChange))+
            scale_x_log10()+
            theme_pubr(base_size = 16, legend = "none") +
            scale_color_manual(values = colors) +
            labs(subtitle = title,
                 y = "Log Fold Change", 
                 x = "Mean of Normalised Counts")
        return(plot)
}

myVolcano <- function(dataframe, title, logcutoff, n=10) {

        dat = dataframe %>%
                dplyr::mutate(sig = ifelse(padj < .05 & log2FoldChange > logcutoff, "UP", 
        ifelse(padj < .05 & log2FoldChange < -logcutoff, "DOWN", "NS")),
        padj = dplyr::case_when(padj == 0 ~ 1e-307, TRUE ~ padj),
        logpadj = -log10(padj))

        top_up = subset(dat, sig == "UP") %>%
                dplyr::arrange((abs(padj))) %>%
                dplyr::slice(1:n)
        
        top_down = subset(dat, sig == "DOWN") %>%
                dplyr::arrange((abs(padj))) %>%
                dplyr::slice(1:n)
        
        colors = c("NS" = "#636363", "UP" = "#66c2a5", "DOWN" = "#fc8d62")

        plot = ggplot(dat,
                      aes(x = log2FoldChange, y = logpadj, color = sig)) +
            geom_point(alpha = 0.6, size = 4) +
            geom_text_repel(data = top_up,
                             aes(label = external_gene_name),
                             show.legend = FALSE,
                             color = "black", 
                             max.overlaps = 5) +
            geom_text_repel(data = top_down,
                             aes(label = external_gene_name),
                             show.legend = FALSE,
                             color = "black", 
                             max.overlaps = 5) +
            theme_thesis + 
            scale_color_manual(values = colors) +
            labs(subtitle = title,
                 y = "padj (-log10)",
                 x = "Log2 Fold Change") +
            scale_x_continuous(guide = "prism_minor", 
                               expand = c(0,0), 
                               limits = c(min(dat$log2FoldChange)-0.1,
                                          max(dat$log2FoldChange)+0.1)) +
            scale_y_continuous(guide = "prism_minor", 
                               limits = c(min(dat$logpadj)-1,
                                          max(dat$logpadj)+1)) +
            geom_vline(xintercept = logcutoff, 
                    linetype = "dashed", 
                    colour = "grey", 
                    linewidth = 1.5, 
                    alpha = 0.5) +
            geom_vline(xintercept = -logcutoff, 
                    linetype = "dashed", 
                    colour = "grey", 
                    linewidth = 1.5, 
                    alpha = 0.5) +
            geom_hline(yintercept = -log10(0.05), 
                    linetype = "dashed", 
                    linewidth = 1.5,
                    alpha = 0.5,  
                    color = "grey") +
            guides(color = "none")
    
        return(plot)
}

