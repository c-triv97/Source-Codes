library(ggplot2)
library(ggpubr)
library(tidyverse)
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

myVolcano <- function(dataframe, title, logcutoff) {

        dat = dataframe %>%
                dplyr::mutate(sig = ifelse(padj < .05 & log2FoldChange > logcutoff, "UP", 
        ifelse(padj < .05 & log2FoldChange < -logcutoff, "DOWN", "NS")),
        padj = dplyr::case_when(padj == 0 ~ 1e-307, TRUE ~ padj),
        logpadj = -log10(padj))

        top_up = subset(dat, sig == "UP") %>%
                dplyr::arrange((abs(padj))) %>%
                dplyr::slice(1:5)
        
        top_down = subset(dat, sig == "DOWN") %>%
                dplyr::arrange((abs(padj))) %>%
                dplyr::slice(1:5)
        
        colors = c("NS" = "#636363", "UP" = "#66c2a5", "DOWN" = "#fc8d62")

        plot = ggplot(dat,
                      aes(x = log2FoldChange, y = logpadj, color = sig)) +
            geom_point(alpha = 0.6, size = 3) +
            geom_label_repel(data = top_up,
                             aes(label = external_gene_name),
                             show.legend = FALSE,
                             max.overlaps = Inf) +
            geom_label_repel(data = top_down,
                             aes(label = external_gene_name),
                             show.legend = FALSE,
                             max.overlaps = Inf) +
            theme_pubr(base_size = 16, legend = "right") +
            scale_color_manual(values = colors) +
            labs(subtitle = title,
                 y = "padj (-log10)",
                 x = "Fold Change (log2)") +
            scale_x_continuous(limits = c(-10, +15))
        return(plot)
}

