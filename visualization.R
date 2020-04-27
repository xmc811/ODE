
deseq_volcano <- function(deseq2, p_co, lfc_co) {
    
    deseq2 %<>%
        deseq_transform(p_co, lfc_co)
    
        ggplot(deseq2) +
        geom_point(aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = significant),
                   alpha = 0.5) +
        scale_color_manual(values = c("#bdbdbd", "#de2d26")) +
        scale_x_continuous(limits = c(max(abs(deseq2$log2FoldChange)) * c(-1, 1))) +
        theme(aspect.ratio = 1)
    
}


deseq_table <- function(deseq2, p_co, lfc_co) {
    
    deseq2 %<>%
        deseq_transform(p_co, lfc_co)
    
    deseq2 %<>%
        filter(significant == TRUE) %>%
        select(symbol:padj) %>%
        arrange(padj, abs(log2FoldChange)) %>%
        datatable() %>%
        formatRound(columns = c(2:5), digits = 3) %>%
        formatSignif(columns = c(6:7), digits = 3)
    
    return(deseq2)
}


