
deseq_transform <- function(deseq2, p_co, lfc_co) {
    
    deseq2 %<>%
        as.data.frame() %>%
        rownames_to_column(var = "symbol") %>%
        as_tibble() %>%
        filter(!is.na(padj)) %>%
        mutate(significant = ifelse(padj <= p_co & abs(log2FoldChange) >= lfc_co, 
                                    TRUE, FALSE))
    return(deseq2)
}