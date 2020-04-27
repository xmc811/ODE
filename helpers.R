
deseq_format_test <- function(rds) {
    
    if ((class(rds[[1]]) == "DESeqDataSet") & (class(rds[[2]]) == "DESeqResults")) {
        
        return(TRUE)
    }
}

num_colors <- function(pal) {
    
    a <- rownames(brewer.pal.info) == pal
    return(brewer.pal.info$maxcolors[a])
}


deseq_transform <- function(res, p_co, lfc_co) {
    
    res %<>%
        as.data.frame() %>%
        rownames_to_column(var = "symbol") %>%
        as_tibble() %>%
        filter(!is.na(padj)) %>%
        mutate(significant = ifelse(padj <= p_co & abs(log2FoldChange) >= lfc_co, 
                                    TRUE, FALSE))
    return(res)
}

