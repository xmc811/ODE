
deseq_heatmap <- function(dds, var, palette, dir) {
    
    vsd <- vst(dds, blind = FALSE)
    
    sampleDists <- dist(t(assay(vsd)))
    
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- vsd[[var]]
    colnames(sampleDistMatrix) <- vsd[[var]]
    
    num <- num_colors(palette)
    
    colors <- brewer.pal(num,palette)
    
    if (dir) {
        colors <- rev(colors)
    }
    
    col_fun <- colorRamp2(seq(from = 0, 
                              to = max(sampleDistMatrix), 
                              length.out = num), colors)

    Heatmap(sampleDistMatrix, 
            col = col_fun,
            rect_gp = gpar(col = "white", lwd = 2),
            width = unit(20, "cm"), height = unit(20, "cm"))
}

deseq_pca <- function(dds, var, palette, dir) {
    
    vsd <- vst(dds, blind = FALSE)
    
    if (is.numeric(dds@colData[[var]])) {
        
        plotPCA(vst(dds, blind = FALSE), intgroup = var) +
            scale_color_distiller(palette = palette, direction = dir) +
            labs(color = var) +
            theme_bw() +
            theme(aspect.ratio = 1)
        
    } else {
        
        plotPCA(vst(dds, blind = FALSE), intgroup = var) +
            scale_color_brewer(palette = palette) +
            labs(color = var) +
            theme_bw() +
            theme(aspect.ratio = 1)
    }
}



deseq_volcano <- function(res, p_co, lfc_co) {
    
    res %<>%
        deseq_transform(p_co, lfc_co)
    
        ggplot(res) +
        geom_point(aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = significant),
                   alpha = 0.5) +
        scale_color_manual(values = c("#bdbdbd", "#de2d26")) +
        scale_x_continuous(limits = c(max(abs(res$log2FoldChange)) * c(-1, 1))) +
        theme_bw() +
        theme(aspect.ratio = 1)
    
}


deseq_table <- function(res, p_co, lfc_co) {
    
    res %<>%
        deseq_transform(p_co, lfc_co)
    
    res %<>%
        filter(significant == TRUE) %>%
        select(symbol:padj) %>%
        arrange(padj, abs(log2FoldChange)) %>%
        datatable() %>%
        formatRound(columns = c(2:5), digits = 3) %>%
        formatSignif(columns = c(6:7), digits = 3)
    
    return(res)
}




seq(from = 1, to = 10, length.out = )

