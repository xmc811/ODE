
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

deseq_ma <- function(res, ylim = c(-4,4), p_cutoff = 0.1) {
    
    resDF <- res %>%
        as.data.frame() %>%
        rownames_to_column(var = "Symbol") %>%
        as_tibble()
    
    stt <- c()
    
    stt[1] <- resDF %>%
        filter(padj < p_cutoff) %>%
        filter(log2FoldChange > 0) %>%
        nrow()
    stt[2] <- resDF %>%
        filter(padj < p_cutoff) %>%
        filter(log2FoldChange < 0) %>%
        nrow()
    stt[3] <- resDF %>%
        filter(is.na(pvalue)) %>%
        nrow()
    
    pct <- paste0(round(stt / nrow(resDF) * 100,1), "%")
    
    resDF %>%
        drop_na() %>%
        mutate(shape = ifelse(log2FoldChange > ylim[2] | log2FoldChange < ylim[1], TRUE, FALSE),
               log2FoldChange = replace(log2FoldChange, log2FoldChange > ylim[2], ylim[2] * 1.05),
               log2FoldChange = replace(log2FoldChange, log2FoldChange < ylim[1], ylim[1] * 1.05)) %>%
        arrange(-padj) %>%
        ggplot() +
        geom_point(aes(x = baseMean, 
                       y = log2FoldChange, 
                       color = ifelse(padj < p_cutoff, TRUE, FALSE), 
                       shape = shape,
                       alpha = factor(ifelse(padj < p_cutoff, 1, 0.8)))) +
        scale_color_manual(values = c("black", "#e31a1c"), 
                           labels = c(paste("\u2265", p_cutoff),paste("<",p_cutoff)), 
                           name = "Adjusted p-value") +
        scale_shape_manual(values = c(16, 17), guide = FALSE) +
        scale_alpha_manual(values = c(0.25, 1), guide = FALSE) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        theme_bw() +
        labs(y = expression(Log[2]~Fold~Change), x = "Mean of Normalized Count") +
        theme(legend.position = "bottom",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14)) +
        geom_hline(yintercept = 0, color = "#80b1d3", size = 1.5, alpha = 0.5) +
        annotation_logticks(sides = "b") +
        annotate("text", 
                 x = 10 ^ (ceiling(log10(max(resDF$baseMean))) * 0.75), 
                 y = -Inf, 
                 hjust = 0,
                 vjust = -0.5,
                 label = paste0("Out of ", nrow(res), " genes:", "\n",
                                "Adjusted p-value < ", p_cutoff, "\n",
                                "LFC > 0 (up): ", stt[1], ", ", pct[1], "\n",
                                "LFC < 0 (down): ", stt[2], ", ", pct[2], "\n",
                                "Gene outliers: ", stt[3], ", ", pct[3]))
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

deseq_gsea <- function(res) {
    
    res <- deseq_to_stat(res)
    
    gsea_res <- fgsea(pathways = hmks_hs, 
                      stats = res, 
                      nperm = 50000)
    
    gsea_res %>%
        mutate(pathway = str_remove(string = pathway, pattern = "HALLMARK_")) %>%
        mutate(color = -log10(padj) * ifelse(padj <= 0.05, 1, 0) * ifelse(NES > 0, 1, -1)) %>%
        ggplot() +
        geom_bar(aes(x = reorder(pathway, NES), y = NES, fill = color), stat = "identity") +
        scale_fill_gradient2(high = "#d7301f", 
                              mid = "#f0f0f0",
                              low = "#0570b0",
                              midpoint = 0) +
        coord_flip() +
        labs(x = "Pathway",
             fill = "-log10 adjusted p-value") +
        theme_bw()
    
}
