
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
        mutate(significant = ifelse(padj <= p_co & log2FoldChange >= lfc_co, 
                                    "Up",
                                    ifelse(padj <= p_co & log2FoldChange <= -lfc_co, 
                                           "Down", 
                                           "Not Sig")))
    return(res)
}

get_mtx_dds <- function(dds, res, raw = F) {
    
    if (raw) {
        vsd <- counts(dds)
    } else {
        vsd <- vst(dds, blind = FALSE)
        vsd <- as.matrix(vsd@assays@data[[1]])
    }
    
    res %<>%
        deseq_transform(p_co = 1, lfc_co = 0)
    
    genes <- res %>%
        arrange(desc(abs(log2FoldChange))) %>%
        head(50) %>%
        pull(symbol)
    
    mtx <- vsd[genes,]
    
    return(mtx)
}

get_nm_count_dds <- function(dds, res, var) {

    res %<>%
        deseq_transform(p_co = 1, lfc_co = 0)
    
    genes <- res %>%
        arrange(desc(abs(log2FoldChange))) %>%
        head(10) %>%
        pull(symbol)
    
    df_list <- list()
    
    for (i in seq_along(1:length(genes))) {
        
        d <- plotCounts(dds, 
                        gene = genes[i], 
                        intgroup = var, 
                        returnData = TRUE)
        
        d %<>%
            rownames_to_column() %>%
            as_tibble() %>%
            mutate(symbol = genes[i])
        
        df_list[[i]] <- d
    }
    
    df <- do.call("rbind", df_list)
    
    return(df)
}


mtx_rescale <- function(mtx) {
    
    mtx2 <- mtx
    
    for (i in 1:nrow(mtx)) {
        mtx2[i,] <- (mtx[i,] - min(mtx[i,]))/(max(mtx[i,]) - min(mtx[i,])) * 2 - 1
    }
    return(mtx2)
}

mtx_logtrans <- function(mtx, b) {
    
    mtx2 <- mtx
    
    for (i in 1:nrow(mtx)) {
        mtx2[i,] <- logb(mtx[i,], b)
    }
    
    return(mtx2)
}

deseq_to_stat <- function(res) {
    
    stat <- res$stat
    names(stat) <- rownames(res)
    return(stat)
}


makeMetaTable <- function(dnaMaf, rnaTsv, rppa){
    metaData <- NULL
    sid <- NULL
    if(!missing(dnaMaf)){
        mat <- read.delim(dnaMaf, sep = "\t", header = TRUE, as.is = TRUE)
        sid <- getSampleIds(mat, what = "DNA")
        metaData <- rbind(metaData, c(type = "DNA", file = dnaMaf))
    }
    if(!missing(rnaTsv)){
        mat <- read.delim(rnaTsv, sep = "\t", header = TRUE, as.is = TRUE)
        sid <- getSampleIds(mat, what = "RNA")
        metaData <- rbind(metaData, c(type = "RNA", file = rnaTsv))
    }
    if(!missing(rppa)){
        mat <- read.delim(rppa, sep = "\t", header = TRUE, as.is = TRUE)
        sid <- getSampleIds(mat, what = "RPPA")
        metaData <- rbind(metaData, c(type = "RPPA", file = rppa))
    }
    if(is.null(sid)){
        stop("DAN or RNA data set is required")
    }
    
    return(list(sid = sid, meta = metaData))
}

getSampleIds <- function(dataMat, what = "DNA"){
    if(toupper(what) == "DNA"){
        return(unique(dataMat[, "Tumor_Sample_Barcode"]))
    }
    if(toupper(what == "RNA")){
        return(colnames(dataMat)[-c(1, 2)])
    }
    if(toupper(what == "RPPA")){
        return(dataMat[, "Sample.Name"])
    }
    stop("DAN or RNA data set is required")
    
    return(invisible())
}