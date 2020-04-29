
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