
library(shiny)
library(shinythemes)
library(tidyverse)
library(magrittr)
library(DESeq2)
library(DT)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(scales)
library(fgsea)
library(maftools)


rnaseq <- readRDS("./data/rnaseq.rds")
hmks_hs <- gmtPathways("./data/h.all.v7.0.symbols.gmt")

# example input files

dnaMaf <- "./data/DNA_filtered__mutect_Pindel_GM_pipe_201002.maf"
rnaTsv <- "./data/RNA_merged_HTseq_count.tsv"
rppa <- "./data/IACS__RPPA.txt"

source("tab_design.R")
source("helpers.R")
source("visualization.R")

# user interface
ui <- navbarPage(
    
    theme = shinytheme("yeti"),
    title = "Genomic Data Explorer",
    id = "tabs",
    
    tab_upload,
    tab_dna,
    tab_rna,
    tab_rppa,
    tab_scrna,
    tab_integrate
)

# server function
server <- function(input, output) {
    
    # I/O
    
    tempMat <- reactiveValues(mat = NULL)
    observeEvent(input$submitFiles, {
        
        meta <- makeMetaTable(dnaMaf = input$dnafile)
        
        output$specifyDnaGenes <- renderUI({
            textInput("genes", "Name of genes", "EGFR;TP53;RB1;ALK;AML1;TERT;KRAS;PI3K")
        })
        
        output$specifyDnaSamples <- renderUI({
            textInput("samples", "Name of samples", paste(meta[["sid"]], sep = "", collapse = ";"))
        })
        
        
        meta <- makeMetaTable(rppa = input$rppafile)
        
        output$specifyRppaGenes <- renderUI({
            textInput("genes", "Name of genes", "EGFR;TP53;RB1;ALK;AML1;TERT;KRAS;PI3K")
        })
        
        output$specifyRppaSamples <- renderUI({
            textInput("samples", "Name of samples", paste(meta[["sid"]], sep = "", collapse = ";"))
        })
        
        output$ph1 <- renderText("Upload datasets successfully!")
    })
    
    
    # DNA
    
    observeEvent(input$resetDnaGenesSamples, {
        f <- as.character(input$dnafile)
        temp <- read.delim(f, sep = "	", header = TRUE, as.is = TRUE)
        temp <- temp[temp[, "Tumor_Sample_Barcode"] %in% unlist(strsplit(as.character(input$samples), ";")), ]
        temp <- temp[temp[, "Hugo_Symbol"] %in% unlist(strsplit(as.character(input$genes), ";")), ]
        tempMat$mat <- read.maf(temp)
        
        output$waterfall <- renderPlot({
            if(is.null(tempMat$mat)){
                tempMat$mat <- NULL
            }else{
                if(class(tempMat$mat) == "MAF"){
                    require(maftools)
                    oncoplot(tempMat$mat, draw_titv = TRUE)
                }
            }
        })
        
        output$mafSummary <- renderPlot({
            if(is.null(tempMat$mat)){
                tempMat$mat <- NULL
            }else{
                if(class(tempMat$mat) == "MAF"){
                    require(maftools)
                    plotmafSummary(maf = tempMat$mat, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
                }
            }
        })
        
        output$vafPlot <- renderPlot({
            if(is.null(tempMat$mat)){
                tempMat$mat <- NULL
            }else{
                if(class(tempMat$mat) == "MAF"){
                    require(maftools)
                    plotVaf(maf = tempMat$mat, vafCol = 'VAF', top = 30)
                }
            }
        })
        
    })
    
    
    # RNA
    
    output$deseq_hm <- renderPlot({
        deseq_heatmap(rnaseq[[1]], 
                      input$pca_var,
                      input$palette_con,
                      input$palette_dir)
    }, height = 700, width = 800)
    
    output$deseq_pca <- renderPlot({
        if (is.numeric(rnaseq[[1]]@colData[[input$pca_var]])) {
            deseq_pca(rnaseq[[1]], 
                      input$pca_var, 
                      input$palette_con,
                      ifelse(input$palette_dir, 1, -1))
        } else {
            deseq_pca(rnaseq[[1]], 
                      input$pca_var, 
                      input$palette_cat)
        }
    }, height = 600)
    
    output$deseq_ma <- renderPlot({
        deseq_ma(rnaseq[[2]],
                 input$p_co, 
                 input$lfc_co,
                 input$lfc_plot_lim)
    }, height = 600)
    
    output$deseq_volcano <- renderPlot({
        deseq_volcano(rnaseq[[2]], 
                      input$p_co, 
                      input$lfc_co,
                      input$p_plot_lim,
                      input$lfc_plot_lim)
    }, height = 600)
    
    output$deseq_table <- DT::renderDataTable({
        deseq_table(rnaseq[[2]], input$p_co, input$lfc_co)
    })
    
    output$deseq_box <- renderPlot({
        deseq_box(rnaseq[[1]], 
                      rnaseq[[2]], 
                      input$pca_var, 
                      input$palette_cat)
    }, height = 700, width = 800)
    
    output$deseq_cluster <- renderPlot({
        deseq_cluster(rnaseq[[1]], 
                      rnaseq[[2]], 
                      input$palette_con, 
                      input$palette_dir)
    }, height = 700)
    
    output$deseq_gsea <- renderPlot({
        deseq_gsea(rnaseq[[2]])
    }, height = 600)
    
    
    # RPPA

    observeEvent(input$resetRppaGenesSamples, {
        f <- as.character(input$rppafile)
        temp <- read.delim(f, sep = "	", header = TRUE, as.is = TRUE)
        rownames(temp) <- rown <- temp[, "Sample.Name"]
        temp <- temp[, -which(colnames(temp) == "Sample.Name")]
        temp <- temp[, rownames(temp) %in% unlist(strsplit(as.character(input$samples), ";"))]
        temp <- apply(apply(temp, 2, function(x) log2(as.numeric(x) + 1)), 2, function(x) (x - mean(x, na.rm = TRUE))/sd(x))
        rownames(temp) <- rown
        temp <- data.matrix(t(temp))
        class(temp) <- "rppa"
        tempMat$mat <- temp

        output$clusterRppa <- renderPlot({
            if(is.null(tempMat$mat)){
                tempMat$mat <- NULL
            }else{
                if(class(tempMat$mat) == "rppa"){
                    require(gplots)
                    cols = colorRampPalette(c("forestgreen", "yellow", "red"))(n =299)
                    plot(hclust(dist(tempMat$mat)))
                    heatmap.2(tempMat$mat, col = cols, tracecol = NA, dendgrogram = "column")
                }
            }
        })
    })

    output$ph4 <- renderText("Under Construction...")
    output$ph5 <- renderText("Under Construction...")
    
}

shinyApp(ui = ui, server = server)
