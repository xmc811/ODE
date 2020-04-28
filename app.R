
library(shiny)
library(shinythemes)
library(tidyverse)
library(magrittr)
library(DESeq2)
library(DT)

# example input files
rnaRds <- "./data/deseq2_res_list.rds"
dnaMaf <- "./data/DNA_filtered__mutect_Pindel_GM_pipe_201002.maf"
rnaTsv <- "./data/RNA_merged_HTseq_count.tsv"
rppa <- "./data/IACS__RPPA.txt"

source("./tab_design.R")
source("./helpers.R")
source("./visualization.R")

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

    tempMat <- reactiveValues(mat = NULL)
    observeEvent(input$submitFiles, {

        meta <- makeMetaTable(dnaMaf = input$dnafile)
        output$specifyDnaGenes <- renderUI({
            textInput("genes", "Name of genes", "EGFR;TP53;RB1;ALK;AML1;TERT;KRAS;PI3K")
        })
        output$specifyDnaSamples <- renderUI({
            textInput("samples", "Name of samples", paste(meta[["sid"]], sep = "", collapse = ";"))
        })

        meta <- makeMetaTable(rnaTsv = input$rnafile)
        output$specifyRnaGenes <- renderUI({
            textInput("genes", "Name of genes", "EGFR;TP53;RB1;ALK;AML1;TERT;KRAS;PI3K")
        })
        output$specifyRnaSamples <- renderUI({
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

    observeEvent(input$resetRnaGenesSamples, {
        f <- as.character(input$rnafile)
        temp <- read.delim(f, sep = "	", header = TRUE, as.is = TRUE)
        temp <- temp[temp[, "symbol"] %in% unlist(strsplit(as.character(input$genes), ";")), , drop = FALSE]
        temp <- temp[!duplicated(temp[, "symbol"]), ]
        rown <- temp[, "symbol"]
        temp <- temp[, colnames(temp) %in% unlist(strsplit(as.character(input$samples), ";"))]
        temp <- apply(apply(temp, 2, function(x) log2(as.numeric(x) + 1)), 2, function(x) (x - mean(x, na.rm = TRUE))/sd(x))
        rownames(temp) <- rown
        temp <- as.matrix(temp)
        class(temp) <- "rna"
        tempMat$mat <- temp

        output$clusterRna <- renderPlot({
            if(is.null(tempMat$mat)){
                tempMat$mat <- NULL
            }else{
                if(class(tempMat$mat) == "rna"){
                    require(gplots)
                    cols = colorRampPalette(c("forestgreen", "yellow", "red"))(n =299)
                    plot(hclust(dist(tempMat$mat)))
                    heatmap.2(tempMat$mat, col = cols, tracecol = NA, dendgrogram = "column")
                }
            }
        })

        deseq2_res <- readRDS(input$rdsfile)[[1]]

        output$volcano <- renderPlot({
            deseq_volcano(deseq2_res, input$p_co, input$lfc_co)
        }, height = 700)

        output$table <- DT::renderDataTable({
            deseq_table(deseq2_res, input$p_co, input$lfc_co)
        })
    })

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
