
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
    tab_integrate,
    tab_about
)

# server function
server <- function(input, output, session) {

    # I/O

    tempMat <- reactiveValues(mat = NULL)
    observeEvent(input$submitFiles, {
        
        # input DNA-Seq dataset
        output$ph1 <- renderText("")
        f <- as.character(input$dnafile)
        temp <- read.delim(f, sep = "	", header = TRUE, as.is = TRUE)
        tempMat$mat <- read.maf(temp, useAll = FALSE)
        geneList_sorted <- unlist(getGeneSummary(tempMat$mat)[, 1])
        sampleList_sorted <- unlist(getSampleSummary(tempMat$mat)[, 1])
        geneList <- paste0(geneList_sorted[1:20], collapse = ";")
        sampleList <- paste0(sampleList_sorted, collapse = ";")
        updateTextInput(session, inputId = "genes", value = geneList)
        updateTextInput(session, inputId = "samples", value = sampleList)
        output$ph1 <- renderText(paste("Upload dataset successfully! -", basename(input$dnafile)))

        # input RPPA dataset
        output$ph3 <- renderText("")
        rppaMeta <- makeMetaTable(rppa = input$rppafile)
        output$specifyRppaGenes <- renderUI({
            textInput("genes", "Name of genes", "EGFR;TP53;RB1;ALK;AML1;TERT;KRAS;PI3K")
        })
        output$specifyRppaSamples <- renderUI({
            textInput("samples", "Name of samples", paste(rppaMeta[["sid"]], sep = "", collapse = ";"))
        })
        output$ph3 <- renderText(paste("Upload dataset successfully! -", basename(input$rppafile)))
    })


    # DNA

    observeEvent(input$resetDnaGenesSamples, {
        selectedGenes <- unlist(strsplit(as.character(input$genes), ";"))
        selectedSamples <- unlist(strsplit(as.character(input$samples), ";"))
        maf <- subsetMaf(tempMat$mat, tsb = selectedSamples)
        
        # 1) Waterfall Plot
        output$waterfall <- renderPlot({
            oncoplot(maf = maf, genes = selectedGenes, gene_mar = 10, draw_titv = TRUE)
        })
        # 2) Summary Plot
        output$mafSummary <- renderPlot({
            plotmafSummary(maf = tempMat$mat, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
        })
        # 3) TiTv Plot
        output$titvPlot <- renderPlot({
            if(is.null(tempMat$mat)){
                tempMat$mat <- NULL
            }else{
                titv(maf = tempMat$mat, plot = TRUE, useSyn = TRUE)
            }
        })
        # 4) Lollipop Plot
        output$lollipopPlot1 <- renderPlot({
            lollipopPlot(maf = tempMat$mat, gene = selectedGenes[1], AACol = 'Protein_Change', showMutationRate = TRUE, showDomainLabel = TRUE,
                         repel = TRUE, labPosAngle = 45, labelPos = 'all')
        })
        output$lollipopPlot2 <- renderPlot({
            lollipopPlot(maf = tempMat$mat, gene = selectedGenes[2], AACol = 'Protein_Change', showMutationRate = TRUE, showDomainLabel = TRUE,
                         repel = TRUE, labPosAngle = 45, labelPos = 'all')
        })
        output$lollipopPlot3 <- renderPlot({
            lollipopPlot(maf = tempMat$mat, gene = selectedGenes[3], AACol = 'Protein_Change', showMutationRate = TRUE, showDomainLabel = TRUE,
                         repel = TRUE, labPosAngle = 45, labelPos = 'all')
        })
        # 5) Rainfall Plot
        output$rainfallPlot1 <- renderPlot({
            rainfallPlot(maf = tempMat$mat, tsb = selectedSamples[1], ref.build = "hg19", detectChangePoints = TRUE, pointSize = 0.6)
        })
        output$rainfallPlot2 <- renderPlot({
            rainfallPlot(maf = tempMat$mat, tsb = selectedSamples[2], ref.build = "hg19", detectChangePoints = TRUE, pointSize = 0.6)
        })
        output$rainfallPlot3 <- renderPlot({
            rainfallPlot(maf = tempMat$mat, tsb = selectedSamples[3], ref.build = "hg19", detectChangePoints = TRUE, pointSize = 0.6)
        })
        # 6) TCGA Compare
        output$tcgaCompare <- renderPlot({
            maf.mutload = tcgaCompare(maf = tempMat$mat, cohortName = "Input")
        })
        # 7) 
        output$somaticInteract <- renderPlot({
            somaticInteractions(maf = tempMat$mat, top = 20, pvalue = c(0.05, 0.1))
        })
        # 8) 
        output$drugInteract <- renderPlot({
            dgi = drugInteractions(maf = tempMat$mat, fontSize = 0.75)
        })
        # 9) 
        output$oncogenicPathway1 <- renderPlot({
            pathwayList <- OncogenicPathways(maf = tempMat$mat)[, 1]
        })
        output$oncogenicPathway2 <- renderPlot({
            PlotOncogenicPathways(maf = tempMat$mat, pathways = unlist(pathwayList[10]))
        })
        # 10) 
        output$vafPlot <- renderPlot({
            plotVaf(maf = tempMat$mat, vafCol = 'VAF', top = 30)
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
        deseq_ma(rnaseq[[2]])
    }, height = 600)

    output$deseq_volcano <- renderPlot({
        deseq_volcano(rnaseq[[2]], input$p_co, input$lfc_co)
    }, height = 600)

    output$deseq_table <- DT::renderDataTable({
        deseq_table(rnaseq[[2]], input$p_co, input$lfc_co)
    })

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

