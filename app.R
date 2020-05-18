

source("loading.R")
source("helpers.R")
source("tab_design.R")
source("visualization.R")

options(shiny.maxRequestSize = 100*1024^2)

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
    tab_about,
    
    tags$head(tags$link(rel="stylesheet", 
                        type="text/css", 
                        href="style.css"))
)

# server function
server <- function(input, output, session) {

    # I/O

    tempMat <- reactiveValues(mat = NULL)
    
    observeEvent(input$dnafile, {
        output$ph1 <- NULL
        output$ph2 <- NULL
        output$ph3 <- NULL
        output$waterfall <- NULL
        output$mafSummary <- NULL
        output$titvPlot <- NULL
        output$lollipopPlot1 <- NULL
        output$lollipopPlot2 <- NULL
        output$lollipopPlot3 <- NULL
        output$rainfallPlot1 <- NULL
        output$rainfallPlot2 <- NULL
        output$rainfallPlot3 <- NULL
        output$tcgaCompare <- NULL
        output$somaticInteract <- NULL
        output$drugInteract <- NULL
        output$oncogenicPathway1 <- NULL
        output$oncogenicPathway2 <- NULL
        output$vafPlot <- NULL
    })
    
    observeEvent(input$rppafile, {
        output$ph1 <- NULL
        output$ph2 <- NULL
        output$ph3 <- NULL
        output$clusterRppa <- NULL
    })
    
    observeEvent(input$submitFiles, {
        
        # input DNA-Seq dataset
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
            pathwayList_sorted <- OncogenicPathways(maf = tempMat$mat)[, 1]
            pathwayList <- paste0(unlist(pathwayList_sorted), collapse = ";")
            updateTextInput(session, inputId = "pathways", value = pathwayList)
        })
        output$oncogenicPathway2 <- renderPlot({
            selectedPathways <- unlist(strsplit(as.character(input$pathways), ";"))
            PlotOncogenicPathways(maf = tempMat$mat, pathways = unlist(selectedPathways[1]))
        })
        # 10) 
        output$vafPlot <- renderPlot({
            plotVaf(maf = tempMat$mat, vafCol = input$vafColumn, top = 30)
        })

    })


    # RNA
    
    observeEvent(input$rna_start, {
        
        library(DESeq2)
        
        rna_input <- if(input$data_source == "Example") {
            reactive({readRDS("./data/rnaseq.rds")})
        } else {
            reactive({
                validate(
                    need(input$rna_input, 
                         "Please Upload Data")
                )
                readRDS(input$rna_input$datapath)})
        }
            
        output$rna_var <- renderUI({
            validate(
                need(try(rna_input()), "")
            )
            list(
            h4("Variable"),
            selectInput(inputId = "pca_var", 
                        label = "Categorical Variable for Heatmap, PCA, and Gene Boxplot",
                        choices = colnames(rna_input()[[1]]@colData))
            )
        })
        
        rna_plot_height <- reactive({
            validate(
                need(input$rna_plot_height < 4000, "Plot height shouldn't exceed 4000px.")
                )
            return(input$rna_plot_height)
        })
        
        rna_plot_width <- reactive({
            validate(
                need(input$rna_plot_width < 4000, "Plot width shouldn't exceed 4000px.")
            )
            return(input$rna_plot_width)
        })

        output$deseq_hm <- renderPlot({
            validate(
                need(input$pca_var, "Please Upload Data")
            )
            deseq_heatmap(rna_input()[[1]], 
                          input$pca_var,
                          input$palette_con,
                          input$palette_dir)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_pca <- renderPlot({
            validate(
                need(input$pca_var, "Please Upload Data")
            )
            if (is.numeric(rna_input()[[1]]@colData[[input$pca_var]])) {
                deseq_pca(rna_input()[[1]], 
                          input$pca_var, 
                          input$palette_con,
                          ifelse(input$palette_dir, 1, -1))
            } else {
                deseq_pca(rna_input()[[1]], 
                          input$pca_var, 
                          input$palette_cat)
            }
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_ma <- renderPlot({
            deseq_ma(rna_input()[[2]],
                     input$p_co, 
                     input$lfc_co,
                     input$lfc_plot_lim)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_volcano <- renderPlot({
            deseq_volcano(rna_input()[[2]], 
                          input$p_co, 
                          input$lfc_co,
                          input$p_plot_lim,
                          input$lfc_plot_lim)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_table <- DT::renderDataTable({
            deseq_table(rna_input()[[2]], 
                        input$p_co, 
                        input$lfc_co)
        })
        
        rna_genes <- eventReactive(
            
            c(input$rna_gene_read1,
              input$rna_gene_read2,
              input$rna_genes_file),
            
            {
            
            if(input$rna_gene_ls_src == 'Use Top Genes') {

                get_rna_genes(rna_input()[[2]])[1:input$rna_gene_num]
                
            } else if (input$rna_gene_ls_src == 'Manual Input') {
                
                validate(
                        need(input$rna_genes_man, "Please Input Gene List")
                        )
                parse_rna_genes(input$rna_genes_man)
                
            } else {
                
                validate(
                        need(input$rna_genes_file, "Please Upload Gene List")
                        )
                readLines(input$rna_genes_file$datapath)
            }})
        
            
        output$deseq_box <- renderPlot({
            validate(
                need(input$pca_var, "Please Upload Data")
            )
            deseq_box(rna_input()[[1]], 
                      rna_genes(),
                      input$pca_var, 
                      input$palette_cat)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_cluster <- renderPlot({
            deseq_cluster(rna_input()[[1]], 
                          rna_genes(),
                          input$palette_con, 
                          input$palette_dir)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        
        rna_pathways <- reactive({
                
                if(input$rna_pathway_src == 'Use Hallmarks') {
                    
                    hmks_hs
                    
                } else {
                    
                    validate(
                        need(input$rna_pathway_file, "Please Upload Pathway File")
                    )
                    read_csv(input$rna_pathway_file$datapath, 
                             col_names = FALSE) %>%
                        df_to_signature()
                    
                }})
        
        output$deseq_gsea <- renderPlot({
            deseq_gsea(rna_input()[[2]],
                       rna_pathways())
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
    })
    
    
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

