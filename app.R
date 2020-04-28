
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


source("tab_design.R")
source("helpers.R")
source("visualization.R")

rnaseq <- readRDS("./data/rnaseq.rds")

hmks_hs <- gmtPathways("./data/h.all.v7.0.symbols.gmt")

ui <- navbarPage(
    
    theme = shinytheme("yeti"),
    title = "Genomic Data Explorer",
    id = "tabs",
    
    tab_upload,
    tab_dna,
    tab_rna,
    tab_scrna
)

server <- function(input, output) {

    output$ph1 <- renderText("Placeholder")
    output$ph2 <- renderText("Placeholder")
    output$ph4 <- renderText("Placeholder")
    
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
}

shinyApp(ui = ui, server = server)
