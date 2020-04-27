
library(shiny)
library(shinythemes)
library(tidyverse)
library(magrittr)
library(DESeq2)
library(DT)


source("tab_design.R")
source("helpers.R")
source("visualization.R")

deseq2_res <- readRDS("./data/deseq2_res_list.rds")[[1]]



ui <- navbarPage(
    
    theme = shinytheme("yeti"),
    title = "Omics Data Explorer",
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
    
    output$volcano <- renderPlot({
        deseq_volcano(deseq2_res, input$p_co, input$lfc_co)
    }, height = 700)
    
    output$table <- DT::renderDataTable({
        deseq_table(deseq2_res, input$p_co, input$lfc_co)
    })
}

shinyApp(ui = ui, server = server)
