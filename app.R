
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
        plot_volcano(deseq2_res, input$p_co, input$lfc_co)
    }, height = 700)
    
    output$table <- DT::renderDataTable({
        deseq2_res %>%
            as.data.frame() %>%
            rownames_to_column(var = "symbol") %>%
            as_tibble() %>%
            mutate(significant = ifelse(padj <= input$p_co & abs(log2FoldChange) >= input$lfc_co, 
                                        TRUE, FALSE)) %>%
            filter(significant == TRUE) %>%
            select(symbol:padj) %>%
            arrange(padj, abs(log2FoldChange)) %>%
            datatable() %>%
            formatRound(columns = c(2:5), digits = 3) %>%
            formatSignif(columns = c(6:7), digits = 3)
    })
}

shinyApp(ui = ui, server = server)
