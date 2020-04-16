
library(shiny)
library(shinythemes)


source("tab_design.R")

ui <- navbarPage(
    
    theme = shinytheme("cosmo"),
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
    output$ph3 <- renderText("Placeholder")
    output$ph4 <- renderText("Placeholder")
}

shinyApp(ui = ui, server = server)
