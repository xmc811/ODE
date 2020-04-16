
tab_upload <- tabPanel(
    
    title = "Upload",
    fluid = TRUE,
    value = "v_up",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(

                    h4("Please upload your data:"),
                    fileInput(inputId = "data",
                              label = "Data",
                              buttonLabel = "Browse...")

        ),
        
        mainPanel(
            textOutput(outputId = "ph1")
        )
    )
    
)

tab_dna <- tabPanel(
    
    title = "DNA-seq",
    fluid = TRUE,
    value = "v_dna",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            
            
            
        ),
        
        mainPanel(
            textOutput(outputId = "ph2")
        )
    )
    
)


tab_rna <- tabPanel(
    
    title = "RNA-seq",
    fluid = TRUE,
    value = "v_rna",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            
            
            
        ),
        
        mainPanel(
            textOutput(outputId = "ph3")
        )
    )
    
)

tab_scrna <- tabPanel(
    
    title = "scRNA-seq",
    fluid = TRUE,
    value = "v_scrna",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            
            
            
        ),
        
        mainPanel(
            textOutput(outputId = "ph4")
        )
    )
    
)
