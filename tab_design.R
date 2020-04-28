
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
            h4("PCA"),
            selectInput(inputId = "pca_var", 
                        label = "PCA Plot Variable",
                        choices = colnames(rnaseq[[1]]@colData)),
            br(),
            
            h4("Differential Gene Expression Parameters"),
            numericInput("p_co", label = "Adjusted p-value Cutoff", value = 0.05),
            numericInput("lfc_co", label = "Log2 Fold Change Cutoff", value = 1),
            br(),
            
            h4("Plotting Parameters"),
            selectInput(inputId = "palette_cat", 
                            label = "Palette (Categorical Variable)",
                            choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                            selected = "Set1"),
            selectInput(inputId = "palette_con", 
                            label = "Palette (Continuous Variable)",
                            choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                            selected = "BuPu"),
            checkboxInput(inputId = "palette_dir",
                          label = "Reverse Scale Color Direction",
                          value = FALSE),
            br()
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel(
                    title = "Heatmap",
                    br(),
                    plotOutput("deseq_hm", width = "100%")
                ),
                tabPanel(
                    title = "PCA",
                    br(),
                    plotOutput("deseq_pca", width = "100%")
                ),
                tabPanel(
                    title = "MA Plot",
                    br(),
                    plotOutput("deseq_ma", width = "100%")
                ),
                tabPanel(
                    title = "Volcano Plot",
                    br(),
                    plotOutput("deseq_volcano", width = "100%")
                ),
                tabPanel(
                    title = "Table",
                    br(),
                    DT::dataTableOutput("deseq_table")
                ),
                tabPanel(
                    title = "Gene-wise Boxplot",
                    br(),
                    plotOutput("deseq_genebox", width = "100%")
                ),
                tabPanel(
                    title = "GSEA",
                    br(),
                    plotOutput("deseq_gsea", width = "100%")
                )
                
            )
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
