
tab_upload <- tabPanel(

    title = "Upload",
    fluid = TRUE,
    value = "v_up",

    sidebarLayout(

        sidebarPanel = sidebarPanel(

            h4("Please upload your data:"),
            fileInput(inputId = "data",
                      label = "Data",
                      buttonLabel = "Browse..."),
            br(),
            h3("Data sets used:"),
            textInput("rdsfile", "RDS data file name"),
            textInput("dnafile", "DNA data file name", dnaMaf),
            textInput("rnafile", "RNA data file name", rnaTsv),
            textInput("rppafile", "RPPA data file name", rppa),

            actionButton("submitFiles", "Submit")
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
            h4("Specify genes and samples:"),
            uiOutput("specifyDnaGenes"),
            uiOutput("specifyDnaSamples"),
            actionButton("resetDnaGenesSamples", "Reset")
        ),

        mainPanel(
            tabsetPanel(
                tabPanel(
                    title = "Waterfall Plot",
                    br(),
                    plotOutput("waterfall")
                ),
                tabPanel(
                    title = "Summary Plot",
                    br(),
                    plotOutput("mafSummary")
                ),
                tabPanel(
                    title = "VAF Plot",
                    br(),
                    plotOutput("vafPlot")
                )
            )
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
            splitLayout(numericInput("p_co", label = "Adjusted P-value Cutoff", value = 0.05),
                        numericInput("lfc_co", label = "Log2 Fold Change Cutoff", value = 1)),
            br(),
            h4("Gene List"),
            fileInput(inputId = "rna_genes",
                      label = NULL,
                      buttonLabel = "Browse..."),
            br(),
            h4("Plotting Parameters"),
            splitLayout(selectInput(inputId = "palette_cat", 
                                    label = "Categorical Palette",
                                    choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                    selected = "Set1"),
                        selectInput(inputId = "palette_con", 
                                    label = "Continuous Palette",
                                    choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                    selected = "BuPu")),
            checkboxInput(inputId = "palette_dir",
                          label = "Reverse Scale Color Direction",
                          value = FALSE),
            splitLayout(numericInput("p_plot_lim", label = "Adjusted P-value Squash", value = 5),
                        numericInput("lfc_plot_lim", label = "Log2 Fold Change Squash", value = 5)),
            tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))),
            br()
        ),
        
        mainPanel(
            tabsetPanel(
                id = "rna_panel",
                tabPanel(
                    title = "Heatmap",
                    br(),
                    plotOutput("deseq_hm")
                ),
                tabPanel(
                    title = "PCA",
                    br(),
                    plotOutput("deseq_pca", width = "100%")
                ),
                tabPanel(
                    title = "MA Plot",
                    br(),
                    plotOutput("deseq_ma", width = "100%"),
                ),
                tabPanel(
                    title = "Volcano Plot",
                    br(),
                    plotOutput("deseq_volcano"),
                    br(),
                    uiOutput("pt_color")
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

tab_rppa <- tabPanel(

    title = "RPPA",
    fluid = TRUE,
    value = "v_rppa",

    sidebarLayout(

        sidebarPanel = sidebarPanel(
            h4("Specify genes and samples:"),
            uiOutput("specifyRppaGenes"),
            uiOutput("specifyRppaSamples"),
            actionButton("resetRppaGenesSamples", "Reset")
        ),

        mainPanel(
            tabsetPanel(
                tabPanel(
                    title = "Clustering",
                    br(),
                    plotOutput("clusterRppa")
                ),
                tabPanel(
                    title = "PCA Plot",
                    br(),
                    h4("Under Construction ...")
                ),
                tabPanel(
                    title = "Pathway Analysis",
                    br(),
                    h4("Under Construction ...")
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

tab_integrate <- tabPanel(

    title = "Integration",
    fluid = TRUE,
    value = "v_integrate",

    sidebarLayout(

        sidebarPanel = sidebarPanel(



        ),

        mainPanel(
            textOutput(outputId = "ph5")
        )
    )

)


