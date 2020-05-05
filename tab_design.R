
tab_upload <- tabPanel(

    title = "Upload",
    fluid = TRUE,
    value = "v_up",

    sidebarLayout(

        sidebarPanel = sidebarPanel(

            h4("Please upload your data:"),
            fileInput(inputId = "data",
                      label = "Data",
                      buttonLabel = "Browers..."),
            br(),
            h3("Data sets used:"),
            textInput("rdsfile", "RDS data file name"),
            textInput("dnafile", "DNA data file name", dnaMaf),
            textInput("rnafile", "RNA data file name", rnaTsv),
            textInput("rppafile", "RPPA data file name", rppa),

            actionButton("submitFiles", "Submit")
        ),

        mainPanel(
            br(),
            textOutput(outputId = "ph1"),
            br(),
            textOutput(outputId = "ph2"),
            br(),
            textOutput(outputId = "ph3")
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
            textInput(inputId = "genes", label = "Name of genes", value = "EGFR;TP53;RB1;ALK;AML1;TERT;KRAS;PI3K"),
            textInput(inputId = "samples", label = "Name of samples"),
            actionButton(inputId = "resetDnaGenesSamples", label = "Run")
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
                    title = "TiTv Plot",
                    br(),
                    plotOutput("titvPlot")
                ),
                tabPanel(
                    title = "Lollipop Plot",
                    br(),
                    plotOutput("lollipopPlot1"),
                    br(),
                    plotOutput("lollipopPlot2"),
                    br(),
                    plotOutput("lollipopPlot3")
                ),
                tabPanel(
                    title = "Rainfall Plot",
                    br(),
                    plotOutput("rainfallPlot1"),
                    br(),
                    plotOutput("rainfallPlot2"),
                    br(),
                    plotOutput("rainfallPlot3")
                ),
                tabPanel(
                    title = "TCGA Compare",
                    br(),
                    plotOutput("tcgaCompare")
                ),
                tabPanel(
                    title = "Somatic Interactions",
                    br(),
                    plotOutput("somaticInteract", width = "100%", height = 700)
                ),
                tabPanel(
                    title = "Drug-Gene Interactions",
                    br(),
                    plotOutput("drugInteract")
                ),
                tabPanel(
                    title = "Oncogenic Signaling Pathways",
                    br(),
                    plotOutput("oncogenicPathway1"),
                    br(),
                    plotOutput("oncogenicPathway2")
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

tab_about <- tabPanel(
    
    title = "About",
    fluid = TRUE,
    value = "v_about",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            
            
            
        ),
        
        mainPanel(
            br(),
            h4("Authors:"),
            br(),
            h4("Mingchu Xu"),
            br(),
            h4("Xiaogang (Sean) Wu"),
            br(),
            h4("Jianhua (John) Zhang")
        )
    )
    
)
