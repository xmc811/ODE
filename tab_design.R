
dnaMaf <- "./data/DNA_filtered__mutect_Pindel_GM_pipe_201002.maf"
rnaTsv <- "./data/RNA_merged_HTseq_count.tsv"
rppa <- "./data/IACS__RPPA.txt"

rnaseq <- readRDS("./data/rnaseq.rds")
hmks_hs <- gmtPathways("./data/h.all.v7.0.symbols.gmt")
rna_genes <- readLines("./data/rna_genes.txt")
rna_pathways <- read_csv("./data/rna_pathways.csv", col_names = FALSE) %>% df_to_signature()

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
            conditionalPanel(condition="input.dnaTabSelected==9",
                             textInput(inputId = "pathways", label = "Name of pathways", value = "RTK-RAS")),
            conditionalPanel(condition="input.dnaTabSelected==10",
                             textInput(inputId = "vafColumn", label = "Column of VAF", value = "VAF")),
            actionButton(inputId = "resetDnaGenesSamples", label = "Run")
        ),

        mainPanel(
            tabsetPanel( id = "dnaTabSelected",
                tabPanel( value = 1,
                    title = "Waterfall Plot",
                    br(),
                    plotOutput("waterfall")
                ),
                tabPanel( value = 2,
                    title = "Summary Plot",
                    br(),
                    plotOutput("mafSummary")
                ),
                tabPanel( value = 3,
                    title = "TiTv Plot",
                    br(),
                    plotOutput("titvPlot")
                ),
                tabPanel( value = 4,
                    title = "Lollipop Plot",
                    br(),
                    plotOutput("lollipopPlot1"),
                    br(),
                    plotOutput("lollipopPlot2"),
                    br(),
                    plotOutput("lollipopPlot3")
                ),
                tabPanel( value = 5,
                    title = "Rainfall Plot",
                    br(),
                    plotOutput("rainfallPlot1"),
                    br(),
                    plotOutput("rainfallPlot2"),
                    br(),
                    plotOutput("rainfallPlot3")
                ),
                tabPanel( value = 6,
                    title = "TCGA Compare",
                    br(),
                    plotOutput("tcgaCompare")
                ),
                tabPanel( value = 7,
                    title = "Somatic Interactions",
                    br(),
                    plotOutput("somaticInteract", width = "100%", height = 700)
                ),
                tabPanel( value = 8,
                    title = "Drug-Gene Interactions",
                    br(),
                    plotOutput("drugInteract")
                ),
                tabPanel( value = 9,
                    title = "Oncogenic Signaling Pathways",
                    br(),
                    plotOutput("oncogenicPathway1"),
                    br(),
                    plotOutput("oncogenicPathway2")
                ),
                tabPanel( value = 10,
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
            h4("RNA-seq Data"),
            splitLayout(checkboxInput(inputId = "rna_use_sample",
                                      label = "Use sample data",
                                      value = TRUE),
                        actionButton(inputId = "rna_start", 
                                     label = "Launch",
                                     icon = icon("chart-bar")),
                        cellWidths = c("50%", "50%")),
            tags$style(type='text/css', 
                       "#rna_start {float:right; margin-right: 5px;}"),
            fileInput(inputId = "rna_input",
                      label = NULL,
                      buttonLabel = "Browse.."),
            
            conditionalPanel(
                condition = "input.rna_panel == 1 || 
                            input.rna_panel == 2 || 
                            input.rna_panel == 6",
                
                uiOutput("rna_var"),
            ),
            
            
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4 || 
                            input.rna_panel == 5",
                
                h4("Differential Gene Expression Parameters"),
                splitLayout(numericInput("p_co", 
                                         label = "Adjusted P-value Cutoff", 
                                         value = 0.05),
                            numericInput("lfc_co", 
                                         label = "Log2 Fold Change Cutoff", 
                                         value = 1))
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 6 || 
                            input.rna_panel == 7",
                
                h4("Gene List"),
                splitLayout(fileInput(inputId = "rna_genes",
                                      label = NULL,
                                      buttonLabel = "Browse.."),
                            checkboxInput(inputId = "rna_top_gene",
                                          label = "Use top genes",
                                          value = TRUE),
                            cellWidths = c("70%", "30%"))
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 8",
                
                h4("Pathway List"),
                splitLayout(fileInput(inputId = "rna_pathways",
                                      label = NULL,
                                      buttonLabel = "Browse.."),
                            checkboxInput(inputId = "rna_hallmark",
                                          label = "Use Hallmarks",
                                          value = TRUE),
                            cellWidths = c("70%", "30%")),
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 1 || 
                            input.rna_panel == 2 ||
                            input.rna_panel == 6 || 
                            input.rna_panel == 7",
                
                splitLayout(selectInput(inputId = "palette_cat", 
                                        label = "Categorical Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                        selected = "Set2"),
                            selectInput(inputId = "palette_con", 
                                        label = "Continuous Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                        selected = "Spectral")),
                
                checkboxInput(inputId = "palette_dir",
                              label = "Reverse Scale Color Direction",
                              value = FALSE)
            ),
            
            
            
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4",
                
                splitLayout(numericInput("p_plot_lim", 
                                         label = "Adjusted P-value Squash", 
                                         value = 5),
                            numericInput("lfc_plot_lim", 
                                         label = "Log2 Fold Change Squash", 
                                         value = 5)),
            ),
            
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
                    value = 1,
                    title = "Heatmap",
                    br(),
                    plotOutput("deseq_hm")
                ),
                tabPanel(
                    value = 2,
                    title = "PCA",
                    br(),
                    plotOutput("deseq_pca", width = "100%")
                ),
                tabPanel(
                    value = 3,
                    title = "MA Plot",
                    br(),
                    plotOutput("deseq_ma", width = "100%"),
                ),
                tabPanel(
                    value = 4,
                    title = "Volcano Plot",
                    br(),
                    plotOutput("deseq_volcano"),
                ),
                tabPanel(
                    value = 5,
                    title = "Table",
                    br(),
                    DT::dataTableOutput("deseq_table")
                ),
                tabPanel(
                    value = 6,
                    title = "Gene Boxplot",
                    br(),
                    plotOutput("deseq_box", width = "100%")
                ),
                tabPanel(
                    value = 7,
                    title = "Gene Clustering",
                    br(),
                    plotOutput("deseq_cluster", width = "100%")
                ),
                tabPanel(
                    value = 8,
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
