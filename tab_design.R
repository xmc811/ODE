
dnaMaf <- "./data/DNA_filtered__mutect_Pindel_GM_pipe_201002.maf"
rnaTsv <- "./data/RNA_merged_HTseq_count.tsv"
rppa <- "./data/IACS__RPPA.txt"

hmks_hs <- fgsea::gmtPathways("./data/h.all.v7.0.symbols.gmt")


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
            h3("RNA-seq Data Analysis"),
            br(),
            
            h4("Data Source"),
            splitLayout(radioGroupButtons(inputId = "data_source",
                                          label = NULL,
                                          choices = c("Example","Upload"),
                                          justified = TRUE),
                        actionButton(
                            inputId = "rna_start",
                            label = "Launch",
                            icon = icon("bar-chart"),
                            style = "color: white; background-color: #0570b0;
                            float:right; margin-right: 5px;"),
                        
                        cellWidths = c("67%", "33%")),
            
            conditionalPanel(
                condition = "input.data_source == 'Upload'",
                fileInput(inputId = "rna_input",
                          label = NULL,
                          buttonLabel = "Browse..")
                ),
            br(),
            
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
                                         value = 1)),
                br()
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 6 || 
                            input.rna_panel == 7",
                
                h4("Gene List"),
                radioGroupButtons(inputId = "rna_gene_ls_src",
                                  label = NULL,
                                  choices = c("Use Top Genes",
                                              "Manual Input",
                                              "Upload File"),
                                  justified = TRUE),
                conditionalPanel(
                    condition = "input.rna_gene_ls_src == 'Use Top Genes'",
                    splitLayout(sliderInput(inputId = "rna_gene_num",
                                            label = "Number of Genes", 
                                            min = 1, max = 50, value = 6),
                                actionButton(
                                    inputId = "rna_gene_read1",
                                    label = "Plot",
                                    icon = icon("check"),
                                    style = "color: white; 
                                    background-color: #737373;
                                    margin-top: 25px;
                                    float:right;
                                    margin-right: 5px;"),
                                cellWidths = c("75%", "25%")
                    )
                ),
                conditionalPanel(
                    condition = "input.rna_gene_ls_src == 'Manual Input'",
                    splitLayout(textInput("rna_genes_man", 
                                          label = NULL, 
                                          value = ""),
                                
                                actionButton(
                                    inputId = "rna_gene_read2",
                                    label = "Plot",
                                    icon = icon("check"),
                                    style = "color: white; background-color: #737373;
                            float:right; margin-right: 5px;"),
                                cellWidths = c("75%", "25%")
                                )
                    
                ),
                conditionalPanel(
                    condition = "input.rna_gene_ls_src == 'Upload File'",
                    fileInput(inputId = "rna_genes_file",
                              label = NULL,
                              buttonLabel = "Browse..")
                ),
                br()
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 8",
                
                h4("Pathway List"),
                splitLayout(radioGroupButtons(inputId = "rna_pathway_src",
                                              label = NULL,
                                              choices = c("Use Hallmarks",
                                                          "Upload File"),
                                              justified = TRUE),
                            cellWidths = "67%"),
                
                conditionalPanel(
                    condition = "input.rna_pathway_src == 'Upload File'",
                    fileInput(inputId = "rna_pathway_file",
                              label = NULL,
                              buttonLabel = "Browse..")
                )
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 1 || 
                            input.rna_panel == 2 ||
                            input.rna_panel == 6 || 
                            input.rna_panel == 7",
                
                h4("Plotting Parameters"),
                splitLayout(selectInput(inputId = "palette_cat", 
                                        label = "Categorical Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                        selected = "Set2"),
                            selectInput(inputId = "palette_con", 
                                        label = "Continuous Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                        selected = "Spectral")),
                materialSwitch(
                    inputId = "palette_dir",
                    label = "Reverse Scale Color Direction",
                    value = FALSE,
                    right = TRUE
                )
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4",
                
                h4("Plotting Parameters"),
                splitLayout(numericInput("p_plot_lim", 
                                         label = "Adjusted P-value Squash", 
                                         value = 5),
                            numericInput("lfc_plot_lim", 
                                         label = "Log2 Fold Change Squash", 
                                         value = 5)),
            ),
            
            conditionalPanel(
                condition = "input.rna_panel != 5",
                
                splitLayout(numericInput("rna_plot_height", 
                                         "Plot Height (px)", 
                                         value = 600),
                            numericInput("rna_plot_width", 
                                         "Plot Width (px)", 
                                         value = 800),
                            cellWidths = c("50%", "50%"))
            ),
            
            tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              .progress-bar {
                              height: 20px;
                              }
                              .shiny-notification {
                              width: 200px;
                              top: 50%;
                              right: 30%;
                              position: fixed;
                              }
                              ")))
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
                    plotOutput("deseq_ma", width = "100%") 
                ),
                tabPanel(
                    value = 4,
                    title = "Volcano Plot",
                    br(),
                    plotOutput("deseq_volcano")
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
