
tab_upload <- tabPanel(

    title = "Upload",
    fluid = TRUE,
    value = "v_up",

    sidebarLayout(

        sidebarPanel = sidebarPanel(

            h4("Please upload your data:"),
            fileInput(inputId = "data",
                      label = "Data",
                      buttonLabel = "Browers...",
                      placeholder = rnaRds),
            br(),
            h3("Data sets used:"),
            textInput("rdsfile", "RDS data file name", rnaRds),
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
            h4("Specify genes and samples:"),
            uiOutput("specifyRnaGenes"),
            uiOutput("specifyRnaSamples"),
            actionButton("resetRnaGenesSamples", "Reset"),
            br(),
            h4("RNA-seq significance cutoff:"),
            numericInput("p_co", label = "Adjusted p-value Cutoff", value = 0.05),
            numericInput("lfc_co", label = "Log2 Fold Change Cutoff", value = 1)
        ),

        mainPanel(
            tabsetPanel(
                tabPanel(
                    title = "Clustering",
                    br(),
                    plotOutput("clusterRna")
                ),
                tabPanel(
                    title = "PCA Plot",
                    br(),
                    h4("Under Construction ..."),
                ),
                tabPanel(
                    title = "Volcano Plot",
                    br(),
                    plotOutput("volcano", width = "100%")
                ),
                tabPanel(
                    title = "Table",
                    br(),
                    DT::dataTableOutput("table")
                ),
                tabPanel(
                    title = "Functional Enrichment",
                    br(),
                    h4("Under Construction ...")
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

makeMetaTable <- function(dnaMaf, rnaTsv, rppa){
    metaData <- NULL
    sid <- NULL
    if(!missing(dnaMaf)){
        mat <- read.delim(dnaMaf, sep = "\t", header = TRUE, as.is = TRUE)
        sid <- getSampleIds(mat, what = "DNA")
        metaData <- rbind(metaData, c(type = "DNA", file = dnaMaf))
    }
    if(!missing(rnaTsv)){
        mat <- read.delim(rnaTsv, sep = "\t", header = TRUE, as.is = TRUE)
        sid <- getSampleIds(mat, what = "RNA")
        metaData <- rbind(metaData, c(type = "RNA", file = rnaTsv))
    }
    if(!missing(rppa)){
        mat <- read.delim(rppa, sep = "\t", header = TRUE, as.is = TRUE)
        sid <- getSampleIds(mat, what = "RPPA")
        metaData <- rbind(metaData, c(type = "RPPA", file = rppa))
    }
    if(is.null(sid)){
        stop("DAN or RNA data set is required")
    }

    return(list(sid = sid, meta = metaData))
}

getSampleIds <- function(dataMat, what = "DNA"){
    if(toupper(what) == "DNA"){
        return(unique(dataMat[, "Tumor_Sample_Barcode"]))
    }
    if(toupper(what == "RNA")){
        return(colnames(dataMat)[-c(1, 2)])
    }
    if(toupper(what == "RPPA")){
        return(dataMat[, "Sample.Name"])
    }
    stop("DAN or RNA data set is required")

    return(invisible())
}
