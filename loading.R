
packages <- c("shiny", "shinythemes", "tidyverse", "magrittr", "DESeq2","DT", 
              "RColorBrewer", "circlize", "ComplexHeatmap", "scales",
              "fgsea", "maftools")

lapply(packages, require, character.only = TRUE)

