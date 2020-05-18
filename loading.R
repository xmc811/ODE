

packages <- c("shiny", 
              "shinythemes", 
              "shinyWidgets",
              "tidyverse", "magrittr","DT", 
              "RColorBrewer", "circlize", "ComplexHeatmap", "scales",
              "fgsea","maftools")

lapply(packages, require, character.only = TRUE)



