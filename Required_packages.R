#devtools::install_github("randel/EnsDeconv", dependencies = TRUE)
#BiocManager::install("limma", character.only = TRUE)
#BiocManager::install("AnnotationDbi", character.only = TRUE)
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

pkgs <- c("shiny",
          "ggplot2",
          "dplyr",
          "tidyverse",
          "devtools",
          "shinythemes",
          "shinycustomloader",
          "shinyhelper",
          "foreach",
          "sparseMatrixStats",
          "ggpubr",
          "DeconRNASeq",
          "formattable",
          "shinysky",
          "shinycssloaders",
          "shinyWidgets",
          "knitr",
          "markdown",
          "rintrojs"
          )
check.packages(pkgs)

required_packages <- c(
  "EnsDeconv"
)

new.packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(new.packages)) {
  devtools::install_github("randel/EnsDeconv", dependencies = TRUE)
}

# library(reticulate)
# py_install("statistics")
# py_install("sys")
# library(knitr)
# library(markdown)
# library(shiny)
# library(shinythemes)
# library(shinycustomloader)
# library(shinyhelper)
# library(foreach)
# library(tidyverse)
# library(EnsDeconv)
# library(sparseMatrixStats)
# #library(scran)
# library(ggpubr)
# library(DeconRNASeq)
# library(dplyr)
# library(formattable)
# library(shinysky)
# library(shinycssloaders)
# library(shinyWidgets)
#library(rintrojs)