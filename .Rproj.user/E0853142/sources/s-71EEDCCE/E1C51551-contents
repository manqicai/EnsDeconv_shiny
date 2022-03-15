library(reticulate)
py_install("statistics")
py_install("sys")

# packages <- c("devtools", "BiocManager","data.table","ggplot2","tidyverse",
#               "Matrix","matrixStats",
#               "gtools",
#               "foreach","doMC","doSNOW", #for parallelism
#               "Seurat","sctransform", #sc-specific normalization
#               "nnls","FARDEEP","MASS","glmnet","ComICS","dtangle") #bulk deconvolution methods
# 
# for (i in packages){ install.packages(i, character.only = TRUE)}
# 
# # Installation using BiocManager:
# # Some packages that didn't work with install.packages (e.g. may not be present in a CRAN repository chosen by the user)
# packages3 = c('limma','edgeR','DESeq2','pcaMethods','BiocParallel','preprocessCore','scater','SingleCellExperiment','Linnorm','DeconRNASeq','multtest','GSEABase','annotate','genefilter','preprocessCore','graph','MAST','Biobase',"AnnotationDbi") #last two are required by DWLS and MuSiC, respectively.
# for (i in packages3){ BiocManager::install(i, character.only = TRUE)}
# 
# # Dependencies for CellMix: 'NMF', 'csSAM', 'GSEABase', 'annotate', 'genefilter', 'preprocessCore', 'limSolve', 'corpcor', 'graph', 'BiocInstaller'
# packages2 = c('NMF','csSAM','limSolve','corpcor')
# for (i in packages2){ install.packages(i, character.only = TRUE)}
# # 
# # # Special instructions for CellMix and DSA
#  install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.7/bioc/")
# system('wget http://web.cbio.uct.ac.za/~renaud/CRAN/src/contrib/CellMix_1.6.2.tar.gz')
# system("R CMD INSTALL CellMix_1.6.2.tar.gz")
# system('wget https://github.com/zhandong/DSA/raw/master/Package/version_1.0/DSA_1.0.tar.gz')
# system("R CMD INSTALL DSA_1.0.tar.gz")
# 
# # Following packages come from Github
# devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE) #requires knitr
# devtools::install_github("xuranw/MuSiC") 
# devtools::install_bitbucket("yuanlab/dwls", ref="default")
# devtools::install_github("meichendong/SCDC")
# devtools::install_github("rosedu1/deconvSeq")
# devtools::install_github("cozygene/bisque")
# devtools::install_github("dviraran/SingleR@v1.0")
devtools::install_github("randel/EnsDeconv", dependencies = TRUE)
#BiocManager::install("limma", character.only = TRUE)
#BiocManager::install("AnnotationDbi", character.only = TRUE)
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

pkgs <- c("shiny",
          "DT",
          "readxl",
          "ggplot2",
          "ggrepel",
          "vegan",
          "RColorBrewer",
          "data.table",
          "stringr",
          "dplyr",
          "scales",
          "plotly",
          "shinyjs",
          "tidyverse",
          "devtools")
check.packages(pkgs)