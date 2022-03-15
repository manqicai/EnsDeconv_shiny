# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#options(shiny.maxRequestSize=50000*1024^2)

#library(BiocManager)
#options(repos = BiocManager::repositories())
########Library Packages
library(shiny)
library(markdown)
library(knitr)
library(shinythemes)
library(shinycustomloader)
library(shinyhelper)
library(foreach)
library(tidyverse)
library(EnsDeconv)
library(sparseMatrixStats)
#library(scran)
library(ggpubr)



# load(file = "data/human_lengths.rda")
# load(file = "data/db.RData")

#source file
# source("methods/Deconv_meth_helper.R")
# source("methods/my_bisque.R")
# source("methods/my_music_construct.R")
# source("methods/my_music_utils.R")
# source("methods/cosg.R")
# source("pkg/DeconRNASeq.R")
# source("pkg/hspe.R")
shinyServer(function(input, output,session){
   
########### Update reference ###########
   refdata <- reactive({
      req(input$metaref)
      if(str_detect(input$metaref$datapath,"csv")){
         metaref <- read.csv(input$metaref$datapath,
                              header = TRUE)
      }else if(str_detect(input$metaref$datapath,"rds")){
         metaref <- readRDS(input$metaref$datapath)
      }else if(str_detect(input$metaref$datapath,"RData")){
         metaref <- load(input$metaref$datapath)
      }
   })
   
   
   observeEvent(refdata(), {
      updateSelectInput(session, "columnsref", choices=colnames(refdata()))
      updateSelectInput(session, "columnssample", choices=colnames(refdata()))
   })


###############CRM output################
#########################################
  observe_helpers()
   dcInput <- eventReactive(input$dcupdate, {
     nct <- input$numofct

     req(input$bulk)
     #load bulk data
     if(str_detect(input$bulk$datapath,"csv")){
        to_deconv <- read.csv(input$bulk$datapath,
                              header = TRUE)
     }else if(str_detect(input$bulk$datapath,"rds")){
        to_deconv <-readRDS(input$bulk$datapath)
     }else if(str_detect(input$bulk$datapath,"RData")){
        load(input$bulk$datapath)
     }
     # load meta ref
     if(str_detect(input$metaref$datapath,"csv")){
        metaref <- read.csv(input$metaref$datapath,
                              header = TRUE)
     }else if(str_detect(input$metaref$datapath,"rds")){
        metaref <- readRDS(input$metaref$datapath)
     }else if(str_detect(input$metaref$datapath,"RData")){
        metaref <- load(input$metaref$datapath)
     }
     # load reference data
     if(str_detect(input$ref$datapath,"csv")){
        ref <- read.csv(input$ref$datapath,
                            header = TRUE)
     }else if(str_detect(input$ref$datapath,"rds")){
        ref <- readRDS(input$ref$datapath)
     }else if(str_detect(input$ref$datapath,"RData")){
        ref <- load(input$ref$datapath)
     }

     #gene <- intersect(rownames(to_deconv), rownames(ref))
     #to_deconv <- to_deconv[gene,]
     #ref <- ref[pmatch(gene,rownames(ref)),]
     #to_deconv <- as.matrix(to_deconv)
     #ref<- as.matrix(ref)
     
     # scaling
     # to_deconv <- Normalization(to_deconv,input$norm)
     # ref <- Normalization(ref,input$norm)
     # 
     # 
     # # Scaling
     # if(input$scale == "log"){
     #    to_deconv = log2(to_deconv+1)
     #    if(class(ref)[[1]] == "dgCMatrix"){
     #       ref@x <- log2(ref@x + 1)
     #    }else{
     #       ref <- log2(ref+1)
     #    }
     # }
      
     bat <- ifelse(input$batchcorrec == "Yes",TRUE,FALSE)
     
     
     ref_list = list()
     ref_list$"ref" = list()
     ref_list$"ref"$ref_matrix = as.matrix(ref)
     ref_list$"ref"$meta_ref = metaref
     
     res <- list()
     #res[[paste0(mth,"_",mrk_method,"_",input$scale,"_",input$norm)]] <- run_deconv_method(method_name = mth, to_deconv = as.matrix(to_deconv),ref_matrix= ref,meta_ref = metaref, markers,data_type = datatype,marker_method = mrk_method,batchcorrec = bat,scale = input$scale,nmrks = input$nmrk,deconv_clust = input$columnsref,samplesname = input$columnssample)
     params = get_params(data_type = input$datatype, data_name = "ref", n_markers = input$nmrk,Marker.Method = input$mrk,TNormalization = input$norm,CNormalization =input$norm ,dmeths = input$Deconv,Scale = input$scale)
     #params
     res <- EnsDeconv(count_bulk = as.matrix(to_deconv), ref_list = ref_list, ncore = 4, parallel_comp = F, params = params)
     res_p = lapply(res[[2]], function(x) x[["a"]][["p_hat"]][[1]][[1]])
     res_p[["Ensemble"]] = res[["EnsDeconv"]][["ensemble_p"]]

     res_new <- list(res_p = res_p, case_ctrl = input$case_ctrl)
     #outt = lapply(res, function(x) x[["estimate"]])
     res_p[sapply(res_p, is.null)] <- NULL
     res_p
     })


   # output$dctable <- renderTable({
   #    dcInput()}, rownames = TRUE,colnames = T)
   
   # output$dcplots <- renderPlot({
   #    res<- dcInput()
   #    
   #    
   #    # Basic piechart
   #    # ggplot(res[[1]], aes(x="", y=prop, fill=group)) +
   #    #   geom_bar(stat="identity", width=1, color="white") +
   #    #   coord_polar("y", start=0) +
   #    #   theme_void() + 
   #    #   theme(legend.position="none") +
   #    #   
   #    #   geom_text(aes(y = ypos, label = group), color = "white", size=6) +
   #    #   scale_fill_brewer(palette="Set1")
   #    
   #    pheatmap(res[[1]], display_numbers = T)
   #    
   #   })
   
   output$plots <- renderPlot({
    
     res<- dcInput()
     n <- length(res)
     print(res)
     plot_output_list <- lapply(1:n, function(i) {

       plotname <- names(res)[[i]]
       plotOutput(plotname, height = 580, width = 550)


     })


     do.call(tagList, plot_output_list)
     df <- list()
     for (i in 1:length(res)) {
        u <- intersect(rownames(res[[1]]),rownames(res[[i]]))
        d <- data.frame(Subject = u,Method = names(res)[i],res[[i]][u,]/rowSums(res[[i]][u,]))
        df[[i]] <-reshape2::melt(d,id.vars = c('Subject','Method'))
     }
     df <- bind_rows(df)
     names(df) <- c("Subject","Method","CellType","p_hat")

     ggboxplot(df,x = "CellType", y ="p_hat",color = "CellType")+facet_wrap(~Method)+ylab("Estimated cell type proportion")+xlab("Cell type")+color_palette("jco")


   })
   # observe({             
   #   
   #   for (i in 1:length(dcInput())) {
   #     local({ 
   #       
   #       
   #       plotname <- paste("plot", i, sep="")
   #       my_i <- i
   #       plotname <- names(dcInput())[[my_i]]
   #       output[[plotname]] <- renderPlot({
   #         
   #         #function_plot is the function generate plot
   #         #pheatmap( dcInput()[[my_i]], display_numbers = T,main = plotname)
   #          boxplot(dcInput()[[my_i]], main = plotname)
   #         
   #       })
   #     })#endlocal
   #   }
   #   
   # })

   # Downloadable csv of selected dataset ----
   output$downloadData <- downloadHandler(
     filename = function() {
       "outcome.csv"
     },
     content = function(file) {
        resout = dcInput()
        resout = lapply(resout, as.data.frame)
        resout <- lapply(resout, function(x) rownames_to_column(x, "sample"))
        write_xlsx(resout, file)
     }
   )
   
   data <- reactive({
      req(input$metabulk)
      if(str_detect(input$metabulk$datapath,"csv")){
         metabulk <- read.csv(input$metabulk$datapath,
                             header = TRUE)
      }else if(str_detect(input$metabulk$datapath,"rds")){
         metabulk <- readRDS(input$metabulk$datapath)
      }else if(str_detect(input$metabulk$datapath,"RData")){
         metabulk <- load(input$metabulk$datapath)
      }
   })
   
   # filtereddata <- eventReactive({
   #    input$update
   #    data()
   # },  {
   #    req(data())
   #    if(is.null(input$columns) || input$columns == "")
   #       data() else 
   #          data()[, colnames(data()) %in% input$columns]
   # })
   
   observeEvent(data(), {
      updateSelectInput(session, "columns", choices=colnames(data()))
   })

})
