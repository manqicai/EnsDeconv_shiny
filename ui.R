#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#devtools::install_github("randel/EnsDeconv", dependencies = TRUE)
########Library Packages#

#options(repos = BiocManager::repositories())
#source file
# source("methods/Deconv_meth_helper.R")
# source("methods/my_bisque.R")
# source("methods/my_music_construct.R")
# source("methods/my_music_utils.R")

library(knitr)
library(markdown)
library(shiny)
library(shinythemes)
library(shinycustomloader)
library(shinyhelper)
library(foreach)
library(tidyverse)
library(EnsDeconv)
library(sparseMatrixStats)
#library(scran)
library(ggpubr)
library(DeconRNASeq)
library(dplyr)
library(formattable)
library(shinysky)
library(shinycssloaders)
library(shinyWidgets)




appCSS <- "
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"
# parameter_tabs ----
parameter_tabs <- tabsetPanel(
  id = "params",
  type = "hidden",
  tabPanel("Select..."),
  tabPanel("custom",
           fluidRow(column(width = 6,fileInput("ref", label = h4("Ref Data"),
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv",
                                                          ".rds",
                                                          ".txt"))%>%
                             helper(colour = "green",type = "inline",size = "m",content = "Upload the file of reference data (.csv, .rds or .txt),rows are genes and columns are samples")),
                    column(width = 6,fileInput("metaref", label = h4("Meta data"),
                                               accept = c("text/csv",
                                                          "text/comma-separated-values,text/plain",
                                                          ".csv",
                                                          ".rds",
                                                          ".txt"))%>%
                             helper(colour = "green",type = "inline",size = "m",content = "Upload the file of meta data for the reference data (.csv, .rds or .txt)"))),
           fluidRow(column(width = 6, selectInput("columnsref", h4("Select cell type variable"), choices = NULL)%>% helper(colour = "green",type = "inline", content = "Select the variable that correspond to cell type cluster")),
                     column(width = 6,selectInput("columnssample", h4("Select sample ID variable"), choices = NULL)%>% helper(colour = "green",type = "inline", content = "Select the variable that correspond to sample ID")))
  ),
  tabPanel("brain", 
           selectInput("localbrain", label = h4("Choose brain reference data"),choices = list("Darmanis","Habib"))%>%
             helper(colour = "green",type = "inline",size = "m",content = "Choose brain reference data")
  ),
  tabPanel("blood",
           selectInput("localblood", label = h4("Choose blood reference data"),choices = list("lm22","skin_signature"))%>%
             helper(colour = "green",type = "inline",size = "m",content = "Choose blood reference data")
  )
)


parallel_parameter_tabs <- tabsetPanel(
  id = "para_params",
  type = "hidden",
  tabPanel("FALSE"),
  tabPanel("TRUE",
           fluidRow(column(width = 12,numericInput("ncore", label = h4("Num. of cores"),value = 2)
                                    %>% helper(colour = "green",type = "inline", content = "The number of cores to use for parallel execution"))
  )
))


navbarPage(title = div("EnsDeconv (Ensemble Deconvolution)",tagList(a(href = "https://publichealth.pitt.edu/biostatistics/", style = "color:#606060", tags$img(src='University_of_Pittsburgh_Logo_CMYK_Primary_3-Color.png',height = 30,width =60)))),

          theme = shinytheme("yeti"),
           fluid = TRUE,

           tabPanel(value = "int",title = "Introduction" ,
                    fluidRow(align="center",
			style="opacity: 1;background-color:white; margin-top: 0px;width: 100%;",
			helpText(strong("- Introduction -" , style="color:green ; font-family: 'times'; font-size:30pt ; font-type:bold" ))),
                    	fluidRow(
			style="opacity: 1;background-color:white; margin-top: 0px;width: 100%;",

			column(6,offset=3,
				includeMarkdown("RMD/intro.Rmd")
				),
				br()
				)
			,


                    tags$head(tags$script(HTML('
          var customHref = function(tabName) {
                                               var dropdownList = document.getElementsByTagName("a");
                                               for (var i = 0; i < dropdownList.length; i++) {
                                               var link = dropdownList[i];
                                               if(link.getAttribute("data-value") == tabName) {
                                               link.click();
                                               }
                                               }
                                               };
                                               '))),

                    fluidRow(align="center",
			style="opacity: 1;background-color:white; margin-top: 0px;width: 100%;",
			column(6,offset=3,
				# Set the style of this page
				#br(),
				helpText(strong("- Get Started -" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ) ,
				helpText(HTML("<a onclick=","customHref('gs')" ,">",
				              "Get Started","</a>")),
				hr(),
				helpText(strong("- Analysis -" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ) ,
				helpText(HTML("<a onclick=","customHref('dc')" ,">",
				              "Analysis","</a>")),
				# helpText(HTML("<a onclick=","customHref('braindc')" ,">",
				#               "Brain Data","</a>")),
				# helpText(HTML("<a onclick=","customHref('blooddc')" ,">",
				#               "Blood Data","</a>")),
	
				hr(),
				p("We summarized following methods:"),
				helpText(strong("- Deconvolution Method -" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ) ,
				helpText(HTML("<a onclick=","customHref('prf')" ,">",
                        "Partial Reference Free Method","</a>")),
				helpText(HTML("<a onclick=","customHref('rb')" ,">",
				              "Reference Based Method","</a>"))
	
			))

                    ),
			tabPanel("Get Started",
			         value = "gs",
			         helpText(strong("Get Started" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ),
			         includeMarkdown("RMD/gs.Rmd")
			),
# General Analysis ----
tabPanel("Analysis",value = "dc",
			                    helpText(strong("Analysis" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ),
			                    # sidebarPanel(
			                    #   selectInput("chooseref",h4("References"), 
			                    #               choices = c("Custom", "local brain references", "local blood references")
			                    #   ),
			                    #   parameter_tabs,
			                    # ),
			                    sidebarPanel(strong("Required"),
			                                 fluidRow(column(width = 6,fileInput("bulk", label = h4("Bulk Data"),
			                                                                     accept = c("text/csv",
			                                                                                "text/comma-separated-values,text/plain",
			                                                                                ".csv",
			                                                                                ".rds",
			                                                                                ".RData",
			                                                                                ".txt"))
			                                                 %>%
			                                                   helper(colour = "green",type = "inline", content = "Upload the file of bulk data you want to deconvolve (.csv, .rds or .txt),rows are genes and columns are samples"))),
			                                 
			                                 fluidRow(column(width = 12,
			                                                   selectInput("chooseref",h4("References"),
			                                                               choices = c("Select...","Custom" = "custom", "brain references"= "brain", "blood references" = "blood")
			                                                   )%>%
			                                                   helper(colour = "green",type = "inline", content = "Custom: upload customized references dataset; brain (blood) references: select existing references"),
			                                                   parameter_tabs))
			                                 #actionButton("updateref", "incorporate external information", class = "btn-info"),
			                                ,fluidRow(column(width = 6, multiInput("Deconv", label = h4("Deconvolution Method"),
			                                                                       choices = list("dtangle", "hspe","DSA","CIBERSORT","EPIC","MuSiC","Bisque", "ICeDT","DeconRNASeq","FARDEEP","DCQ"),
			                                                                       selected = c("Bisque","CIBERSORT"))
			                                                 %>% helper(colour = "green",type = "inline", content = "Select the deconvolution methods that you want to apply")),
			                                          column(width = 6,  multiInput("mrk", label = h4("Marker Gene Approach"),
			                                                                        choices = list("none"  , "t","wilcox","combined","p.value","regression"),
			                                                                        selected = "none")
			                                                 %>% helper(colour = "green",type = "inline", content = "Choose the marker gene selection methods that you want to apply"))),
			                                 
			                                
			                               
			                                fluidRow(column(width = 6,numericInput("nmrk", label = h4("Num. of Marker Gene"),value = 50)
			                                                %>% helper(colour = "green",type = "inline", content = "Enter the number of markers")),
			                                         column(width = 6,selectInput("datatype", label = h4("Type of reference data"),choices = list("singlecell-rna","microarray"))%>%
			                                                  helper(colour = "green",type = "inline",size = "m",content = "Choose the type of reference data"))),
			                                fluidRow(column(width = 6,			                                 checkboxGroupButtons("scale", label = h4("Type of scaling approach"),choices = c("log","linear"),
			                                                                                                                   status = "primary",
			                                                                                                                   checkIcon = list(
			                                                                                                                     yes = icon("ok", 
			                                                                                                                                lib = "glyphicon"),
			                                                                                                                     no = icon("remove",
			                                                                                                                               lib = "glyphicon")),
			                                                                                                                   selected = "linear")%>%
			                                                  helper(colour = "green",type = "inline",size = "m",content = "Choose the scaling approach of bulk & reference data")),
			                                         column(width = 6,  checkboxGroupButtons("norm", label = h4("Type of normalization approach"),choices = c("CPM","TPM","QN","none"),
			                                                                                 status = "primary",
			                                                                                 checkIcon = list(
			                                                                                   yes = icon("ok", 
			                                                                                              lib = "glyphicon"),
			                                                                                   no = icon("remove",
			                                                                                             lib = "glyphicon")),
			                                                                                 selected = "none")%>%
			                                                  helper(colour = "green",type = "inline",size = "m",content = "Choose the Normalization approach of bulk & reference data")))
			                                 ,fluidRow(column(width = 12,
			                                                  selectInput("choosepara",h4("Parallel Computing"),
			                                                              choices = c("False"= "FALSE","True" = "TRUE")
			                                                  )%>%
			                                                    helper(colour = "green",type = "inline", content = "Use parallel computing or not"),
			                                                  parallel_parameter_tabs))
			                                #actionButton("updateref", "incorporate external information", class = "btn-info"),
			                                ,
			                                 actionButton("dcupdate", "Run", class = "btn-info"),
			                                 downloadButton("downloadData", "Download results")
			                                 # ,hr(),
			                                 # strong("Selected"),
			                                 # switchInput(
			                                 #   inputId = "case_ctrl",
			                                 #   label = "Compare bulk data by case and control",
			                                 #   onLabel = "Yes",
			                                 #   offLabel = "No")%>%
			                                 #   helper(colour = "green",type = "inline",size = "m",content = "Choose whether plot boxplot based on case and control"),
			                                 # fileInput("metabulk", label = h4("Meta Data for Bulk data"),
			                                 #           accept = c("text/csv",
			                                 #                      "text/comma-separated-values,text/plain",
			                                 #                      ".csv",
			                                 #                      ".rds",
			                                 #                      ".RData"))
			                                 # %>%
			                                 #   helper(colour = "green",type = "inline", content = "Upload the file of bulk data you want to deconvolve (.csv, .rds or .txt)"),
			                                 # actionButton("update", "incorporate external information", class = "btn-info"),
			                                 # selectInput("columns", h4("Select case control variable"), choices = NULL) %>%
			                                 #   helper(colour = "green",type = "inline", content = "Select the case control indicator variable"),
			                                 # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
			                                 #                  tags$div("Loading...",id="loadmessage"))
			                                ),
			                    mainPanel(strong("Summary "),
			                              p("DeconvolutionMethods_MarkerGeneSelection_Scale_Normalization: "),
			                              #tableOutput("dctable")%>% withSpinner(),
			                              # uiOutput("plots") %>% withSpinner(),
			                              plotOutput("plots") %>% withSpinner(),
			                              
			                              p(""),
			                              withMathJax()
			                              #,
			                              #strong("References")
			                              
			                    )
			           )
	# # Brain Data ----
	# 		           tabPanel("Brain Data",
	# 		                    value = "braindc",
	# 		                    helpText(strong("Brain Data" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ),
	# 		                    # sidebarPanel(
	# 		                    #   selectInput("chooseref",h4("References"), 
	# 		                    #               choices = c("Custom", "local brain references", "local blood references")
	# 		                    #   ),
	# 		                    #   parameter_tabs,
	# 		                    # ),
	# 		                    strong("References"),
	# 		           ),
	# 		           tabPanel("Blood Data",
	# 		                    value = "blooddc",
	# 		                    helpText(strong("Blood Data" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ),
	# 		                    strong("References")
	# 		           )
			,navbarMenu("Deconvolution Method",
                   tabPanel("Partial Reference Free Method (Marker guided)",
			                    value = "prf",
			                    helpText(strong("Partial Reference Free Method" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ),
			                    includeMarkdown("RMD/prf.Rmd"),
			                    strong("References")
			           ),
			           tabPanel("Reference Based Method",
			                    value = "rb",
			                    helpText(strong("Reference Based Method" , style="color:green ; font-family: 'times'; font-size:20pt ; font-type:bold" ) ),
			                    includeMarkdown("RMD/rb.Rmd"),
			                    #includeMarkdown("Tabls.Rmd"),
			                    strong("References")

			           )
			),
          tabPanel("About",
                     fluidRow(align="center",
                              style="opacity:0.9; background-color: white ;margin-top: 0px; width: 100%;",
                              column(6,offset=3,
                                     br(),
                                     helpText( strong("About", style="color:Green ; font-family: 'times'; font-size:30pt ; font-type:bold" )) ,
                                     hr()
                              )),
                   fluidRow(
                     style="opacity:0.9; background-color: white ;margin-top: 0px; width: 100%;",
                     column(6,offset=3, align="center",
                            img(src="University_of_Pittsburgh_Logo_CMYK_Primary_3-Color.png" ,  width = 300),
                            br(),
                            br(),
                            #p("ShinyApp created by ", strong("Manqi Cai")," of "),
                            p(strong("University of Pittsburgh")),
                            p(("Copyright (C) 2021, code licensed "))
                     ))
                     )

)
