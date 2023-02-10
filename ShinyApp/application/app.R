library(shiny)
# Library to read genes.gtf
library(GenomicFeatures)
# Library to implement Random Forest Classifier
library(randomForest)
# Library to implement XGBoost Classifier
library(xgboost)
# Libraries to prepare data and evaluate results
library(ROCR)
library(MLmetrics)
library(caret)
library(magrittr)
library(dplyr)
library(Matrix)
library(pracma)
library(data.table)
# Source functions
source("./sleuth_prep.R")
source("./cufflinks_prep.R")
source("./hisat_prep.R")
source("./bitseq_prep.R")
source("./rsem_prep.R")
source("./ebseq_prep.R")
source("./simdata_prep.R")
source("./randomForestModel.R")
source("./xGBoostModel.R")
source("./combineResults.R")
source("./xgbPrediction.R")
source("./rfPrediction.R")

# genes.gtf
txdb <- makeTxDbFromGFF(file = "./data/genes.gtf", format = "gtf")
genes <- as.data.frame(transcripts(txdb, columns = c("GENEID", "TXNAME"))@elementMetadata)
colnames(genes) <- c("gene_id", "transcript_id")

# Sleuth Data
sleuth_data <- read.csv(file = "./data/simulated_data.transcript_results.sleuth.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sleuth_data <- sleuth_prep(sleuth_data, genes)

# Tophat2 - Cufflinks Data
cufflinksIsoformsDiff <- read.table(file = "./data/isoform_exp.diff", sep = "", header = TRUE, stringsAsFactors = FALSE)
cufflinksIsoformsFpkm <- read.table(file = "./data/isoforms.fpkm_tracking", sep = "", header = TRUE, stringsAsFactors = FALSE)
tophat_data <- cufflinks_prep(cufflinksIsoformsDiff, cufflinksIsoformsFpkm)

# Hisat2 - DESeq2 Data
hisat2_data <- read.table(file = "./data/simulated_data.transcript_results.Deseq2.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
hisat2_data <-  hisat_prep(hisat2_data, genes)

# Bitseq Data
bitseq_data <- read.table(file = "./data/data.pplr", sep = " ", header = FALSE, stringsAsFactors = FALSE)
bitseq_transcripts <- read.table(file = "./data/sample_01.bam.data.tr", sep = " ", header = FALSE, stringsAsFactors = FALSE)[2]
bitseq_data <- bitseq_prep(bitseq_data, bitseq_transcripts, genes)

# RSEM Data
rsem_data <- read.table(file = "./data/IsoMat.de.txt", fill = TRUE , sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rsem_data <- setDT(rsem_data, keep.rownames = 'transcript_id')[]
rsem_data <- rsem_prep(rsem_data, genes)

# EBSeq Data
ebseq_data <- read.table(file = "./data/EBSeq.Simulated.All.Out.FC.txt", fill = TRUE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ebseq_data <- setDT(ebseq_data, keep.rownames = TRUE)
colnames(ebseq_data) <- c("trx", "PostFC")

ebseq_data2 <- read.table(file = "./data/EBSeq.Simulated.DE.Out.txt", fill = TRUE , sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ebseq_data2 <- setDT(ebseq_data2, keep.rownames = TRUE)[]
colnames(ebseq_data2) = col.names = c("trx", "PPEE", "PPDE")
ebseq_data <- ebseq_prep(ebseq_data, ebseq_data2, genes)

# Simulated Data
simData <- read.table(file = "./data/sim_tx_info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
simData <- simdata_prep(simData, genes)

# Custom sidebarpanel
sidebarPanel2 <- function (..., out = NULL, width = 4) 
{
  div(class = paste0("col-sm-", width), 
      tags$form(class = "well", ...),
      out
  )
}

# UI side
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(50% - 150px);
             left: calc(50% - 150px);
             width: 300px;
             border: 1px solid black;
             }"
      ),
    )
  ),
  uiOutput(outputId = "trainUI"), # UI used to train model
  uiOutput(outputId = "predictUI") # UI used for DE status prediction
)

# Server side
server <- function(input, output, session) {
  query_modal <- modalDialog(
    title = tags$h3(tags$b("Select workflow to continue.")),
    selectInput(inputId = "workflow", label = "Workflow", choices = c("Prediction", "Model Training"), selected = NULL),
    size = "l", fade = FALSE,
    actionButton(inputId = "actionWork", label = "Submit", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
  # on application startup show pop-up window
  showModal(query_modal)
  # on pop-up submit button
  observeEvent(eventExpr = input$actionWork, handlerExpr = {
    # remove pop-up
    removeModal()
    m <- isolate(input$workflow)
    if(m == "Model Training"){
      output$trainUI <- renderUI({
        tabsetPanel(id = "inTab", type = "tabs",
                    # File Management Panel
                    tabPanel(
                      title = "File Management",
                      value = "selectFilePanel",
                      tags$br(),
                      sidebarLayout(
                        sidebarPanel(
                          tags$b("Select Pipelines for processing."),
                          tags$br(),
                          checkboxInput(inputId = "file1", label = tags$i("Kallisto - Sleuth Pipeline"), value = FALSE, width = "200px"),
                          checkboxInput(inputId = "file2", label = tags$i("Tophat2 - Cufflinks Pipeline"), value = FALSE, width = "200px"),
                          checkboxInput(inputId = "file3", label = tags$i("Hisat2 - DESeq2 Pipeline"), value = FALSE, width = "200px"),
                          checkboxInput(inputId = "file4", label = tags$i("BitSeq Pipeline"), value = FALSE, width = "200px"),
                          checkboxInput(inputId = "file5", label = tags$i("RSEM Pipeline"), value = FALSE, width = "200px"),
                          checkboxInput(inputId = "file6", label = tags$i("EBSeq Pipeline"), value = FALSE, width = "200px"),
                          tags$br(),
                          actionButton(inputId = "selected", label = "Submit", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                          , width = 3),
                        mainPanel(
                          tags$h2("Structured Results", align = "center"),
                          tabsetPanel(
                            tabPanel(title = "Kallisto - Sleuth",
                                     tags$br(),
                                     DT::dataTableOutput(outputId = "kallisto"),
                                     tags$br(),
                                     tags$hr(),
                                     verbatimTextOutput(outputId = "kallistoSummary")),
                            tabPanel(title = "Tophat2 - Cufflinks",
                                     tags$br(),
                                     DT::dataTableOutput(outputId = "tophat2"),
                                     tags$br(),
                                     tags$hr(),
                                     verbatimTextOutput(outputId = "tophat2Summary")),
                            tabPanel(title = "Hisat2 - DESeq2",
                                     tags$br(),
                                     DT::dataTableOutput(outputId = "hisat2"),
                                     tags$br(),
                                     tags$hr(),
                                     verbatimTextOutput(outputId = "hisat2Summary")),
                            tabPanel(title = "BitSeq",
                                     tags$br(),
                                     DT::dataTableOutput(outputId = "bitseq"),
                                     tags$br(),
                                     tags$hr(),
                                     verbatimTextOutput(outputId = "bitseqSummary")),
                            tabPanel(title = "RSEM",
                                     tags$br(),
                                     DT::dataTableOutput(outputId = "rsem"),
                                     tags$br(),
                                     tags$hr(),
                                     verbatimTextOutput(outputId = "rsemSummary")),
                            tabPanel(title = "EBSeq",
                                     tags$br(),
                                     DT::dataTableOutput(outputId = "ebseq"),
                                     tags$br(),
                                     tags$hr(),
                                     verbatimTextOutput(outputId = "ebseqSummary"))
                          )
                        )

                      ),
                      # Footer
                      tags$hr(),
                      div(actionButton(inputId = "next1", label = "Next", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          style="float:right"
                      )
                    ),
                    # Variable Selection Panel
                    tabPanel(
                      title = "Variable Selection",
                      value = "selectVariablePanel",
                      tags$br(),
                      tags$h4(tags$b("NOTE:"),"Pre-selected variables are strong predictors."),
                      tags$h5(tags$b("Id"),"should be included."),
                      tags$hr(),
                      fluidRow(
                        conditionalPanel(condition = "input.next1 > 0",
                                         column(width = 2,
                                                conditionalPanel(condition = "input.file1 == true",
                                                                 wellPanel(
                                                                   checkboxGroupInput(inputId = "kallistoCols",
                                                                                      label = "Sleuth Variables:",
                                                                                      choices = colnames(sleuth_data),
                                                                                      selected = c("Id", "qval", "pval", "b")
                                                                   )
                                                                 )
                                                )

                                         ),
                                         column(width = 2,
                                                conditionalPanel(condition = "input.file2 == true",
                                                                 wellPanel(
                                                                   checkboxGroupInput(inputId = "tophat2Cols",
                                                                                      label = "Tophat2 Variables:",
                                                                                      choices = colnames(tophat_data),
                                                                                      selected = c("Id", "p_value", "q_value", "log2.fold_change.")
                                                                   )
                                                                 )
                                                )
                                         ),
                                         column(width = 2,
                                                conditionalPanel(condition = "input.file3 == true",
                                                                 wellPanel(
                                                                   checkboxGroupInput(inputId = "hisat2Cols",
                                                                                      label = "Hisat2 Variables:",
                                                                                      choices = colnames(hisat2_data),
                                                                                      selected = c("Id", "pvalue", "padj", "log2FoldChange")
                                                                   )
                                                                 )
                                                )
                                         ),
                                         column(width = 2,
                                                conditionalPanel(condition = "input.file4 == true",
                                                                 wellPanel(
                                                                   checkboxGroupInput(inputId = "bitseqCols",
                                                                                      label = "BitSeq Variables:",
                                                                                      choices = colnames(bitseq_data),
                                                                                      selected = c("Id", "PPLR", "Log2FC")
                                                                   )
                                                                 )
                                                )
                                         ),
                                         column(width = 2,
                                                conditionalPanel(condition = "input.file5 == true",
                                                                 wellPanel(
                                                                   checkboxGroupInput(inputId = "rsemCols",
                                                                                      label = "RSEM Variables:",
                                                                                      choices = colnames(rsem_data),
                                                                                      selected = c("Id", "PPDE", "PostFC")
                                                                   )
                                                                 )

                                                )
                                         ),
                                         column(width = 2,
                                                conditionalPanel(condition = "input.file6 == true",
                                                                 wellPanel(
                                                                   checkboxGroupInput(inputId = "ebseqCols",
                                                                                      label = "EBSeq Variables:",
                                                                                      choices = colnames(ebseq_data),
                                                                                      selected = c("Id", "PPEE_", "PPDE_", "PostFC.IsoEBOut..PostFC")
                                                                   )
                                                                 )
                                                )
                                         )
                        )
                      ),
                      fluidRow(
                        div(actionButton(inputId = "back1", label = "Back", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                            style="float:left"
                            ),
                        div(actionButton(inputId = "next2", label = "Next", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                            style="float:right"
                            )
                      )
                    ),
                    # Random Forest and XGBoost Models
                    tabPanel(
                      title = "Model Training",
                      value = "modelPredictionPanel",
                      tags$br(),
                      sidebarLayout(
                        sidebarPanel(
                          tags$h3("Select Model"),
                          selectInput(inputId = "modelSelection",
                                      label = "",
                                      choices = c("Random Forest", "XGBoost"),
                                      selected = NULL,
                                      selectize = TRUE),
                          sliderInput(inputId = "numModels", label = "Number of models to train", value = 20, min = 1, max = 1000, step = 1),
                          actionButton(inputId = "submitModel", label = "Submit", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                          , width = 3),
                        mainPanel(
                          tags$h2("Results", align = "center"),
                          conditionalPanel(condition = "input.next2 > 0",
                                           htmlOutput(outputId = "modelMetrics"),
                                           tags$br(),
                                           verbatimTextOutput(outputId = "model"))
                          ),
                      ),
                      # Footer
                      tags$hr(),
                      div(actionButton(inputId = "back2", label = "Back", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          style="float:left"
                          ),
                      div(downloadButton(outputId = "downloadModel", label = "Download", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          style="float:right"
                          )
                    )
        )
      })
    }else if(m == "Prediction"){
      output$predictUI <- renderUI({
        tabsetPanel(id = "inTab2", type = "tabs",
                    tabPanel(
                      title = "Data Selection",
                      value = "modelInit",
                      tags$br(),
                      sidebarLayout(
                        sidebarPanel(
                          tags$br(),
                          fileInput(inputId = "newdata", label = "Upload Data File - Header needed."),
                          tags$hr(),
                          radioButtons(inputId = "sep", label = "Separator", choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = character(0)),
                          tags$br(),
                          actionButton(inputId = "check", label = "Submit", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                        ),
                        mainPanel(
                          DT::dataTableOutput(outputId = "userData"),
                          verbatimTextOutput(outputId = "userDataSummary")
                        )
                      ),
                      div(actionButton(inputId = "nextpan", label = "Next", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          style="float:right")
                      ),
                    tabPanel(
                      title = "Model Prediction",
                      value = "modelPrediction",
                      tags$br(),
                      sidebarLayout(
                        sidebarPanel2(fluid = FALSE,
                          tags$br(),
                          tags$h5(tags$b("Upload a model or select an existing one.")),
                          tags$br(),
                          fileInput(inputId = "usermodel", label = "Import model"),
                          radioButtons(inputId = "model", label = "Already available", choices = c("Random Forest", "XGBoost"), selected = character(0)),
                          actionButton(inputId = "done", label = "Submit", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                          out = conditionalPanel(condition = "input.done > 0",
                                                 verbatimTextOutput(outputId = "pd"),
                                                 verbatimTextOutput(outputId = "modelDetails")
                                                 )
                        ),
                        mainPanel(
                          DT::dataTableOutput(outputId = "outTrueDE"),
                          conditionalPanel(condition = "input.done > 0", tags$hr()),
                          DT::dataTableOutput(outputId = "outFalseDE")
                        )
                      ),
                      # Footer
                      conditionalPanel(condition = "input.done > 0",
                                       div(downloadButton(outputId = "downloadTrue", label = "Download DE genes", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                           style="float:right"
                                       ),
                                       div(downloadButton(outputId = "downloadFalse", label = "Download non DE genes", width = "100px", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                           style="float:right"
                                       )
                                       )
                      )
                    )
      }) 
    }
  })
  
  # Data for Prediction Upload
  observeEvent(eventExpr = input$check, handlerExpr = {
    output$userData <- DT::renderDataTable({
      if(is.null(input$newdata)){
        return(NULL)
      }
      read.table(file = isolate(input$newdata$datapath), sep = isolate(input$sep), header = TRUE)
    })
    output$userDataSummary <- renderPrint({
      summary(read.table(file = isolate(input$newdata$datapath), sep = isolate(input$sep), header = TRUE))
    })
  })
  # Next Tab button in Prediction UI
  observeEvent(eventExpr = input$nextpan, handlerExpr = {
    if(input$check > 0){
      updateTabsetPanel(session, inputId = "inTab2", selected = "modelPrediction")    
    }else{
      showNotification(tags$h4("Please select a file and try again."), 
                       type = "warning",
                       duration = 5,
      )
    }
  })
  # Make Prediction Using the model
  observeEvent(eventExpr = input$done, handlerExpr = {
      indata <- read.table(file = isolate(input$newdata$datapath), sep = isolate(input$sep), header = TRUE)
      if(is.null(isolate(input$usermodel))){
        if(isolate(input$model) == "Random Forest"){
          # import model
          model <- readRDS(file = "./data/RFmodel.rds")
          output$pd <- renderText({"Predictor Variables."})
          output$modelDetails <- renderPrint({
            return(row.names(model[['importance']]))
          })
          # call function for random forest prediction
          predData <- rfPrediction(indata, model)
        }else if(isolate(input$model) == "XGBoost"){
          # import model
          model <- readRDS(file = "./data/XGBmodel.rds")
          output$pd <- renderText({"Predictor Variables."})
          output$modelDetails <- renderPrint({
            return(model[["feature_names"]])
          })
          predData <- xgbPrediction(indata, model)
        }
      }else{
        model <- readRDS(file = input$usermodel$datapath)
        if(length(model) == 9){
          output$pd <- renderText({"Predictor Variables."})
          output$modelDetails <- renderPrint({
            return(model[["feature_names"]])
          })
          # XGBoost case
          predData <- xgbPrediction(indata, model)
        }else{
          output$pd <- renderText({"Predictor Variables."})
          output$modelDetails <- renderPrint({
            return(row.names(model[['importance']]))
          })
          predData <- rfPrediction(indata, model)
        }
      }
      predData <- predData[,c(ncol(predData), 1:ncol(predData)-1)]
      output$outTrueDE <- DT::renderDataTable({
        return(predData[which(predData$DE == TRUE),])
      })
      output$outFalseDE <- DT::renderDataTable({
        return(predData[which(predData$DE == FALSE),])
      })  
      write.csv(x = predData[which(predData$DE == TRUE),], file = "./data/true_pred.csv", row.names = FALSE)
      write.csv(x = predData[which(predData$DE == FALSE),], file = "./data/false_pred.csv", row.names = FALSE)
  })
  
  output$downloadTrue <- downloadHandler(
    filename = function(){
      "TrueDE.csv"
    },
    content = function(file){
      write.csv(x = read.csv(file = "./data/true_pred.csv", header = TRUE), file = file)
    }
  )
  
  output$downloadFalse <- downloadHandler(
    filename = function(){
      "FalseDE.csv"
    },
    content = function(file){
      write.csv(x = read.csv(file = "./data/false_pred.csv", header = TRUE), file = file)
    }
  )
  
  # Get selected results reactively in UI for training models
  files_selection <- eventReactive(eventExpr = input$selected, {
    datalist <- list()
    if(input$file1 == TRUE){
      datalist[['kallisto']] <- sleuth_data
    }else{
      datalist[['kallisto']] <- data.frame(matrix(nrow = 0, ncol = ncol(sleuth_data)))
      colnames(datalist[['kallisto']]) <- colnames(sleuth_data)
    }
    if(input$file2 == TRUE){
      datalist[['tophat2']] <- tophat_data
    }else{
      datalist[['tophat2']] <- data.frame(matrix(nrow = 0, ncol = ncol(tophat_data)))
      colnames(datalist[['tophat2']]) <- colnames(tophat_data)
    }
    if(input$file3 == TRUE){
      datalist[['hisat2']] <- hisat2_data
    }else{
      datalist[['hisat2']] <- data.frame(matrix(nrow = 0, ncol = ncol(hisat2_data)))
      colnames(datalist[['hisat2']]) <- colnames(hisat2_data)
    }
    if(input$file4 == TRUE){
      datalist[['bitseq']] <- bitseq_data
    }else{
      datalist[['bitseq']] <- data.frame(matrix(nrow = 0, ncol = ncol(bitseq_data)))
      colnames(datalist[['bitseq']]) <- colnames(bitseq_data)
    }
    if(input$file5 == TRUE){
      datalist[['rsem']] <- rsem_data
    }else{
      datalist[['rsem']] <- data.frame(matrix(nrow = 0, ncol = ncol(rsem_data)))
      colnames(datalist[['rsem']]) <- colnames(rsem_data)
    }
    if(input$file6 == TRUE){
      datalist[['ebseq']] <- ebseq_data
    }else{
      datalist[['ebseq']] <- data.frame(matrix(nrow = 0, ncol = ncol(ebseq_data)))
      colnames(datalist[['ebseq']]) <- colnames(ebseq_data)
    }
    return(datalist)
  })
  # Print output tables
  output$kallisto <- DT::renderDataTable({
    datalist <- files_selection()
    datalist[['kallisto']]
  })
  output$tophat2 <- DT::renderDataTable({
    datalist <- files_selection()
    datalist[['tophat2']]
  })
  output$hisat2 <- DT::renderDataTable({
    datalist <- files_selection()
    datalist[['hisat2']]
  })
  output$bitseq <- DT::renderDataTable({
    datalist <- files_selection()
    datalist[['bitseq']]
  })
  output$rsem <- DT::renderDataTable({
    datalist <- files_selection()
    datalist[['rsem']]
  })
  output$ebseq <- DT::renderDataTable({
    datalist <- files_selection()
    datalist[['ebseq']]
  })
  
  # Print Summary
  output$kallistoSummary <- renderPrint({
    datalist <- files_selection()
    if(!is.null(datalist[['kallisto']])){
      summary(datalist[['kallisto']])    
    }
  })
  output$tophat2Summary <- renderPrint({
    datalist <- files_selection()
    if(!is.null(datalist[['tophat2']])){
      summary(datalist[['tophat2']])    
    }
  })
  output$hisat2Summary <- renderPrint({
    datalist <- files_selection()
    if(!is.null(datalist[['hisat2']])){
      summary(datalist[['hisat2']])    
    }
  })
  output$bitseqSummary <- renderPrint({
    datalist <- files_selection()
    if(!is.null(datalist[['bitseq']])){
      summary(datalist[['bitseq']])    
    }
  })
  output$rsemSummary <- renderPrint({
    datalist <- files_selection()
    if(!is.null(datalist[['rsem']])){
      summary(datalist[['rsem']])    
    }
  })
  output$ebseqSummary <- renderPrint({
    datalist <- files_selection()
    if(!is.null(datalist[['ebseq']])){
      summary(datalist[['ebseq']])    
    }
  })
  
  # Update Tab Panel when Next button is activated in selectFilePanel
  observeEvent(eventExpr = input$next1, handlerExpr = {
    if(input$selected > 0){
      updateTabsetPanel(session, inputId = "inTab", selected = "selectVariablePanel")    
    }else{
      showNotification(tags$h4("No Pipelines submitted."), 
                       type = "error",
                       duration = 5,
                       )
    }
    
  })
  
  # Update Tab Panel when Back button is avtivated in selectVariablePanel
  observeEvent(eventExpr = input$back1, handlerExpr = {
    updateTabsetPanel(session, inputId = "inTab", selected = "selectFilePanel")
  })
  
  # Update Tab Panel when Next button is activated in selectVariablePanel
  observeEvent(eventExpr = input$next2, handlerExpr = {
    l = length(input$kallistoCols) + length(input$tophat2Cols) + length(input$hisat2Cols) + 
      length(input$bitseqCols) + length(input$rsemCols) + length(input$ebseqCols)
    if(l == 0){
      showNotification(tags$h4("Number of selected predictors is equal to zero."), type = "error")
    }else{
      updateTabsetPanel(session, inputId = "inTab", selected = "modelPredictionPanel")  
    }
  })
  
  # Model Selection
  observeEvent(eventExpr = input$submitModel, handlerExpr = {
      # select data given from user
      datalist <- files_selection()
      if(nrow(datalist[['kallisto']]) == 0){
        sleuth_data <- sleuth_data[0,]
      }
      if(nrow(datalist[['tophat2']]) == 0){
        tophat_data <- tophat_data[0,]
      }
      if(nrow(datalist[['hisat2']]) == 0){
        hisat2_data <- hisat2_data[0,]
      }
      if(nrow(datalist[['bitseq']]) == 0){
        bitseq_data <- bitseq_data[0,]
      }
      if(nrow(datalist[['rsem']]) == 0){
        rsem_data <- rsem_data[0,]
      }
      if(nrow(datalist[['ebseq']]) ==0){
        ebseq_data <- ebseq_data[0,]
      }
      # Select user-defined columns
      sleuth_data <- sleuth_data[, names(sleuth_data)[names(sleuth_data) %in% input$kallistoCols]]
      tophat_data <- tophat_data[, names(tophat_data)[names(tophat_data) %in% input$tophat2Cols]]
      hisat2_data <- hisat2_data[, names(hisat2_data)[names(hisat2_data) %in% input$hisat2Cols]]
      bitseq_data <- bitseq_data[, names(bitseq_data)[names(bitseq_data) %in% input$bitseqCols]]
      rsem_data <- rsem_data[, names(rsem_data)[names(rsem_data) %in% input$rsemCols]]
      ebseq_data <- ebseq_data[, names(ebseq_data)[names(ebseq_data) %in% input$ebseqCols]]
      # Merge Data
      finalData <- compineResults(sleuth_data, tophat_data, hisat2_data, bitseq_data, rsem_data, ebseq_data, simData)
      m <- isolate(input$modelSelection)
      nModels <- isolate(input$numModels)
      if(m == "Random Forest"){
        # RF
        mlModel <- randomForestModel(finalData, nModels)
      }else if(m == "XGBoost"){
        # XGBoost
        mlModel <- xGBoostModel(finalData, nModels)
      }
      saveRDS(object = mlModel[['MODEL']], file = "./data/model.rds")
      output$modelMetrics <- renderUI({
        str1 <- tags$b(paste("Model Accuracy = ", mlModel[['ACCURACY']]))
        str2 <- tags$b(paste("Area Under the Curve = ", mlModel[['AUC']]))
        str3 <- tags$b(paste("F1 Score of Classifier = ", mlModel[['F1SCORE']]))
        HTML(paste(str1, str2, str3, sep = '<br/>'))
      })
      output$model <- renderPrint({
        return(mlModel)
      }) 
  })
  
  output$downloadModel <- downloadHandler(
    filename = function(){
      paste0(input$modelSelection, ".rds")
    },
    content = function(file){
      saveRDS(object = readRDS(file = "./data/model.rds"), file = file)
    }
  )
  
  # Update Tab Panel when Back button is activated in modelPredictionPanel
  observeEvent(input$back2,{
    updateTabsetPanel(session, inputId = "inTab", selected = "selectVariablePanel")
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
