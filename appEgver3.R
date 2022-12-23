
rm(list = ls())

# library needed

library(shiny)
library(plotly)
library(DT)   #
#library(shinydashboard)  # for the box
##library required for models 
library(tidyverse)
library(investr)



# Define UI ----
ui <- fluidPage(
  titlePanel("BCA Automation"),
  
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", label = "Data for model training", accept = ".txt", multiple = TRUE),
      
      fileInput("file2", label = "Data for prediction", accept = ".xlsx", multiple = TRUE),
      downloadButton('bca_report')
      
    ),
    
    mainPanel(
      
      fluidPage (
        
        fluidRow(
          column(7, h3("Estimates & PredictIntervals"),plotOutput("modelPlot")), 
          column(5,h3("Heat Map "), plotOutput("heatPlot"))
            
             )
       
         )
    
    )
 
    
     ), 
  
  tabsetPanel(
    tabPanel("Prediction Data", DT::dataTableOutput("tablePredict")),
    tabPanel("Reference Data", DT::dataTableOutput("tableRef"))
     
    )
     
  
  )

# Define server logic ----
server <- function(input, output) {
  observe({
    file1 = input$file1
    file2 = input$file2
    if (is.null(file1) || is.null(file2)) {
      return(NULL)
    }
    #1a. process the first data set
    #----------------------------
    raw_BCA = read.delim(file1$datapath,skip = 2, fileEncoding = "utf-16", 
                         sep="\t", as.is = T )
    raw_BCA1 <- raw_BCA%>% dplyr::select(X1:X12) %>% slice(1:8)
    
    #RESHAPING OF THE DATA FROM A PLATE VIEW TO A SIMPLE TABLE WITH POSITIONS
    raw_BCA <- raw_BCA1 %>% mutate(row=LETTERS[1:8]) %>%
      gather("column","abs",X1:X12) %>%
      mutate(column = parse_number(column)) %>%
      unite(position, c(row,column), sep="")
    
    #1b. process the second data set
    #--------------------------------------
    plan = xlsx::read.xlsx(file2$datapath, 
                           sheetIndex = 1) %>% dplyr::select(X1:X12)
    #RESHAPING OF THE PLATE TO MAKE A SIMPLE TABLE
    plan <- plan %>% mutate(row=LETTERS[1:8]) %>% 
      gather("column","sample",X1:X12) %>% 
      mutate(column = parse_number(column)) %>% 
      unite(position, c(row,column), sep="")
    #
    #MATCHING THE NAMEs from the USER INPUT and THE INTENSITY FROM THE RAW DATA
    BCA <- left_join(plan,raw_BCA, by="position") %>%
      separate(sample, c("sample","dilution"),sep="@") 
    
    #OUTPUT - plot1: heatplot
    #-------------------------------
    
    output$heatPlot <- renderPlot({
      pheatmap::pheatmap(as.matrix(raw_BCA1), 
                         cluster_rows = F,cluster_cols = F, 
                         color = viridis::viridis(96, option = "D"))
      
    })
    
    #2. build the model
    #-----------------------------------------------------------      
    # create the prediction interval range 
    BCA$lower.ci <- rep(NA, length(BCA$sample))
    BCA$upper.ci <- rep(NA, length(BCA$sample))
    BCA$est <- rep(NA, length(BCA$sample))
    
    BCA$dilution[is.na(BCA$dilution)]<-1
    BCA$dilution <- as.numeric(BCA$dilution)
    
    
    #indicates a reference sample
    BCA_model <- BCA %>% filter(grepl("STD",sample)) %>% 
      mutate(sample = parse_number(sample))
    
    BCA_model2 <- filter(BCA_model, sample != 0)
    
    # Range of values considered in regression
    min_abs = min(BCA_model2$abs)
    max_abs = max(BCA_model2$abs)
    
    # Modelling: 
    #-------------------------------------------------------
    #attach(BCA_model)
    
    fit3 <- lm(abs ~ sample + I(sample^2) + I(sample^3), data = BCA_model )
    summary(fit3)
    
    # Using the inverse function
    for (i in 1: length(BCA$abs)){
      if (min_abs <= BCA$abs[i] & BCA$abs[i] <= max_abs){
        
        possibleError <- tryCatch(invest(fit3, y0= BCA$abs[i], interval = "Wald"),
                                  error=function(e) e)
        
        if(inherits(possibleError, "error")) next
        inveVal <- invest(fit3, y0= BCA$abs[i], interval = "Wald")
        BCA$est[i] = inveVal[1]
        BCA$upper.ci[i] = inveVal[3]
        BCA$lower.ci[i] = inveVal[2]
      }
    }
    
    # provide the plot
    y0 <- seq(min_abs, max_abs, by = 0.01)
    y0lower.ci <- rep(NA, length(y0))
    y0upper.ci <- rep(NA, length(y0))
    y0est <- rep(NA, length(y0))
    for (i in 1: length(y0)){
      possibleError <- tryCatch(invest(fit3, y0= y0[i], interval = "Wald"),
                                error=function(e) e)
      
      if(inherits(possibleError, "error")) next
      res <- invest(fit3, y0 = y0[i], interval = "Wald")
      y0est[i] <- res[1]
      y0lower.ci[i] <- res[2]
      y0upper.ci[i] <- res[3]
    }
    
    
    
    # create means and compare
    BCA$weighEst<- as.numeric(BCA$est) * BCA$dilution
    BCA$range <- rapply(BCA$est, function(x){
      ifelse(is.na(x),"out of range", "valid") })
    
    
    Results_summary <- BCA %>%
      group_by(sample) %>%
      summarise(mean=mean(weighEst),
                sd=sd(weighEst),
                RSD = sd/mean*100,
                RSD_abs = sd(abs)/mean(abs)*100) %>% 
      mutate(ODvalid = case_when(RSD_abs<15 ~"valid",
                                 RSD_abs>=15 ~"invalid")) %>% ungroup()
    BCA_result <- merge(BCA[!is.na(BCA$sample),],Results_summary, 
                        by.x= "sample", by.y = "sample")
    
    BCA_ref <- BCA_result %>% filter(grepl("STD",sample))
    BCA_pred <- BCA_result %>% filter(grepl("Epool",sample))
    
    
    #OUTPUT - plot2: plot(BCA_model$sample, BCA_model$abs)
    #---------------------------------------------------------
    output$modelPlot <- renderPlot({ plot(BCA_model$sample, BCA_model$abs,pch=2,col = "blue", 
                                          xlab = "Concentration", ylab = "Absorbance" )
      points(BCA_pred$est,BCA_pred$abs, pch=16,col = "green")
      lines(y0upper.ci,y0, col = "green")
      lines(y0lower.ci,y0, col = "green")
      #abline(fit1, col = "red")
      regline <- function(x) return(fit3$coefficients[1] + fit3$coefficients[2] * x + fit3$coefficients[3] * x^2+
                                      fit3$coefficients[4] * x^3)
      x <- seq(0, 2000, by = 0.1)
      lines(x, regline(x), col = "blue")
      
    })
    
    # create a list for Rmarkdown 
    #listBCA = list(BCA_model,BCA_pred,y0upper.ci,y0lower.ci,y0)
    
    #############################################################
    #  table output 
    output$tablePredict <- DT::renderDataTable({BCA_pred})
    output$tableRef <- DT::renderDataTable({BCA_ref}) 
    
    
    #####Fragment Report Download ----
    output$bca_report <- downloadHandler(
      filename = function() {
        paste('bca-report', sep = '.', "docx")
      },
      
      content = function(file) {
        fragmentReport <- file.path(tempdir(), "report.Rmd") 
        file.copy("report.Rmd", fragmentReport, overwrite = TRUE)
        
        #Définir les variables qui seront dans le rapport (with params!!! is a list ???? yes list (x = 1, year = 2021,))
        params = list(
          rmd_heatmap = raw_BCA1,
          rmd_prediction = list(BCA_model, BCA_pred, y0upper.ci, y0lower.ci,y0),
          rmd_summary = fit3)
        
        rmarkdown::render(fragmentReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    ) # end of downloaderhandler
    
    
    
  })
  
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
