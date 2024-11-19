library("shiny")
library("ggplot2")
library("gridExtra")
library("mvtnorm")
library("doParallel")
library("parallel")
library("cubature")
library("plotly")
library("tidyr")
library("drugdevelopR")

###############################################
# Potentially more than two arms in phase III #
###############################################
# n2: total sample size in phase II 
# n3: total sample size in phase III

# 1. Strategy: Only best promising treatment goes to phase III
# -> Phase III is always 2 arm trial (1:1 sample size allocation)
# 2. Strategy: All promising treatments go to phase III
# -> Phase III is 2 or 3 arm trial (1:1 or 1:1:1 sample size allocation)

mainPath <- "/opt/shiny-server/samplesizr/multiarm/"
# mainPath <- "./results/"

shinyServer(function(input, output,session) {
   
     output$table <- renderTable({
        
         input$go
         Select = isolate(input$Select)
         input$go
         N2 = isolate(input$N2)
         input$go
         HRgo = isolate(input$HRgo)
         input$go
         stepN2 = isolate(input$stepN2)
         input$go
         stepHRgo = isolate(input$stepHRgo)
         input$go
         alpha = isolate(input$alpha)
         input$go
         beta = isolate(input$beta)
         input$go
         ec = isolate(input$ec)
         input$go
         R = isolate(input$R)
         input$go
         c02 = isolate(input$c02)
         input$go
         c03 = isolate(input$c03)
         input$go
         c2 = isolate(input$c2)
         input$go
         c3 = isolate(input$c3)
         input$go
         b1 = isolate(input$b1)
         input$go
         b2 = isolate(input$b2)
         input$go
         b3 = isolate(input$b3)
         input$go
         HR1 = isolate(input$HR1)
         input$go
         HR2 = isolate(input$HR2)
         input$go
         refresh = isolate(input$refresh)
         input$go
         Nc = isolate(input$Nc)
         
         if(Select==1){
           strategy = 1
           
           }
         if(Select==2){
           strategy = 2
           }
         if(Select==3){
           strategy = 3
           }
         
         y <- function(){
           result <- optimal_multiarm(
             hr1 = HR1,
             hr2 = HR2,
             ec = ec,
             n2min = N2[1],
             n2max = N2[2],
             stepn2 = stepN2,
             hrgomin = HRgo[1],
             hrgomax = HRgo[2],
             stephrgo = stepHRgo,
             alpha = alpha,
             beta = beta,
             c2 = c2,
             c3 = c3,
             c02 = c02,
             c03 = c03,
             K = Inf,
             N = Nc,
             S = -Inf,
             steps1 = 1,
             stepm1 = 0.95,
             stepl1 = 0.85,
             b1 = b1,
             b2 = b2,
             b3 = b3,
             strategy = strategy,
             num_cl = 1
           )
           return(result)
         }
         
         

         result <- progressr::withProgressShiny({ y() },
          message = "Optimization progress", 
          detail = paste("for strategy", strategy))

         if(refresh==0){
           load(file=paste0(mainPath, "optimizationresults.RData"))
         } else {
             DF = NULL
         }
         
         DF <- rbind(DF, result)
         save(result, file = paste0(mainPath, "optimizationresult_last.RData"))
         save(DF, file = paste0(mainPath, "optimizationresults.RData"))
         DF_out <- DF
         DF_out <- DF[, !names(DF) %in% c("K", "S")]
         return(DF_out)
         })

     output$plot <- renderPlotly({
           input$go
           Select = isolate(input$Select)
           input$go
           Plot = isolate(input$Plot)
           
           if(Plot==1){
             load(file = paste0(mainPath, "optimizationresult_last.RData"))
             xid <- "hrgo"
             yid <- "n2"
             xlab <- list(title = "HRgo")
             ylab <- list(title = "n2")
             trace <- attr(result, "trace")
             zid <- "ufkt"
             zlab <- list(title = "expected utility")
             
              if(Select==1){
                bests = 1
              }
              if(Select==2){
                bests = 2
              }
              if(Select==3){
                # Show only the plot with the larger utility
                if (result[which(result$Strategy == 1),"u"] > 
                    result[which(result$Strategy == 2),"u"]) {
                  trace <- trace[, which(trace["strat",]==1)]
                  bests <- 1
                } else {
                  trace <- trace[,which(trace["strat",]==2)]
                  bests <- 2
                }
                 
              }
             zmat <- t(trace[c(xid, yid, zid), ]) %>% 
               as.data.frame() %>% 
               pivot_wider(names_from = all_of(yid), values_from = all_of(zid))
             x <- zmat[[xid]]
             y <- as.numeric(colnames(zmat)[-1])
             xmat <- matrix(x, nrow = length(y), ncol = length(x), byrow=TRUE)
             ymat <- matrix(y, nrow = length(y), ncol = length(x))
             zmat <- t(as.matrix(select(zmat, -any_of(xid))))
             rownames(zmat) <- NULL
              plot_ly(
                x = xmat,
                y = ymat,
                z = zmat,
                type = "surface"
              )  %>%
                layout(
                  title =
                    paste0("Optimization region of strategy ", bests),
                  scene = list(
                    xaxis = xlab,
                    yaxis = ylab,
                    zaxis = zlab
                  )
                )
              
           }

           
        })   
        
     
})
