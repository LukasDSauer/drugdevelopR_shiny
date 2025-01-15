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
# Potentially more than one trial phase III #
###############################################
# d2: total sample size in phase II 
# d3: total sample size in phase III

# 1. Strategy: Only best promising treatment goes to phase III
# -> Phase III is always 2 arm trial (1:1 sample size allocatiob)
# 2. Strategy: All promising treatments go to phase III
# -> Phase III is 2 or 3 arm trial (1:1 or 1:1:1 sample size allocatiob)

# mainPath <- "/opt/shiny-server/samplesizr/multitrial/"
mainPath <- "./results/"

shinyServer(function(input, output,session) {
   
     output$table <- renderTable({
        
         input$go
         case = isolate(input$Case)
         input$go
         d2 = isolate(input$D2)
         input$go
         HRgo = isolate(input$HRgo)
         input$go
         stepd2 = isolate(input$stepD2)
         input$go
         stepHRgo = isolate(input$stepHRgo)

         input$go
         beta = isolate(input$beta)
         input$go
         p2 = isolate(input$p2)
         input$go
         p3 = isolate(input$p3)
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
         HR = isolate(input$HR)
         input$go
         refresh = isolate(input$refresh)
         
         
         if(case==1){
            # Strategy 1alpha vs. Strategy 1/2,
            STRATEGY = c(1, 2)
            }
         if(case==2){
            # Strategy 1alpha^2 vs. Strategy 2/2 vs. Strategy 2/3 vs. Strategy 2/2( + 1)
            STRATEGY = c(1, 2, 3)
            }
         if(case==3){
            # Strategy 1alpha^3 vs. Strategy 3/3 vs. Strategy 3/4
            STRATEGY = c(1, 3, 4)
            }
         
         
         for(strategy in STRATEGY){
           input$go
           alpha = isolate(input$alpha)
           y <- function(x){
             result <- optimal_multitrial(
               w = NULL,
               hr1 = HR,
               hr2 = NULL,
               id1 = NULL,
               id2 = NULL,
               d2min = d2[1],
               d2max = d2[2],
               stepd2 = stepd2,
               hrgomin = HRgo[1],
               hrgomax = HRgo[2],
               stephrgo = stepHRgo,
               alpha = alpha,
               beta = beta,
               xi2 = p2,
               xi3 = p3,
               c2 = c2,
               c3 = c3,
               c02 = c02,
               c03 = c03,
               K = Inf,
               N = Inf,
               S = -Inf,
               b1 = b1,
               b2 = b2,
               b3 = b3,
               case = case,
               strategy = strategy,
               fixed = TRUE,
               num_cl = 1
             )
             
             return(result)
           }
         
         result <- progressr::withProgressShiny({
           y()
         }, message = "Optimization progress", 
         detail = paste("for strategy", strategy))
      
         
         if(refresh==1 & strategy==1){
            DF = NULL
         } else {
           load(file=paste0(mainPath, "optimizationresults.RData"))
         }
         DF <- rbind(DF, result)
         save(result, file = paste0(mainPath, "optimizationresult_strat", 
                                    strategy,".RData"))
         save(DF, file = paste0(mainPath, "optimizationresults.RData"))
         }
         # Remove unnecessary columns
         DF_out <- DF
         DF_out <- DF[, !names(DF) %in% c("K", "N", "S", "gamma")]
         return(DF_out)
    })

     output$plot <- renderPlotly({
           
           input$go
           case = isolate(input$Case)
           input$go
           Plot = isolate(input$Plot)
           if(Plot==1){
             load(file = paste0(mainPath, "optimizationresults.RData"))
             xid <- "hrgo"
             yid <- "d2"
             xlab <- list(title = "HRgo")
             ylab <- list(title = "d2")
             
             zid <- "ufkt"
             zlab <- list(title = "expected utility")
             ufkt_vec <- DF[,"u"]
             bests <- DF[which.max(ufkt_vec),"Strategy"]
             load(file = paste0(mainPath, 
                                paste0("optimizationresult_strat", 
                                       bests,".RData")))
             trace <- attr(result, "trace")
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
