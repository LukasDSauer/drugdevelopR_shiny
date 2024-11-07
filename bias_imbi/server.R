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


# mainPath <- "/opt/shiny-server/samplesizr/bias/"

mainPath <- "./results/"

shinyServer(function(input, output, session) {
  output$table <- renderTable({
    input$go
    Select = isolate(input$Select)
    input$go
    D2 = isolate(input$D2)
    input$go
    R = isolate(input$R)
    input$go
    alphaL = isolate(input$alphaL)
    input$go
    HRgo = isolate(input$HRgo)
    input$go
    stepD2 = isolate(input$stepD2)
    input$go
    stepHRgo = isolate(input$stepHRgo)
    input$go
    stepR = isolate(input$stepR)
    input$go
    stepalphaL = isolate(input$stepalphaL)
    input$go
    alpha = isolate(input$alpha)
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
    
    
    if (Select == 1) {
      meth <- "multiplicative"
      meth_lab <- paste(meth, "method")
      y = function() {
        result <- optimal_bias(
          w = NULL,
          hr1 = HR,
          hr2 = NULL,
          id1 = NULL,
          id2 = NULL,
          d2min = D2[1],
          d2max = D2[2],
          stepd2 = stepD2,
          hrgomin = HRgo[1],
          hrgomax = HRgo[2],
          stephrgo = stepHRgo,
          adj = meth,
          lambdamin = R[1],
          lambdamax = R[2],
          steplambda = stepR,
          alphaCImin = NULL,
          alphaCImax = NULL,
          stepalphaCI = NULL,
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
          steps1 = 1,
          stepm1 = 0.95,
          stepl1 = 0.85,
          b1 = b1,
          b2 = b2,
          b3 = b3,
          fixed = TRUE,
          num_cl = 1
        )
        return(result)
      }
    }
    
    if (Select == 2) {
      meth <- "additive"
      meth_lab <- paste(meth, "method")
      y = function() {
        result <- optimal_bias(
          w = NULL,
          hr1 = HRgo,
          hr2 = NULL,
          id1 = NULL,
          id2 = NULL,
          d2min = D2[1],
          d2max = D2[2],
          stepd2 = stepD2,
          hrgomin = HRgo[1],
          hrgomax = HRgo[2],
          stephrgo = stepHRgo,
          adj = meth,
          lambdamin = NULL,
          lambdamax = NULL,
          steplambda = NULL,
          alphaCImin = alphaL[1],
          alphaCImax = alphaL[2],
          stepalphaCI = stepalphaL,
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
          steps1 = 1,
          stepm1 = 0.95,
          stepl1 = 0.85,
          b1 = b1,
          b2 = b2,
          b3 = b3,
          fixed = TRUE,
          num_cl = 1
        )
        return(result)
      }
    }
    
    if (Select == 3) {
      meth <- "both"
      meth_lab <- paste(meth, "methods")
      y = function() {
        result <- optimal_bias(
          w = NULL,
          hr1 = HRgo,
          hr2 = NULL,
          id1 = NULL,
          id2 = NULL,
          d2min = D2[1],
          d2max = D2[2],
          stepd2 = stepD2,
          hrgomin = HRgo[1],
          hrgomax = HRgo[2],
          stephrgo = stepHRgo,
          adj = meth,
          lambdamin = R[1],
          lambdamax = R[2],
          steplambda = stepR,
          alphaCImin = alphaL[1],
          alphaCImax = alphaL[2],
          stepalphaCI = stepalphaL,
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
          steps1 = 1,
          stepm1 = 0.95,
          stepl1 = 0.85,
          b1 = b1,
          b2 = b2,
          b3 = b3,
          fixed = TRUE,
          num_cl = 1
        )
        return(result)
      }
    }
    
    result <- progressr::withProgressShiny({
      y()
    }, message = "Optimization progress", detail = paste("for", meth_lab))
    
    
    if (refresh == 0) {
      load(file = paste0(mainPath, "optimizationresults.RData"))
    } else {
      DF = NULL
    }
    
    DF <- rbind(DF, result)
    save(result, file = paste0(mainPath, "optimizationresult_last.RData"))
    save(DF, file = paste0(mainPath, "optimizationresults.RData"))
    # Remove unnecessary columns
    DF_out <- DF
    # DF_out <- DF[, !names(DF) %in% c("skipII", "N", "S", "gamma")]
    return(DF_out)
  })
  
  output$plot <- renderPlotly({
    input$go
    Select = isolate(input$Select)
    input$go
    Plot = isolate(input$Plot)
    browser()
    if (Plot == 1) {
      load(file = paste0(mainPath, "optimizationresult_last.RData"))
      xid <- "hrgo"
      yid <- "d2"
      xlab <- list(title = "HRgo")
      ylab <- list(title = "d2")
      trace <- attr(result, "trace")
      zid <- "ufkt"
      zlab <- list(title = "expected utility")
      if (Select == 1) {
        showplot = "multiplicatively"
        adj <- result[,"Adj"]
      }
      if (Select == 2) {
        showplot = "additively"
        adj <- result[,"Adj"]
      }
      if (Select == 3) {
 
        # Show only the plot with the larger utility
        if (result[which(result$Method == "multipl."),"u"] > 
            result[which(result$Method == "add."),"u"]) {
          trace <- trace[,which(trace["strat",]=="multipl.")]
          showplot <- "multiplicatively"
          adj <- result[which(result$Method == "multipl."),"Adj"]
        } else{
          trace <- trace[,which(trace["strat",]=="add.")]
          showplot <- "additively"
          adj <- result[which(result$Method == "add."),"Adj"]
        }
        
      }
      # Only select the best adjustment value
      trace <- trace[,,which(trace["adj",]==adj)]
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
            paste0("Optimization region of ", showplot, " adjusted setting"),
          scene = list(
            xaxis = xlab,
            yaxis = ylab,
            zaxis = zlab
          )
        )
    }
  })
  
  
})
