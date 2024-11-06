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

#mainPath <- "/opt/shiny-server/samplesizr/basic/"
mainPath <- "./results/"

shinyServer(function(input, output, session) {
  output$table <- renderTable({
    Select = input$Select
    if (Select == 1) {
      input$go1
      Num1 = isolate(input$Num1)
      input$go1
      HRgo = isolate(input$HRgo)
      input$go1
      stepNum1 = isolate(input$stepNum1)
      input$go1
      stepHRgo = isolate(input$stepHRgo)
      input$go1
      HR = isolate(input$HR)
      input$go1
      alpha = isolate(input$alpha1)
      input$go1
      beta = isolate(input$beta1)
      input$go1
      p2 = isolate(input$p2)
      input$go1
      p3 = isolate(input$p3)
      input$go1
      c02 = isolate(input$c021)
      input$go1
      c03 = isolate(input$c031)
      input$go1
      c2 = isolate(input$c21)
      input$go1
      c3 = isolate(input$c31)
      input$go1
      steps1 = isolate(input$small1)
      input$go1
      stepm1 = isolate(input$medium1)
      input$go1
      stepl1 = isolate(input$large1)
      input$go1
      b1 = isolate(input$b11)
      input$go1
      b2 = isolate(input$b21)
      input$go1
      b3 = isolate(input$b31)
      input$go1
      refresh = isolate(input$refresh1)
      input$go1
      K = isolate(input$K1)
    }
    if (Select == 2) {
      input$go2
      p0 = isolate(input$p0)
      input$go2
      p1 = isolate(input$p1)
      input$go2
      Num2 = isolate(input$Num2)
      input$go2
      RRgo = isolate(input$RRgo)
      input$go2
      stepNum2 = isolate(input$stepNum2)
      input$go2
      stepRRgo = isolate(input$stepRRgo)
      input$go2
      alpha = isolate(input$alpha2)
      input$go2
      beta = isolate(input$beta2)
      input$go2
      p2 = isolate(input$p22)
      input$go2
      p3 = isolate(input$p32)
      input$go2
      c02 = isolate(input$c022)
      input$go2
      c03 = isolate(input$c032)
      input$go2
      c2 = isolate(input$c22)
      input$go2
      c3 = isolate(input$c32)
      input$go2
      steps1 = isolate(input$small2)
      input$go2
      stepm1 = isolate(input$medium2)
      input$go2
      stepl1 = isolate(input$large2)
      input$go2
      b1 = isolate(input$b12)
      input$go2
      b2 = isolate(input$b22)
      input$go2
      b3 = isolate(input$b32)
      input$go2
      refresh = isolate(input$refresh2)
      input$go2
      K = isolate(input$K2)
    }
    if (Select == 3) {
      input$go3
      Delta = isolate(input$Delta)
      input$go3
      Num3 = isolate(input$Num3)
      input$go3
      kappa = isolate(input$kappa)
      input$go3
      stepNum3 = isolate(input$stepNum3)
      input$go3
      stepkappa = isolate(input$stepkappa)
      input$go3
      alpha = isolate(input$alpha3)
      input$go3
      beta = isolate(input$beta3)
      input$go3
      p2 = isolate(input$p23)
      input$go3
      p3 = isolate(input$p33)
      input$go3
      c02 = isolate(input$c023)
      input$go3
      c03 = isolate(input$c033)
      input$go3
      c2 = isolate(input$c23)
      input$go3
      c3 = isolate(input$c33)
      input$go3
      steps1 = isolate(input$small3)
      input$go3
      stepm1 = isolate(input$medium3)
      input$go3
      stepl1 = isolate(input$large3)
      input$go3
      b1 = isolate(input$b13)
      input$go3
      b2 = isolate(input$b23)
      input$go3
      b3 = isolate(input$b33)
      input$go3
      refresh = isolate(input$refresh3)
      input$go3
      K = isolate(input$K3)
    }
    if (Select == 1) {
      y = function() {
        result <- optimal_tte(
          w = NULL,
          hr1 = HR,
          hr2 = NULL,
          id1 = NULL,
          id2 = NULL,
          d2min = Num1[1],
          d2max = Num1[2],
          stepd2 = stepNum1,
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
          K = K,
          N = Inf,
          S = -Inf,
          steps1 = steps1,
          stepm1 = stepm1,
          stepl1 = stepl1,
          b1 = b1,
          b2 = b2,
          b3 = b3,
          gamma = 0,
          fixed = TRUE,
          skipII = FALSE,
          num_cl = 1
        )
        return(result)
      }
    }
    if (Select == 2) {
      y = function() {
        result <- optimal_binary(
          w = NULL,
          p0 = p0,
          p11 = p1,
          p12 = NULL,
          in1 = NULL,
          in2 = NULL,
          n2min = Num2[1],
          n2max = Num2[2],
          stepn2 = stepNum2,
          rrgomin = RRgo[1],
          rrgomax =  RRgo[2],
          steprrgo = stepRRgo,
          alpha = alpha,
          beta = beta,
          c2 = c2,
          c3 = c3,
          c02 = c02,
          c03 = c03,
          K = K,
          N = Inf,
          S = -Inf,
          steps1 = steps1,
          stepm1 = stepm1,
          stepl1 = stepl1,
          b1 = b1,
          b2 = b2,
          b3 = b3,
          gamma = 0,
          fixed = TRUE,
          skipII = FALSE,
          num_cl = 1
        )
        return(result)
      }
    }
    if (Select == 3) {
      y = function() {
        result <- optimal_normal(
          w = NULL,
          Delta1 = Delta,
          Delta2 = NULL,
          in1 = NULL,
          in2 = NULL,
          a = NULL,
          b = NULL,
          n2min = Num3[1],
          n2max = Num3[2],
          stepn2 = stepNum3,
          kappamin = kappa[1],
          kappamax = kappa[2],
          stepkappa = stepkappa,
          alpha = alpha,
          beta = beta,
          c2 = c2,
          c3 = c3,
          c02 = c02,
          c03 = c03,
          K = K,
          N = Inf,
          S = -Inf,
          steps1 = steps1,
          stepm1 = stepm1,
          stepl1 = stepl1,
          b1 = b1,
          b2 = b2,
          b3 = b3,
          gamma = 0,
          fixed = TRUE,
          skipII = FALSE,
          num_cl = 1
        )
        return(result)
      }
    }
    
    result <- progressr::withProgressShiny({
      y()
    }, message = "Optimization progress")
    
    if (refresh == 0) {
      load(file = paste0(mainPath, "optimizationresults.RData"))
    } else {
      DF = NULL
    }
    DF <- rbind(DF, result)
    save(result, file = paste0(mainPath, "optimizationresult_last.RData"))
    save(DF, file = paste0(mainPath, "optimizationresults.RData"))
    # Remove unnecessary columns
    DF_out <- DF[, !names(DF) %in% c("skipII", "N", "S", "gamma")]
    return(DF_out)
  })
  
  output$plot <- renderPlotly({
    Select = input$Select
    if (Select == 1) {
      input$go1
      Plot = input$Plot1
    }
    if (Select == 2) {
      input$go2
      Plot = input$Plot2
    }
    if (Select == 3) {
      input$go3
      Plot = input$Plot3
    }
    if (Plot == 1) {
      load(file = paste0(mainPath, "optimizationresult_last.RData"))
      
      if (Select == 1) {
        xid <- "hrgo"
        yid <- "d2"
        xlab <- list(title = "HRgo")
        ylab <- list(title = "d2")
      } else {
        if (Select == 2) {
          xid <- "rrgo"
          yid <- "n2"
          xlab <- list(title = "RRgo")
          ylab <- list(title = "n2")
        } else {
          xid <- "kappa"
          yid <- "n2"
          xlab <- list(title = "kappa")
          ylab <- list(title = "n2")
        }
      }
      trace <- attr(result, "trace")
      zid <- "ufkt"
      zmat <- t(trace[c(xid, yid, zid), ]) %>% 
        as.data.frame() %>% 
        pivot_wider(names_from = all_of(yid), values_from = all_of(zid))
      x <- zmat[[xid]]
      y <- as.numeric(colnames(zmat)[-1])
      xmat <- matrix(x, nrow = length(y), ncol = length(x), byrow=TRUE)
      ymat <- matrix(y, nrow = length(y), ncol = length(x))
      zmat <- t(as.matrix(select(zmat, -any_of(xid))))
      rownames(zmat) <- NULL
      save(zmat, file = "results/test.RData")
      zlab <- list(title = "expected utility")
      collab <- zlab
      plot_ly(
        x = xmat,
        y = ymat,
        z = zmat,
        type = "surface"
      ) %>%
        layout(
          title =
            paste0("Optimization region"),
          scene = list(
            xaxis = xlab,
            yaxis = ylab,
            zaxis = zlab
          )
        ) %>% 
        return(.)
    }
  })
})
