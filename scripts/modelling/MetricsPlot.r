library(plotly)
# read files
rf <- readRDS(file = "./test/rfAccAucF1.rds")
xgb <- readRDS(file = "./test/xgbAccAucF1.rds")
# random forest accuracy
figrfAcc <- plot_ly(data = rf, x = ~Accuracy, y = ~trainingLength, z = ~testingLength, 
                  type = "scatter3d", mode = "lines+markers", opacity = 1,
                  line = list(width = 6, color = "purple4", reverscale = F)) %>%
  layout(                       
    title = "Accuracy",
    scene = list(
      xaxis = list(           
        title = "Accuracy",     
        showgrid = T),       
      yaxis = list(           
        title = "Training Data length",
        showgrid = T),
      zaxis = list(
        title = "Testing Data length",
        showgrid = T)
    )
  )
# orca(p = figAcc, file = "./test/rfAccFig.png")
# random forest auc
figrfAuc <- plot_ly(data = rf, x = ~AUC, y = ~trainingLength, z = ~testingLength, 
                  type = "scatter3d", mode = "lines+markers", opacity = 1,
                  line = list(width = 6, color = "purple4", reverscale = F)) %>%
  layout(                       
    title = "Area Under the Curve",
    scene = list(
      xaxis = list(           
        title = "AUC",     
        showgrid = T),       
      yaxis = list(           
        title = "Training Data length",
        showgrid = T),
      zaxis = list(
        title = "Testing Data length",
        showgrid = T)
    )
  )
# random forest f1 score
figrfF1 <- plot_ly(data = rf, x = ~F1score, y = ~trainingLength, z = ~testingLength, 
                  type = "scatter3d", mode = "lines+markers", opacity = 1,
                  line = list(width = 6, color = "purple4", reverscale = F)) %>%
  layout(                       
    title = "F1 Measure",
    scene = list(
      xaxis = list(           
        title = "F1 score",     
        showgrid = T),       
      yaxis = list(           
        title = "Training Data length",
        showgrid = T),
      zaxis = list(
        title = "Testing Data length",
        showgrid = T)
    )
  )

# xgboost accuracy
figxgbAcc <- plot_ly(data = xgb, x = ~Accuracy, y = ~trainingLength, z = ~testingLength, 
                    type = "scatter3d", mode = "lines+markers", opacity = 1,
                    line = list(width = 6, color = "purple4", reverscale = F)) %>%
  layout(                       
    title = "Accuracy",
    scene = list(
      xaxis = list(           
        title = "Accuracy",     
        showgrid = T),       
      yaxis = list(           
        title = "Training Data length",
        showgrid = T),
      zaxis = list(
        title = "Testing Data length",
        showgrid = T)
    )
  )
# orca(p = figAcc, file = "./test/rfAccFig.png")
# xgb auc
figxgbAuc <- plot_ly(data = xgb, x = ~AUC, y = ~trainingLength, z = ~testingLength, 
                    type = "scatter3d", mode = "lines+markers", opacity = 1,
                    line = list(width = 6, color = "purple4", reverscale = F)) %>%
  layout(                       
    title = "Area Under the Curve",
    scene = list(
      xaxis = list(           
        title = "AUC",     
        showgrid = T),       
      yaxis = list(           
        title = "Training Data length",
        showgrid = T),
      zaxis = list(
        title = "Testing Data length",
        showgrid = T)
    )
  )
# xgb f1 score
figxgbF1 <- plot_ly(data = xgb, x = ~F1score, y = ~trainingLength, z = ~testingLength, 
                   type = "scatter3d", mode = "lines+markers", opacity = 1,
                   line = list(width = 6, color = "purple4", reverscale = F)) %>%
  layout(                       
    title = "F1 Measure",
    scene = list(
      xaxis = list(           
        title = "F1 score",     
        showgrid = T),       
      yaxis = list(           
        title = "Training Data length",
        showgrid = T),
      zaxis = list(
        title = "Testing Data length",
        showgrid = T)
    )
  )
