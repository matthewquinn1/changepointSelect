getBestFitBySegment <- function(series, changepoints, thresholdLinear=1.5, thresholdSeasonal=1.5, numHarmonics=2, verbose=T, plotFits=T){
  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }
  if(length(changepoints) < 1){
    stop("The changepoints must be in a vector and there must be at least 1 changepoint.")
  }

  #Include 0 and n as endpoints.
  changepoints <- unique(c(0, changepoints, length(series)))

  #Iterate over each segment, checking which model is best.
  results <- list()

  for(i in 1:(length(changepoints)-1)){

    message("Working on segment ", i, ".")

    segment <- series[(changepoints[i]+1):changepoints[i+1]]

    #Constant fit (i.e. mean of segment)
    constantFit <- lm(segment ~ 1)
    constantRMSE <- sqrt(mean(resid(constantFit)^2))
    segmentResult <- list(start=changepoints[i]+1, end=changepoints[i+1], type="constant", model=constantFit)

    #Get the linear regression
    if(length(segment) >= 3){
      linearFit <- getLinearFit(segment)
      linearRMSE <- sqrt(mean(resid(linearFit)^2))

      bestRMSE <- linearRMSE
      threshold <- thresholdLinear

      segmentResult <- list(start=changepoints[i]+1, end=changepoints[i+1], type="linear", model=linearFit)
    }
    else if(verbose){
      message("A linear regression is not appropriate for segment ", i, " due to the small number of observations in this segment.")
    }

    #Get the harmonic regression
    if(length(segment) >= (2+numHarmonics*2)){
      harmonicFit <- getHarmonicFit(segment)
      harmonicRMSE <- sqrt(mean(resid(harmonicFit)^2))

      if(linearRMSE > harmonicRMSE){
        bestRMSE <- harmonicRMSE
        threshold <- thresholdSeasonal
        segmentResult <- list(start=changepoints[i]+1, end=changepoints[i+1], type="harmonic", model=harmonicFit)
      }
    }
    else if(verbose){
      message("A harmonic regression is not appropriate for segment ", i, " due to the small number of observations in this segment.")
    }


    #If the better fitting model is not an improvement by a factor of linearThreshold or seasonalThreshold, record constant function as best.
    if(segmentResult[["type"]] != "constant"){
      #summaryInfo <- summary(segmentResult[["model"]])
      #pValue <- pf(summaryInfo$fstatistic[1],summaryInfo$fstatistic[2],summaryInfo$fstatistic[3],lower.tail=FALSE)

      #if(((constantRMSE/bestRMSE) < threshold) | (pValue >= alpha)){
       if((constantRMSE/bestRMSE) < threshold){
        segmentResult <- list(start=changepoints[i]+1, end=changepoints[i+1], type="constant", model=constantFit)
      }
    }

    results[[i]] <- segmentResult
  }

  if(plotFits){
    plot(y=series, x=1:length(series), main="Series", xlab="Index", ylab="Value", pch=19)
    for(i in 2:(length(changepoints)-1)){
      abline(v=changepoints[i], col="red")
    }
    for(i in 1:(length(changepoints)-1)){
      points(y=results[[i]][["model"]]$fitted.values,  x=results[[i]][["start"]]:results[[i]][["end"]], col="dodgerblue", type="l", lwd=3)
    }
  }

  return(results)
}
