basicSegmentCorrection <- function(series, changepoints = c(), bestFitsBySegment, segmentsToCorrect = c(), referenceSegment=NA, verbose=T, plotResults=T){

  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }
  if(length(changepoints) < 1){
    stop("The changepoints must be in a vector and there must be at least 1 changepoint.")
  }
  if(!is.na(referenceSegment)){
    if(referenceSegment > length(bestFitsBySegment)){
      stop("The index for the reference segment is larger than the number of segments in the series.")
    }
    if(referenceSegment < 0){
      stop("The index for the reference segment is negative.")
    }
  }


  #If changepoints is NA, find them based on the best fits by segment.
  if(length(changepoints) == 0){
    changepoints <- c()

    for(i in 1:length(bestFitsBySegment)){
      changepoints <- c(changepoints, bestFitsBySegment[[i]]$end)
    }
  }
  #Include 0 and n as endpoints.
  changepoints <- sort(unique(c(0, changepoints, length(series))))

  #Iterate over each segment, recording the corrected series
  correctedSeries <- series

  #If the segments to correct is NA, correct all by default.
  if(length(segmentsToCorrect) == 0){
    segmentsToCorrect <- 1:(length(changepoints)-1)
  }

  #Iterate over segments, performing de-trending and de-seasonalizing where desired.
  for(i in segmentsToCorrect){
    #De-trend or de-seasonalize based on the best fit
    subSeries <- bestFitsBySegment[[i]]$model$residuals
    correctedSeries[(changepoints[i]+1):changepoints[i+1]] <- subSeries
  }

  #Get the mean and standard deviation of the reference segment
  if(is.na(referenceSegment)){
    referenceSegment <- which.max(diff(changepoints))
    subSeries <- series[bestFitsBySegment[[referenceSegment]]$start:bestFitsBySegment[[referenceSegment]]$end]
    if(length(subSeries) == 1){
      stop("The reference segment has length 1. Please choose a different reference segment.")
    }
    refMean <- mean(subSeries)
    refSD <- summary(bestFitsBySegment[[referenceSegment]]$model)$sigma
  }
  else if(referenceSegment == 0){
    refMean <- 0
    refSD <- 1
  }
  else{
    subSeries <- series[(bestFitsBySegment[[referenceSegment]]$start):(bestFitsBySegment[[referenceSegment]]$end)]
    if(length(subSeries) == 1){
      stop("The reference segment has length 1. Please choose a different reference segment.")
    }
    refMean <- mean(subSeries)
    refSD <- summary(bestFitsBySegment[[referenceSegment]]$model)$sigma
  }

  if(refSD == 0){
    stop("The reference segment has a standard deviation of 0. Please choose a different reference segment.")
  }

  if(is.na(refSD)){
    stop("The reference segment has an NA residual standard error. Please choose a different reference segment.")
  }

  #Iterate over segments, shifting and scaling to match the reference segment.
  for(i in segmentsToCorrect){
    subSeries <- correctedSeries[(changepoints[i]+1):changepoints[i+1]]
    subSeries <- (subSeries - mean(subSeries))/sd(subSeries)
    subSeries <- refSD*subSeries + refMean
    correctedSeries[(changepoints[i]+1):changepoints[i+1]] <- subSeries
  }

  #Plot the corrected series, if desired.
  if(plotResults){
    plot(y=correctedSeries, x=1:length(series), main="Corrected Series", xlab="Index", ylab="Value", type="l")
    for(i in 1:(length(changepoints)-1)){
      abline(v=changepoints[i+1], col="red")
    }
  }

  return(correctedSeries)

}
