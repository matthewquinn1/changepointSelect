#Find the log-likelihood for a vector of observations given a vector of changepoints.
getLogLik <- function(series, changepoints, means=NA, sds=NA, segLengths=NA){

  if(any(is.na(series)) | any(is.na(changepoints))){
    stop("NAs cannot be present when trying to find the log-likelihood in this package. Changepoint detection isn't appropriate in the presence of missing values.")
  }

  changepoints <- unique(c(0, changepoints, length(series)))


  #If means and sds for each segment are not already given, find them.
  if(is.na(means) || is.na(sds) || is.na(segLengths)){

    #Record the means and sds associated with each segment
    means <- sds <- segLengths <- rep(NA, length(changepoints) - 1)

    #The mean() and sd() functions are avoided to reduce runtime by not performing error checks.
    #They are therefore written out manually below.
    for(i in 1:(length(changepoints)-1)){
      subSeries <- series[(changepoints[i]+1):changepoints[i+1]]
      segLengths[i] <- changepoints[i+1] - changepoints[i]
      means[i] <- sum(subSeries)/segLengths[i] #mean(subSeries)
      sds[i] <- sqrt((1/(segLengths[i]-1))*sum((subSeries - means[i])^2)) #sd(subSeries)
    }
  }

  #Prepare means and sds by obs
  meansVec <- rep(means, times=segLengths)
  sdsVec <- rep(sds, times=segLengths)

  #Get the logLikelihood of the entire series, ignoring segments that only have one observation (i.e. sd is NA) or are constant (i.e. sd = 0).
  moreThanOne <- (!is.na(sdsVec) & sdsVec != 0)

  #dnorm() is avoided to reduce runtime by not performing error checks.
  #It is therfefore written out manually below.
  #logLik <- sum(dnorm(series[moreThanOne], mean=meansVec[moreThanOne], sd=sdsVec[moreThanOne], log=T))
  logLik <- sum(-(1/2)*log(2*pi) - (1/2)*log(sdsVec[moreThanOne]^2) - (1/(2*sdsVec[moreThanOne]^2))*((series[moreThanOne] - meansVec[moreThanOne])^2))

  return(logLik)
}

#Given a vector of segment lengths, a vector of means, and vector of sds, returns
#sequence randomly generated from corresponding normals for each segment.
generateSample <- function(means, sds, segLengths){
  sample <- c()

  for(i in 1:length(segLengths)){
    sample <- c(sample, rnorm(n=segLengths[i], mean=means[i], sd=sds[i]))
  }

  return(sample)
}


#Takes in a vector of observations, first vector of changepoints, second vector of changepoints
#Get empirical p-value for the osberved improvement in fit when segmenting by the
#two different sets of changepoints.
getPValue <- function(series, changepoints1, changepoints2, numTrials=10000){
  if(any(is.na(series))){
    stop("NAs cannot be present when trying to find the empirical p-value in this package. Changepoint detection isn't appropriate in the presence of missing values.")
  }

  #Check if the previous changepoints are strictly a subset of the new changepoints
  previousAreSubset <- all(changepoints1 %in% changepoints2)

  #Include first and last observations as changepoints.
  changepoints1 <- unique(c(0, changepoints1, length(series)))
  changepoints2 <- unique(c(0, changepoints2, length(series)))

  #If previous changepoints are a subset of the new ones, only simulate segments
  #neighboring the new changepoints for getting the p-value (saves on computation).
  #Else, simulate all segments.
  if(previousAreSubset){
    #Grab the segments to the left and right of each new changepoint.
    newChangepointIndices <- which(!(changepoints2 %in% changepoints1))

    #Subset the series only to the pertinent segments
    #Some computation is purposefully allowed to be redundant. Redundancies could be avoided by
    #checking which segments neighbor one another in the original series, but this would also complicate
    #the code substantially.
    pertinentSeries <- rep(NA, length(series))
    shiftedChangepoints1 <- shiftedChangepoints2 <- c(0)
    for(k in 1:length(newChangepointIndices)){
      changeIndex <- newChangepointIndices[k]

      #Record segment left of changepoint, record shifted changepoint
      pertinentSeries[(changepoints2[changeIndex-1]+1):changepoints2[changeIndex]] <- series[(changepoints2[changeIndex-1]+1):changepoints2[changeIndex]]
      shiftedNewChangepoint <- sum(!is.na(pertinentSeries))
      shiftedChangepoints2 <- c(shiftedChangepoints2, shiftedNewChangepoint)

      #Record segment right of changepoint, record shifted changepoint
      pertinentSeries[(changepoints2[changeIndex]+1):changepoints2[changeIndex+1]] <- series[(changepoints2[changeIndex]+1):changepoints2[changeIndex+1]]
      shiftedNewChangepoint <- sum(!is.na(pertinentSeries))
      shiftedChangepoints2 <- c(shiftedChangepoints2, shiftedNewChangepoint)
      if(changepoints2[changeIndex+1] %in% changepoints1){
        shiftedChangepoints1 <- c(shiftedChangepoints1, shiftedNewChangepoint)
      }
    }

    #Redefine the original series as its pertinent segments and
    #redefine the changepoints accordingly
    series <- pertinentSeries[which(!is.na(pertinentSeries))]
    changepoints1 <- sort(unique(c(0, shiftedChangepoints1, length(series))))
    changepoints2 <- sort(unique(c(0, shiftedChangepoints2, length(series))))
  }

  #Record the means and sds and segments lengths associated with each segment.
  means1 <- sds1 <- segLengths1 <-  rep(NA, length(changepoints1) - 1)
  means2 <- sds2 <- segLengths2 <- rep(NA, length(changepoints2) - 1)

  #Get the means, standard deviations, and segment lengths under the first set of changepoints.
  for(i in 1:(length(changepoints1)-1)){
    means1[i] <- mean(series[(changepoints1[i]+1):changepoints1[i+1]])
    sds1[i] <- sd(series[(changepoints1[i]+1):changepoints1[i+1]])
    segLengths1[i] <- changepoints1[i+1] - changepoints1[i]
  }

  #Get the means, standard deviations, and segment lengths under the second set of changepoints.
  for(i in 1:(length(changepoints2)-1)){
    means2[i] <- mean(series[(changepoints2[i]+1):changepoints2[i+1]])
    sds2[i] <- sd(series[(changepoints2[i]+1):changepoints2[i+1]])
    segLengths2[i] <- changepoints2[i+1] - changepoints2[i]
  }

  #In the case of segments with constant values, replace missing sd with 0.
  sds1[which(is.na(sds1))] <- 0
  sds2[which(is.na(sds2))] <- 0

  #Record the observed change in the log-likelihood.
  obsDiff <- getLogLik(series, changepoints2, means2, sds2, segLengths2) - getLogLik(series, changepoints1, means1, sds1, segLengths1)

  #Generate differences under the null (i.e. the additional changepoints in changepoints2 are not true changepoints).
  nullDiffs <- rep(NA, numTrials)

  #if(serial){
  for(i in 1:numTrials){
    randomSamp <- generateSample(means=means1, sds=sds1, segLengths = segLengths1)
    nullDiffs[i] <- getLogLik(randomSamp, changepoints2) - getLogLik(randomSamp, changepoints1)
  }
  #}

  #Formerly used for parallelism, but no longer used.
  #else{
  #  nullDiffs <- foreach(i=1:numTrials, .combine=c) %dopar% {
  #    randomSamp <- generateSample(means=means1, sds=sds1, segLengths = segLengths1)
  #    getLogLik(randomSamp, changepoints2) - getLogLik(randomSamp, changepoints1)
  #  }
  #}

  return(sum(nullDiffs >= obsDiff)/numTrials)
}


#Given a vector of observations, returns the optimal set of changepoints based on a significance level.
getChangepoints <- function(series, alpha=0.01, numTrials=10000, serial=T, numCores=NA, minPenalty=0, maxPenalty=10e12, verbose=T){
  if(any(is.na(series))){
    stop("Changepoint detection isn't appropriate in the presence of missing values. The series cannot have any NAs.")
  }

  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }

  if(alpha < 0){
    stop("The significance level, alpha, must be a nonnegative value (e.g. 0.01, 0.05, 0.1, etc).")
  }

  if(numTrials < 1){
    stop("numTrials must be a positive integer.")
  }

  #Make cluster if running in parallel.
  if(!serial){
    if(is.na(numCores)){
      cl <- makeCluster(max(1, detectCores()-1))
    }
    else{
      cl <- makeCluster(max(1, numCores))
    }
    registerDoParallel(cl)
    #clusterExport(cl, list("generateSample", "getLogLik"))
  }

  #Run CROPS on PELT to detect changepoints based on changes in mean.
  #capture.output prevents progress messages that are printed by running CROPS.
  capture.output(results <- cpt.mean(series, penalty="CROPS", pen.value=c(minPenalty, maxPenalty), method="PELT", test.stat="Normal", class=F, minseglen=1)$changepoints)
  pValue <- 0

  #Start at end of "results", corresponding to no changepoints. Iterate backwards, including more changepoints.
  maxIndex <- length(results)
  if(maxIndex == 1){
    stop("Only one set of changepoints found. Please try decreasing minPenalty towards 0 and/or increasing maxPenalty.")
  }

  index <- maxIndex

  #Run p-value assessments in serial if requested.
  if(serial){
    while((pValue < alpha) & (index > 1)){
      pValue <- getPValue(series, changepoints1 = results[[index]], changepoints2 = results[[index - 1]], numTrials = numTrials)
      index <- index - 1

      if(verbose){
        if(pValue < alpha){
          if(pValue == 0){
            message(paste("Changepoint set", maxIndex - index, "significant with empirical p-value below:", 1/numTrials))
          }
          else{
            message(paste("Changepoint set", maxIndex - index, "significant with empirical p-value:", pValue))
          }
        }
        else{
          message(paste("Changepoint set", maxIndex - index, "insignificant with empirical p-value:", pValue))
        }
      }

    }

    #If the last set of changepoints is reached without finding an insignificant result, return the last significant result.
    if((pValue < alpha) & (index == 1)){
      if(verbose){
        message("No insignificant set of changepoints considered by PELT. Returning last significant set as final result.")
      }
      return(results[[index]])
    }

    return(results[[index+1]])
  }

  #Run p-value assessments in parallel if requested. Will likely do excessive computation since the stopping point
  #can't be predicted, but should generally still be faster than running in serial.
  else{
    #Track the indices of which set of changepoints should be assessed for significance
    indices <- maxIndex:(maxIndex-numCores+1)
    pValues <- rep(0, length(indices))

    #Get results in parallel for next group of sets of changepoints.
    while((max(pValues) < alpha) & (max(indices) > 1)){

      #Handle whether or not indices goes below 1.
      if(min(indices) <= 1){
        indices <- indices[1]:2
        pValues <- rep(0, length(indices))
      }

      pValues <- foreach(i=indices, .combine=c, .packages = "changepointSelect") %dopar% {
        getPValue(series, changepoints1 = results[[i]], changepoints2 = results[[i - 1]], numTrials = numTrials)
      }


      if(verbose){
        for(i in 1:length(indices)){
          if(pValues[i] < alpha){
            if(pValues[i] == 0){
              message(paste("Changepoint set", maxIndex - indices[i] + 1, "significant with empirical p-value below:", 1/numTrials))
            }
            else{
              message(paste("Changepoint set", maxIndex - indices[i] + 1, "significant with empirical p-value:", pValues[i]))
            }
          }
          else{
            message(paste("Changepoint set", maxIndex - indices[i] + 1, "insignificant with empirical p-value:", pValues[i]))
            break
          }
        }
      }

      #Decrement the indices associated with the sets of changepoints
      indices <- indices - numCores
    }

    #Record which set of changepoints is the first to be insignificant, if such occurs
    if(max(pValues) >= alpha){
      index <- indices[min(which(pValues >= alpha))] + numCores
      stopCluster(cl)
      return(results[[index]])
    }


    #If the last set of changepoints is reached without finding an insignificant result, return the last significant result.
    else{
      if(verbose){
        message("No insignificant set of changepoints considered by PELT. Returning last significant set as final result.")
      }
      stopCluster(cl)
      return(results[[indices[length(indices)] + numCores - 1]])
    }

  }

}
