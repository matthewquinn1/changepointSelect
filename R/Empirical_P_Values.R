#Find the log-likelihood for a vector of observations given a vector of changepoints.
getLogLik <- function(series, changepoints, means=NA, sds=NA, segLengths=NA){
  changepoints <- unique(c(0, changepoints, length(series)))


  #If means and sds for each segment are not already given, find them.
  if(is.na(means) || is.na(sds) || is.na(segLengths)){

    #Record the means and sds associated with each segment
    means <- sds <- segLengths <- rep(NA, length(changepoints) - 1)

    for(i in 1:(length(changepoints)-1)){
      subSeries <- series[(changepoints[i]+1):changepoints[i+1]]
      means[i] <- mean(subSeries, na.rm=T)
      sds[i] <- sd(subSeries, na.rm=T)
      segLengths[i] <- changepoints[i+1] - changepoints[i]
    }
  }

  #Prepare means and sds by obs
  meansVec <- rep(means, times=segLengths)
  sdsVec <- rep(sds, times=segLengths)

  #Get the logLikelihood of the entire series, ignoring segments that only have constant segments (i.e. sd is NA)
  nonConstants <- which(!is.na(sdsVec) & sdsVec != 0)
  logLik <- sum(dnorm(series[nonConstants], mean=meansVec[nonConstants], sd=sdsVec[nonConstants], log=T))

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
getPValue <- function(series, changepoints1, changepoints2, numTrials=10000, serial=T){
  #Record the means and sds and segments lengths associated with each segment.
  means1 <- sds1 <- segLengths1 <-  rep(NA, length(changepoints1) + 1)
  means2 <- sds2 <- segLengths2 <- rep(NA, length(changepoints2) + 1)

  #Include first and last observations as changepoints.
  changepoints1 <- unique(c(0, changepoints1, length(series)))
  changepoints2 <- unique(c(0, changepoints2, length(series)))

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

  #If serial processing is requested, run serial. Otherwise, get empirical null distribution in parallel.
  if(serial){
    for(i in 1:numTrials){
      randomSamp <- generateSample(means=means1, sds=sds1, segLengths = segLengths1)
      nullDiffs[i] <- getLogLik(randomSamp, changepoints2) - getLogLik(randomSamp, changepoints1)
    }
  }

  else{
    nullDiffs <- foreach(i=1:numTrials, .combine=c) %dopar% {
      randomSamp <- generateSample(means=means1, sds=sds1, segLengths = segLengths1)
      getLogLik(randomSamp, changepoints2) - getLogLik(randomSamp, changepoints1)
    }
  }

  return(sum(obsDiff < nullDiffs)/numTrials)
}


#Given a vector of observations, returns the optimal set of changepoints based on a significance level.
getChangepoints <- function(series, alpha=0.01, numTrials=10000, serial=T, numCores=NA, verbose=T){
  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }

  if(alpha < 0){
    stop("alpha must be nonnegative.")
  }

  if(numTrials < 1){
    stop("numTrials must be positive.")
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
  capture.output(results <- cpt.mean(series, penalty="CROPS", pen.value=c(0, 10e12), method="PELT", test.stat="Normal", class=F, minseglen=1)$changepoints)
  pValue <- 0

  #Start at end of "results", corresponding to no changepoints. Iterate backwards, including more changepoints.
  maxIndex <- length(results)
  index <- maxIndex

  while((pValue < alpha) & (index > 1)){
    pValue <- getPValue(series, changepoints1 = results[[index]], changepoints2 = results[[index - 1]], numTrials, serial=serial)
    index <- index - 1

    if(verbose){
      if(pValue < alpha){
        message(paste("Changepoint set", maxIndex - index, "significant with empirical p-value:", pValue))
      }
      else{
        message(paste("Changepoint set", maxIndex - index, "insignificant with empirical p-value:", pValue))
      }
    }

  }

  if(!serial){
    stopCluster(cl)
  }

  return(results[[index+1]])
}
