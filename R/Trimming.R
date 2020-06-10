#Given a series of observations and set of changepoints, trim them in the event of seasonality or trends.

getPiecewiseRMSE <- function(series, changepoint, trendsOnly=F, seasonalityOnly=F, numHarmonics=2){
  #Retrieve the best piecewise RMSE based on looking for trends and/or seasonality
  if(trendsOnly){
    RMSE1 <- getLinearRMSE(series[1:changepoint])
    RMSE2 <- getLinearRMSE(series[(changepoint+1):length(series)])
  }

  if(seasonalityOnly){
    RMSE1 <- getHarmonicRMSE(series[1:changepoint], numHarmonics)
    RMSE2 <- getHarmonicRMSE(series[(changepoint+1):length(series)], numHarmonics)
  }

  # if(trendedSeasonality){
  #   RMSE1 <- getTrendedHarmonicRMSE(series[1:changepoint], numHarmonics)
  #   RMSE2 <- getTrendedHarmonicRMSE(series[(changepoint+1):length(series)], numHarmonics)
  # }

  #Fit both linear and harmonic regressions to each subsegment, recording the piecewise RMSE as the minimum
  else{
    RMSE1 <- min(getLinearRMSE(series[1:changepoint]), getHarmonicRMSE(series[1:changepoint], numHarmonics))
    RMSE2 <- min(getLinearRMSE(series[(changepoint+1):length(series)]), getHarmonicRMSE(series[(changepoint+1):length(series)], numHarmonics))
  }

  return(sqrt((1/length(series))*(changepoint*RMSE1^2 + (length(series)-changepoint)*RMSE2^2)))
}

getLinearRMSE <- function(series){
  if(length(series) < 3){
    return(0)
  }

  x <- 1:length(series)

  fitLinear <- lm(series ~ x)
  return(sqrt(mean(resid(fitLinear)^2)))
}

getHarmonicRMSE <- function(series, numHarmonics=2){
  if(length(series) < 3){
    return(0)
  }

  x <- 1:length(series)

  #General description of spectral analysis/periodogram
  #http://www.ams.sunysb.edu/~zhu/ams586/Periodogram.pdf

  #Fit model with first two harmonics
  #https://stats.stackexchange.com/questions/60500/how-to-find-a-good-fit-for-semi-sinusoidal-model-in-r
  #https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
  ssp <- spectrum(series, plot=F, method="pgram")
  per <- 1/ssp$freq[ssp$spec==max(ssp$spec)][1]

  if(numHarmonics==1){
    fitHarmonic <- lm(series~ sin(2*pi*x/per)+cos(2*pi*x/per))
  }

  else{
    fitHarmonic <- lm(series~ sin(2*pi*x/per)+cos(2*pi*x/per) + sin(4*pi*x/per)+cos(4*pi*x/per))
  }
  return(sqrt(mean(resid(fitHarmonic)^2)))
}

# getTrendedHarmonicRMSE <- function(series, numHarmonics=2){
#   if(length(series) < 3){
#     return(0)
#   }
#
#   x <- 1:length(series)
#
#   #Fit model with first two harmonics
#   #https://stats.stackexchange.com/questions/60500/how-to-find-a-good-fit-for-semi-sinusoidal-model-in-r
#   #https://stats.stackexchange.com/questions/60994/fit-a-sinusoidal-term-to-data
#   ssp <- spectrum(series, plot=F, method="pgram")
#   per <- 1/ssp$freq[ssp$spec==max(ssp$spec)][1]
#
#
#   if(numHarmonics==1){
#     fitHarmonic <- lm(series~ x + sin(2*pi*x/per)+cos(2*pi*x/per))
#   }
#
#   else{
#     fitHarmonic <- lm(series~x + sin(2*pi*x/per)+cos(2*pi*x/per) + sin(4*pi*x/per)+cos(4*pi*x/per))
#   }
#
#   return(sqrt(mean(resid(fitHarmonic)^2)))
# }

trimLinearTrends <- function(series, changepoints, threshold=1.15){
  if(threshold < 1){
    stop("The trimming threshold should be no smaller than 1.")
  }

  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }

  if(length(changepoints) < 1){
    stop("The changepoints must be in a vector and there must be at least 1 changepoint.")
  }

  #Include 0 and n as endpoints.
  changepoints <- unique(c(0, changepoints, length(series)))

  #Iteratively combine segments if a linear regression RMSE is no more than threshold times that for piecewise fits.
  while(length(changepoints) > 2){

    linearResults <- piecewiseResults <- rep(NA, length(changepoints) - 2)

    #Get the overall linear fit and piecewise fit RMSEs corresponding to each changepoint.
    for(i in 1:(length(changepoints)-2)){
      linearResults[i] <- getLinearRMSE(series[(changepoints[i]+1):changepoints[i+2]])
      piecewiseResults[i] <- getPiecewiseRMSE(series[(changepoints[i]+1):changepoints[i+2]], changepoint=changepoints[i+1]-changepoints[i], trendsOnly = T)
    }

    ratios <- linearResults/piecewiseResults

    #Once no linear RMSE is less than threshold times that for piecewise functions, stop.
    if(min(ratios) > threshold){
      return(changepoints[2:(length(changepoints)-1)])
    }

    #Remove changepoint for which linear regression is smallest improvement over piecewise fits.
    removeInd <- which.min(ratios) + 1
    changepoints <- changepoints[-removeInd]
  }

  if(length(changepoints) > 2){
    return(changepoints[2:(length(changepoints)-1)])
  }
  return(c())
}


trimSeasonality <- function(series, changepoints, threshold=1.15, numHarmonics=2){
  if(threshold < 1){
    stop("The trimming threshold should be no smaller than 1.")
  }

  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }

  if(length(changepoints) < 1){
    stop("The changepoints must be in a vector and there must be at least 1 changepoint.")
  }

  #Include 0 and n as endpoints.
  changepoints <- unique(c(0, changepoints, length(series)))

  #Iteratively combine segments if a harmonic regression RMSE is no more than threshold times that for piecewise fits.
  while(length(changepoints) > 2){

    seasonalResults <- piecewiseResults <- rep(NA, length(changepoints) - 2)

    #Get the overall harmonic fit and piecewise fit RMSEs corresponding to each changepoint.
    for(i in 1:(length(changepoints)-2)){
      seasonalResults[i] <- getHarmonicRMSE(series[(changepoints[i]+1):changepoints[i+2]], numHarmonics)
      piecewiseResults[i] <- getPiecewiseRMSE(series[(changepoints[i]+1):changepoints[i+2]], changepoint=changepoints[i+1]-changepoints[i], seasonalityOnly = T, numHarmonics)
    }

    ratios <- seasonalResults/piecewiseResults

    #Once no harmonic RMSE is less than threshold times that for piecewise functions, stop.
    if(min(ratios) > threshold){
      return(changepoints[2:(length(changepoints)-1)])
    }

    #Remove changepoint for which harmonic regression is smallest improvement over piecewise fits.
    removeInd <- which.min(ratios) + 1
    changepoints <- changepoints[-removeInd]
  }

  if(length(changepoints) > 2){
    return(changepoints[2:(length(changepoints)-1)])
  }
  return(c())
}


# trimTrendedSeasonality <- function(series, changepoints, threshold=1.15, numHarmonics=2){
#   #Include 0 and n as endpoints.
#   changepoints <- unique(c(0, changepoints, length(series)))
#
#   #Iteratively combine segments if a trended harmonic regression RMSE is no more than threshold times that for piecewise fits.
#   while(length(changepoints) > 2){
#
#     trendedSeasonalResults <- piecewiseResults <- rep(NA, length(changepoints) - 2)
#
#     #Get the overall trended harmonic fit and piecewise fit RMSEs corresponding to each changepoint.
#     for(i in 1:(length(changepoints)-2)){
#       trendedSeasonalResults[i] <- getTrendedHarmonicRMSE(series[(changepoints[i]+1):changepoints[i+2]], numHarmonics)
#       piecewiseResults[i] <- getPiecewiseRMSE(series[(changepoints[i]+1):changepoints[i+2]], changepoint=changepoints[i+1]-changepoints[i], trendedSeasonality = T, numHarmonics)
#     }
#
#     ratios <- trendedSeasonalResults/piecewiseResults
#
#     #Once no trended harmonic RMSE is less than threshold times that for piecewise functions, stop.
#     if(min(ratios) > threshold){
#       return(changepoints[2:(length(changepoints)-1)])
#     }
#
#     #Remove changepoint for which trended harmonic regression is smallest improvement over piecewise fits.
#     removeInd <- which.min(ratios) + 1
#     changepoints <- changepoints[-removeInd]
#   }
#
#   if(length(changepoints) > 2){
#     return(changepoints[2:(length(changepoints)-1)])
#   }
#   return(c())
# }



trimLinearTrendsAndSeasonality <- function(series, changepoints, thresholdLinear=1.15, thresholdSeasonal=1.15, numHarmonics=2){
  if((thresholdLinear < 1) | (thresholdSeasonal < 1)){
    stop("The trimming threshold should be no smaller than 1.")
  }

  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }

  if(length(changepoints) < 1){
    stop("The changepoints must be in a vector and there must be at least 1 changepoint.")
  }

  #Include 0 and n as endpoints.
  changepoints <- unique(c(0, changepoints, length(series)))

  #Iteratively combine segments if a linear regression or harmonic regression RMSE is no more than threshold times that for piecewise fits.
  while(length(changepoints) > 2){

    linearResults <- seasonalResults <- piecewiseResults <- rep(NA, length(changepoints) - 2)

    #Get the overall linear fit, harmonic fit, and piecewise fit RMSEs corresponding to each changepoint.
    for(i in 1:(length(changepoints)-2)){
      linearResults[i] <- getLinearRMSE(series[(changepoints[i]+1):changepoints[i+2]])
      seasonalResults[i] <- getHarmonicRMSE(series[(changepoints[i]+1):changepoints[i+2]], numHarmonics)
      piecewiseResults[i] <- getPiecewiseRMSE(series[(changepoints[i]+1):changepoints[i+2]], changepoint=changepoints[i+1]-changepoints[i], numHarmonics)
    }

    ratiosLinear <- linearResults/piecewiseResults
    ratiosSeasonal <- seasonalResults/piecewiseResults
    ratios <- pmin(ratiosLinear, ratiosSeasonal) #Take smallest ratio across linear trends and seasonality


    #Once no linear RMSE or harmonic RMSE is less than threshold times that for piecewise functions, stop.
      if((min(ratiosLinear) > thresholdLinear) && (min(ratiosSeasonal) > thresholdSeasonal)){
        return(changepoints[2:(length(changepoints)-1)])
      }

      #Remove changepoint for which linear regression or harmonic regression is smallest improvement over piecewise fits.
      removeInd <- which.min(ratios) + 1
      changepoints <- changepoints[-removeInd]

    #When using dynamic trimming.
    # else{
    #   segLengths <- diff(changepoints)
    #   thresholds <- 1 + dynamicMultiplier*runmin(segLengths,k=2)/runmax(segLengths,k=2)
    #
    #   #Once no ratio of RMSEs is below a threshold, stop.
    #   if(all(ratios > thresholds)){
    #     return(changepoints[2:(length(changepoints)-1)])
    #   }
    #
    #   #Remove changepoint for which linear regression or harmonic regression is smallest improvement over piecewise fits.
    #   belowThreshold <- which(ratios <= thresholds)
    #   removeInd <- belowThreshold[which.min(ratios[belowThreshold]/thresholds[belowThreshold])] + 1
    #   changepoints <- changepoints[-removeInd]
    # }

  }

  if(length(changepoints) > 2){
    return(changepoints[2:(length(changepoints)-1)])
  }
  return(c())
}


#Wrapper for other trimming functions
trimChangepoints <- function(series, changepoints, thresholdLinear=1.15, thresholdSeasonal=1.15, numHarmonics=2, verbose=T){
  if((thresholdLinear <= 1) & (thresholdSeasonal <= 1)){
    stop("Both trimming thresholds are no more than 1. No trimming done.")
  }
  if(length(series) < 1){
    stop("The series must be a vector containing at least 1 observation.")
  }
  if(length(changepoints) < 1){
    stop("The changepoints must be in a vector and there must be at least 1 changepoint.")
  }

  #Call other trimming functions
  if((thresholdLinear > 1) & (thresholdSeasonal > 1)){
    if(verbose) {message("Trimming for both linear trends and seasonality.")}
    return(trimLinearTrendsAndSeasonality(series, changepoints, thresholdLinear, thresholdSeasonal, numHarmonics))
  }

  if((thresholdLinear > 1) & (thresholdSeasonal <= 1)){
    if(verbose) {message("Trimming for linear trends only.")}
    return(trimLinearTrends(series, changepoints, thresholdLinear))
  }

  if((thresholdLinear <= 1) & (thresholdSeasonal > 1)){
    if(verbose) {message("Trimming for seasonality only.")}
    return(trimSeasonality(series, changepoints, thresholdSeasonal, numHarmonics))
  }
}

#Helper functions for trimTrendsAndSeasonality if using dynamic trimming
# runmin <- function(values, window){
#   mins <- rep(NA, length(values)-window+1)
#
#   for(i in 1:(length(values)-window+1)){
#     mins[i] <- min(values[i:(i+window-1)])
#   }
#
#   return(mins)
# }
#
# runmax <- function(values, window){
#   maxes <- rep(NA, length(values)-window+1)
#
#   for(i in 1:(length(values)-window+1)){
#     maxes[i] <- max(values[i:(i+window-1)])
#   }
#
#   return(maxes)
# }
