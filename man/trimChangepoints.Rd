\name{trimChangepoints}
\alias{trimChangepoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trim proposed changepoints arising from linear trends or seasonality
}
\description{
Wrapper function that will call \code{\link{trimLinearTrends}}, \code{\link{trimSeasonality}}, or \code{\link{trimLinearTrendsAndSeasonality}} to remove proposed changepoints that appear to be false positives due to linear trends and/or seasonality.
}
\usage{
trimChangepoints(series, changepoints, thresholdLinear=1.15, thresholdSeasonal=1.15, numHarmonics=2, verbose=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis.}
  \item{changepoints}{A vector indicating the indices at which there are proposed changepoints. A changepoint is located immediately prior to a shift in the series.}
  \item{thresholdLinear}{Parameter for determining when to trim changepoints based on linear trends. If the ratio of the linear regression RMSE to the optimal piecewise RMSE is below this threshold, then the changepoint is subject to removal. That is, a larger threshold makes it easier to remove proposed changepoints. Reasonable values may be between 1 and 1.5, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold. It is recommended to use the same threshold as for thresholdSeasonal. If this parameter is less than or equal to 1, then trimming for linear trends is not performed.}
    \item{thresholdSeasonal}{Parameter for determining when to trim changepoints based on seasonality. If the ratio of the harmonic regression RMSE to the optimal piecewise RMSE is below this threshold, then the changepoint is subject to removal. That is, a larger threshold makes it easier to remove proposed changepoints. Reasonable values may be between 1 and 1.5, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold. It is recommended to use the same value as used for thresholdLinear. If this parameter is less than or equal to 1, then trimming for seasonality is not performed.}
  \item{numHarmonics}{Indicates the number of harmonics to be modeled when assessing seasonality (currently restricted to being 1 or 2).}
  \item{verbose}{If TRUE, prints out message indicating status.}
}
\details{In order to trim changepoints for both linear trends and seasonality, \code{thresholdLinear} and \code{thresholdSeaonal} must both be larger than 1. If  \code{thresholdLinear > 1} and  \code{thresholdSeasonal <= 1}, then the changepoints will only be trimmed for linear trends. If  \code{thresholdLinear <= 1} and  \code{thresholdSeasonal > 1}, then the changepoints will only be trimmed for seasonality. However, it is generally recommended to trim for both linear trends and seasonality.
}
\value{
Returns a vector of changepoints that are retained. This will be the subset of the input changepoints vector, but excluding those that were trimmed.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Matthew Quinn}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getChangepoints}}
}
\examples{
#Obtain changepoints for the simulated data. 10000 simulation trials may take some time.
simChangepoints <- getChangepoints(series=simSeries, alpha=0.01, numTrials=10000, serial=F, numCores = 2)

#Trim changepoints for linear trends and seasonality
simChangepointsTrimmed <- trimChangepoints(series=simSeries, changepoints=simChangepoints, thresholdLinear=1.2, thresholdSeasonal=1.2)

#Plot
plot(y=simSeries, x=1:length(simSeries), type="l", ylab="Simulated Value", xlab="Index", main="Simulated Data with Trimmed Changepoints")
abline(v=simChangepointsTrimmed, col="red")


#Trim changepoints for linear trends only
simChangepointsTrimmed <- trimChangepoints(series=simSeries, changepoints=simChangepoints, thresholdLinear=1.2, thresholdSeasonal=1)

#Plot
plot(y=simSeries, x=1:length(simSeries), type="l", ylab="Simulated Value", xlab="Index", main="Simulated Data with Trimmed Changepoints")
abline(v=simChangepointsTrimmed, col="red")

#Trim changepoints for seasonality only
simChangepointsTrimmed <- trimChangepoints(series=simSeries, changepoints=simChangepoints, thresholdLinear=1, thresholdSeasonal=1.2)

#Plot
plot(y=simSeries, x=1:length(simSeries), type="l", ylab="Simulated Value", xlab="Index", main="Simulated Data with Trimmed Changepoints")
abline(v=simChangepointsTrimmed, col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
