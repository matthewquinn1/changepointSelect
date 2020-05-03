\name{trimLinearTrends}
\alias{trimLinearTrends}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trim proposed changepoints arising from linear trends}
\description{
Removes proposed changepoints that appear to be false positives due to linear trends.}
\usage{
trimLinearTrends(series, changepoints, threshold=1.15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis.}
  \item{changepoints}{A vector indicating the indices at which there are proposed changepoints. A changepoint is located immediately prior to a shift in the series.}
  \item{threshold}{Parameter for determining when to trim changepoints. If the ratio of the linear regression RMSE to the optimal piecewise RMSE is below this threshold, then the changepoint is subject to removal. That is, a larger threshold makes it easier to remove proposed changepoints. Reasonable values may be between 1 and 1.5, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold.}

}
\details{
Iterates over proposed changepoints. For a given changepoint, the two neighboring segments on either side of the changepoint are considered. If the ratio of the linear regression RMSE (across both segments together) to the optimal piecewise linear regression RMSE (each segment separately) is below this threshold, then the changepoint is subject to removal. That is, a larger threshold makes it easier to remove proposed changepoints.

After iterating over all changepoints, the changepoint with the smallest ratio (i.e. closest to 1) is removed and the process repeats until no changepoints are subject to removal.}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Obtain changepoints for the simulated data. 10000 simulation trials may take some time.
simChangepoints <- getChangepoints(series=simSeries, alpha=0.01, numTrials=10000, serial=F, numCores = 2)

#Trim changepoints for linear trends
simChangepointsTrimmed <- trimLinearTrends(series=simSeries, changepoints=simChangepoints, threshold=1.2)

#Plot
plot(y=simSeries, x=1:length(simSeries), type="l", ylab="Simulated Value", xlab="Index", main="Simulated Data with Trimmed Changepoints")
abline(v=simChangepointsTrimmed, col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
