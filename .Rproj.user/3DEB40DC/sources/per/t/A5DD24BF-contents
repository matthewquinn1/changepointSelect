\name{trimTrendsAndSeasonality}
\alias{trimTrendsAndSeasonality}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trim proposed changepoints arising from linear trends or seasonality
}
\description{
Removes proposed changepoints that appear to be false positives due to linear trends and/or seasonality.
}
\usage{
trimTrendsAndSeasonality(series, changepoints, thresholdLinear=1.25, thresholdSeasonal=1.25, numHarmonics=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis.}
  \item{changepoints}{A vector indicating the indices at which there are proposed changepoints. A changepoint is located immediately prior to a shift in the series.}
  \item{thresholdLinear}{Parameter for determining when to trim changepoints based on linear trends. If the ratio of the linear regression RMSE to the optimal piecewise RMSE is below this threshold, then the changepoint is subject to removal. That is, a larger threshold makes it easier to remove proposed changepoints. Reasonable values may be between 1 and 1.5, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold. It is recommended to use the same threshold as for thresholdSeasonal.}
    \item{thresholdSeasonal}{Parameter for determining when to trim changepoints based on seasonality. If the ratio of the harmonic regression RMSE to the optimal piecewise RMSE is below this threshold, then the changepoint is subject to removal. That is, a larger threshold makes it easier to remove proposed changepoints. Reasonable values may be between 1 and 1.5, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold. It is recommended to use the same value as used for thresholdLinear.}
  \item{numHarmonics}{Indicates the number of harmonics to be modeled when assessing seasonality (currently restricted to being 1 or 2).}
}
\details{
Iterates over proposed changepoints. For a given changepoint, the two neighboring segments on either side of the changepoint are considered. If the ratio of the linear regression RMSE (across both segments together) to the optimal piecewise regression RMSE (each segment separately) is below thresholdLinear, then the changepoint is subject to removal. If the ratio of the harmonic regression RMSE (across both segments together) to the optimal piecewise regression RMSE (each segment separately) is below thresholdSeasonal, then the changepoint is subject to removal.  That is, a larger threshold makes it easier to remove proposed changepoints.

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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
