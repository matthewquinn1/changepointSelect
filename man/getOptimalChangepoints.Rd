\name{getOptimalChangepoints}
\alias{getOptimalChangepoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the optimal set of changepoints}
\description{Considers results from running CROPS on PELT in order to assess which set of changepoints is the best. Does so using an empirical p-value process. It is \strong{strongly recommended} that the results from this process are then trimmed for linear trends and/or seasonlity.}
\usage{
getOptimalChangepoints(series, alpha=0.01, adjustAlpha=F, numTrials=10000, serial=T, numCores=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis.}
  \item{alpha}{The significance level to use. Typical values include 0.01, 0.05, and 0.10.}
  \item{adjustAlpha}{If TRUE, uses a bonferroni correction for alpha in order to adjust for the number of new changepoints being assessed. It is currently recommended not to use this adjustment, keeping this FALSE.}
  \item{numTrials}{The number of simulations to use in order to obtain an empirical p-value. Recommended to at least be several thousand.}
  \item{serial}{Boolean indicating whether to run this function in serial or in parallel. Running in parallel will offer only a modest speed-up for short series.}
    \item{numCores}{The number of cores to use if serial=F. If not specified, then it will be taken as max(1, detectCores()-1).}

}
\details{
This function does not consider the possibility that false positives could arise due to trends or seasonality. Therefore, it is strongly recommended that the results from this function are trimmed for false positives. It is recommended to use the \code{\link{trimTrendsAndSeasonality}} function.}
\value{
Returns a vector of changepoints. A changepoint is located immediately prior to a shift in the series.
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
\code{\link{trimLinearTrends}}
\code{\link{trimSeasonality}}
\code{\link{trimTrendsAndSeasonality}}
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
