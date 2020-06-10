\name{getChangepoints}
\alias{getChangepoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find the optimal set of changepoints based on a significance level}
\description{Considers results from running CROPS on PELT in order to assess which set of changepoints is the best. Does so using an empirical p-value process. It is \strong{strongly recommended} that the results from this process are then trimmed for linear trends and/or seasonlity.}
\usage{
getChangepoints(series, alpha=0.01, numTrials=10000, serial=T, numCores=NA,  minPenalty=0, maxPenalty=10e12, verbose=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis.}
  \item{alpha}{The significance level to use. Typical values include 0.01, 0.05, and 0.10.}
  \item{numTrials}{The number of simulations to use in order to obtain an empirical p-value. Recommended to at least 1,000.}
  \item{serial}{Boolean indicating whether to run this function in serial or in parallel. Running in parallel will offer only a modest speed-up for short series or when \code{numTrials} is relatively small.}
    \item{numCores}{The number of cores to use if serial=F. If not specified, then it will be taken as \code{max(1, detectCores()-1)}.}
    \item{minPenalty}{The minimum penalty to be used when running CROPS on PELT. This corresponds to the first of two values fed to the \code{pen.value} argument in \code{\link[changepoint]{cpt.mean}}. It is suggested to keep this set to 0.}
    \item{maxPenalty}{The maximum penalty to be used when running CROPS on PELT. This corresponds to the second of two values fed to the \code{pen.value} argument in \code{\link[changepoint]{cpt.mean}}. This can be made arbitrarily large with relatively little impact on runtime. It is important that this value is large enough that no changepoints are identified under such a penalty in PELT. In rare cases, it might be necessary to increase this beyond the default value in order to get valid results.}
    \item{verbose}{If TRUE, prints out messages indicating progress.}

}
\details{
This function does not consider the possibility that false positives could arise due to trends or seasonality. Therefore, it is strongly recommended that the results from this function are trimmed for false positives. It is recommended to use the \code{\link{trimChangepoints}} function.}
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
\code{\link{trimChangepoints}}
}

\examples{
#Obtain changepoints for the simulated data. 10000 simulation trials may take some time.
simChangepoints <- getChangepoints(series=simSeries, alpha=0.01, numTrials=10000, serial=F, numCores = 2)

#Must then go on to trim changepoints.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
