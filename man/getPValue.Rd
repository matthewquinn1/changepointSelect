\name{getPValue}
\alias{getPValue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the empirical p-value for a new set of changepoints}
\description{
Obtains the p-value associated with introducing new changepoints beyond those already present for a series of observations.}
\usage{
getPValue(series, changepoints1, changepoints2, numTrials=10000, serial=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis.}
  \item{changepoints1}{The initial vector of changepoints. A changepoint is located immediately prior to a shift in the series.}
  \item{changepoints2}{A vector of changepoints to assess the significance of, relative to the initial vector.}
  \item{numTrials}{The number of simulations to use in order to obtain an empirical p-value. Recommended to at least be several thousand.}
  \item{serial}{Boolean indicating whether to run this function in serial or in parallel. Running in parallel will offer only a modest speed-up for short series and/or for relatively small numbers of simulations.}
}
\details{
This function is called by \code{\link{getChangepoints}}. If running in parallel, the number of cores to use is set by \code{\link{getChangepoints}}, which sets up the cluster.}
\value{
Returns the empirical p-value (a numeric value between 0 and 1 inclusive).
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
