\name{getPiecewiseRMSE}
\alias{getPiecewiseRMSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the RMSE using piecewise regressions}
\description{
Internal function that retrieves the optimal RMSE from performing piecewise linear regression and/or harmonic regression fits on two consecutive segments.}
\usage{
getPiecewiseRMSE(series, changepoint, trendsOnly=F, seasonalityOnly=F, numHarmonics=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis. For this particular function, this series should only contain two segments, separated by a single changepoint.}
  \item{changepoint}{The index indicating the changepoint. A changepoint is located immediately prior to a shift in the series.}
  \item{trendsOnly}{Boolean that, when TRUE, indicates that only linear regressions should be used.}
  \item{seasonalityOnly}{Boolean that, when TRUE, indicates that only harmonic regressions should be used.}
  \item{numHarmonics}{If harmonic regressions are used, this indicates the number of harmonics to be modeled (currently restricted to being 1 or 2).}
}
\details{
When both trendsOnly and seasonalityOnly are FALSE (the default), then both linear regressions and harmonic regressions are considered for the piecewise fits. The best fits on either side of the changepoint contribute to the piecewise RMSE.

For all regressions, the dependent variable is the series of observations and the independent variable is the index associated with each observation as part of that series.}
\value{
Returns the optimal piecewise RMSE (a numeric value).
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
