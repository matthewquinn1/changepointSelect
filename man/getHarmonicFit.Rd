\name{getHarmonicFit}
\alias{getHarmonicFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the RMSE using harmonic regression}
\description{Internal function that retrieves the model associated with fitting a single harmonic regression over a series of observations.}
\usage{
getHarmonicFit(series, numHarmonics=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis. For this particular function, this is just representative of a single segment or two neighboring segments.}
  \item{numHarmonics}{Indicates the number of harmonics to be modeled (currently restricted to being 1 or 2).}
}
\details{
For all regressions, the dependent variable is the series of observations and the independent variable is the index associated with each observation as part of that series.}
\value{Returns a harmonic regression model resulting from \code{\link[stats]{lm}}.}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
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
