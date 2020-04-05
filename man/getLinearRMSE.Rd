\name{getLinearRMSE}
\alias{getLinearRMSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the RMSE using linear regression}
\description{Internal function that retrieves the RMSE associated with fitting a single linear regression over a series of observations.}
\usage{
getLinearRMSE(series)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis. For this particular function, this is just representative of a single segment or two neighboring segments.}

}
\details{
For all regressions, the dependent variable is the series of observations and the independent variable is the index associated with each observation as part of that series.}
\value{Returns the linear regression RMSE (a numeric value).}
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
