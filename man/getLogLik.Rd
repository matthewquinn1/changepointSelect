\name{getLogLik}
\alias{getLogLik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the log-likelihood}
\description{
Takes in a vector of observations and a vector of changepoints. It returns the log-likelihood assuming that each segment between two changepoints follows a normal distribution. Optionally, the means, standard deviations, and lengths of each segment can be provided. If not provided, they are calculated.
}
\usage{
getLogLik(series, changepoints, means=NA, sds=NA, segLengths=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which to run changepoint analysis. This vector must not contain any missing values.}
  \item{changepoints}{A vector indicating the indices at which there are changepoints. A changepoint is located immediately prior to a shift in the series.}
  \item{means}{A optional vector where each entry corresponds to the mean of a segment.}
  \item{sds}{A optional vector where each entry corresponds to the standard deviation of a segment.}
  \item{segLengths}{A optional vector where each entry corresponds to the length of a segment.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the log likelihood (a numeric value).
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
