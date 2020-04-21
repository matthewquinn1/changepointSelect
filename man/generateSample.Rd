\name{generateSample}
\alias{generateSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a random sample of observations}
\description{
Generates a random sample of a series of observations. Each segment of the random sample corresponds to sampling from a normal with the given mean and standard deviation.}
\usage{
generateSample(means, sds, segLengths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{means}{A vector where each entry corresponds to the mean of a segment.}
  \item{sds}{A vector where each entry corresponds to the standard deviation of a segment.}
  \item{segLengths}{A vector where each entry corresponds to the length of a segment.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector of the same length as the sum of all segment lengths.
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
