\name{basicSegmentCorrection}
\alias{basicSegmentCorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs a basic segment correction.}

\description{Taking in a vector of observations and at least a list of best fitting models, this function will iterate over segments and perform de-trending or de-seasonlizing, followed by shifting and scaling to a reference. This is simply a tool for preliminary analysis after identifying changepoints and should not replace subject-matter knowledge, formal model building, and correction for changepoints.}
\usage{
basicSegmentCorrection(series, changepoints = c(), bestFitsBySegment, segmentsToCorrect = c(), referenceSegment = NA, verbose=T, plotResults=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which changepoint analysis has been run. This vector must not contain any missing values.}
  \item{changepoints}{A vector indicating the indices at which there are proposed changepoints. A changepoint is located immediately prior to a shift in the series. If left as an empty vector (the default), the changepoints will be detected from \code{bestFitsBySegment} based on the end of each segment.}
   \item{bestFitsBySegment}{A list of lists with each sublist corresponding to a segment. The structure of this sublist should match the formatting of the results from \code{\link{getBestFitBySegment}}.}
    \item{segmentsToCorrect}{A vector denoting which segments to perform de-trending, de-seasonlizing, shifting, and scaling on. If left an empty vector (default), all segments will be corrected.}
    \item{referenceSegment}{The results in each corrected segment will separately be shifted and scaled to match the mean and standard deviation of this segment. If left \code{NA} (the default), the largest segment in \code{series} is chosen. If set to \code{0}, then the observations will be standardized to have mean 0 and standard deviation 1. This process is performed after any de-trending and de-seasonalizing. The reference mean will be taken as the mean of \code{referenceSegment} before de-trending or de-seasonlizing while the reference standard deviation will be estimated based on the residual standard error of the best fitting model for \code{referenceSegment}, even if \code{referenceSegment} is not in \code{segmentsToCorrect}. This shifting and scaling is done only to segments in \code{segmentsToCorrect}. Other segments are left unaltered}
    \item{verbose}{If \code{TRUE}, prints out messages indicating progress.}
    \item{plotResults}{If \code{TRUE}, produces a plot of the time series after correction.}
}
\details{
This function should only be run after \code{\link{getBestFitBySegment}}. It iterates over segments. If the segment's best fitting model is a linear trend, the segment is de-trended (i.e. the segment is replaced by the residuals from the linear model). If the best fitting model is a harmonic regression, the segment is de-seasonalized (i.e. the segment is replaced by the residuals from the harmonic regression). All segments are then shifted and scaled to match the mean and standard deviation of the \code{referenceSegment}. The reference mean will be taken as the mean of \code{referenceSegment} before de-trending or de-seasonlizing while the reference standard deviation will be estimated based on the residual standard error of the best fitting model for \code{referenceSegment}, even if \code{referenceSegment} is not in \code{segmentsToCorrect}. The standard deviation is estimated using residual standard error rather than the sample standard deviation of the original segment to avoid having trends and seasonality influence this estimate.

In practice, one should not use this function to replace subject-matter knowledge, formal model building, and correction for changepoints. This is simply an additional tool that provides a preliminary analysis of segments resulting from identified changepoints.
}
\value{
Returns a vector containing the segment-corrected observations.
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
