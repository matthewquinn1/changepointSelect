\name{getBestFitBySegment}
\alias{getBestFitBySegment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Preliminarily find the most appropriate model (constant, linear, harmonic) for each segment}

\description{Taking in a vector of observations and a vector of changepoints, this function returns whether it is most appropriate to model each segment with a constant function (i.e. the mean of the segment), a linear regression, or a harmonic regression. This is simply a tool for preliminary analysis after identifying changepoints and should not replace subject-matter knowledge, formal model building, and correction for changepoints.}
\usage{
getBestFitBySegment(series, changepoints, thresholdLinear=1.25, thresholdSeasonal=1.25, numHarmonics=2, verbose=T, plotFits=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{A vector of observations on which changepoint analysis has been run. This vector must not contain any missing values.}
  \item{changepoints}{A vector indicating the indices at which there are proposed changepoints. A changepoint is located immediately prior to a shift in the series.}
   \item{thresholdLinear}{Parameter for determining when to adopt a linear trend for a segment over a constant function (i.e. the mean of the segment). If the ratio of the constant function's RMSE to the linear regression's RMSE is below this threshold, then the constant function is adopted. That is, a larger threshold makes it easier to adopt the constant function. Reasonable values may be between 1 and 2, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold. If this parameter is less than or equal to 1, then the constant function will not be adopted except for very short segments.}
    \item{thresholdSeasonal}{Parameter for determining when to adopt a harmonic model for a segment over a constant function (i.e. the mean of the segment). If the ratio of the constant function's RMSE to the harmonic regression's RMSE is below this threshold, then the constant function is adopted. That is, a larger threshold makes it easier to adopt the constant function. Reasonable values may be between 1 and 2, but will depend on the context. It is worthwhile for the investigator to try multiple values for the threshold. If this parameter is less than or equal to 1, then the constant function will not be adopted except for very short segments.}
    \item{numHarmonics}{Indicates the number of harmonics to be modeled when assessing seasonality (currently restricted to being 1 or 2).}
    \item{verbose}{If \code{TRUE}, prints out messages indicating progress.}
    \item{plotFits}{If \code{TRUE}, produces a plot of the time series with the best fit by segment.}
}
\details{
For each segment, the constant function (i.e. the mean of the segment), a linear regression, and a harmonic regression are fit. The RMSEs are recorded for each model. Among the linear regression and harmonic regression, the one with the smaller RMSE is considered a candidate model and compared to a constant fit. If the ratio of the constant fit's RMSE to the candidate model's RMSE is less than the corresponding threshold, then the constant model is adopted over the candidate model. Otherwise, if the ratio is at least as large as the threshold, then the candidate model is adopted.

For a linear regression to be fit to a segment, this function requires that the segment contain at least 3 points. If a harmonic regression is to be fit, this function requires the segment to have at least 2 + 2*\code{numHarmonics} points.

Please note that this function makes comparisons between constant, linear, and harmonic models for on segment at a time. This is very different from \code{\link{trimChangepoints}}, which compares piecewise models to cross-segment models for two neighborhing segments at a time.

Also note that this function does not consider anything other than the values within a segment when fitting a model. In practice, one should not use this function to replace subject-matter knowledge, formal model building, and correction for changepoints. This is simply an additional tool that provides a preliminary analysis of segments resulting from identified changepoints.
}
\value{
Returns a list containing a sublist for each segment. Each sublist contains the starting index (inclusive) of the segment, ending index (inclusive) of the segment, the type of model that is most appropriate, and the \code{\link[stats]{lm}} object describing that model.
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
