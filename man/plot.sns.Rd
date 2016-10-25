\name{plot.sns}
\alias{plot.sns}

\title{
Plotting "sns" Objects
}

\description{
Method for visualizing the output of \code{\link{sns.run}}.
}

\usage{
\method{plot}{sns}(x, nburnin = max(nrow(x)/2, attr(x, "nnr"))
  , select = if (length(x) <= 10) 1:5 else 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class "sns", typically the output of \code{\link{sns.run}}.}
  \item{nburnin}{Number of burn-in iterations to discard before generating effective sample size, histograms, and autocorrelation plots.}
  \item{select}{Which plot types must be generated. See below for description.}
  \item{...}{Arguments passed to/from other functions.}
}

\value{
\code{plot.sns} produces the following types of plots: 1) log-probability trace plot (vertical line, if present, indicates transition from nr to mcmc mode), 2) trace plot of state variables (one per coordinate; vertical line has same meaning as 1), 3) effective sample size by coordinate (horizontal line indicates maximum effective size possible, equal to number of samples after discarding nburnin initial iterations), 4) post-burnin state vector histograms (one per coordinate, vertical line indicates post-burnin average, 5) autocorrelation plots, one per coordinate.
}

\references{
Mahani A.S., Hasan A., Jiang M. &  Sharabiani M.T.A. (2016). Stochastic Newton Sampler: The R Package sns. Journal of Statistical Software, Code Snippets, 74(2), 1-33. doi:10.18637/jss.v074.c02
}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}

\seealso{
\code{\link{sns.run}}
}
