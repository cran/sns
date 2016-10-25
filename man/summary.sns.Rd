\name{summary.sns}
\alias{summary.sns}
\alias{print.summary.sns}

\title{
Summarizing "sns" Objects
}

\description{
Methods for summarizing the output of \code{\link{sns.run}}, and for printing the summary.
}

\usage{
\method{summary}{sns}(object, quantiles = c(0.025, 0.5, 0.975)
  , pval.ref = 0.0, nburnin = max(nrow(object)/2, attr(object, "nnr"))
  , end = nrow(object), thin = 1, ess.method = c("coda", "ise"), ...)
\method{print}{summary.sns}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class "sns", typically the output of \code{\link{sns.run}}.}
  \item{quantiles}{Values for which sample-based quantiles are calculated.}
  \item{pval.ref}{Reference value for state space variables, used for calculating sample-based p-values.}
  \item{nburnin}{Number of initial iterations to discard before calculating the sample statistics. A warning is issued if this number is smaller than the initial iterations run in NR mode.}
  \item{end}{Last iteration to use for calculating sample statistics. Defaults to last iteration.}
  \item{thin}{One out of \code{thin} samples are kept for calculating sample statistics. Default is \code{1}, using all samples within specified range.}
  \item{ess.method}{Method used for calculating effective sample size. Default is to call \code{effectiveSize} from package \code{coda}.}
  \item{x}{An object of class "summary.sns", typically the output of \code{summary.sns}.}
  \item{...}{Arguments passed to/from other functions.}
}

\value{
\code{summary.sns} returns a list with these elements:
  \item{K}{Dimensionality of state space.}
  \item{nnr}{Number of NR (Newton-Raphson) iterations performed at the beginning.}
  \item{nburnin}{Number of burn-in iterations. These are discarded before calculating sample statistics.}
  \item{end}{Last iteration to use for calculating sample statistics.}
  \item{thin}{One out of every \code{thin} iterations within the specified range is used for calculating sample statistics.}
  \item{niter}{Total iterations, including NR and MCMC modes.}
  \item{nsmp}{Number of samples within specified range (before applying thinning).}
  \item{nseq}{Number of samples used for calculating sample statistics (after applying thinning).}
  \item{npart}{Number of subsets used in state space partitioning. If no partitioning is done, the value is \code{1}.}
  \item{accept.rate}{Acceptance rate for the MH transition proposals, calculated over \code{nsmp} iterations.}
  \item{reldev.mean}{Mean relative deviation from quadratic approximation, defined as difference between actual log-density change and the value predicted from quadratic fit at density maximum, divided by the actual change. The location of density maximum is assumed to be the value at the end of the last NR iteration. Therefore, for this measure to be accurate, users must ensure \code{nnr} is sufficiently large to allow for convegrence of the optimization phase.}
  \item{pval.ref}{Same as input.}
  \item{ess.method}{Same as input.}
  \item{smp}{A list with elements \code{mean}, \code{sd}, \code{ess}, \code{quantiles}, \code{pval} representing sample-based mean, standard deviation, effective size, quantiles and sample-based p-values, based on specified range and using thinning (if specified).}
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
