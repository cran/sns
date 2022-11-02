\name{sns.run}
\alias{sns.run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Drawing multiple samples using Stochastic Newton Sampler
}
\description{
This is a wrapper around \code{sns}, allowing one to draw multiple samples from a distribution while collecting diagnostic information.
}
\usage{
sns.run(init, fghEval, niter = 100, nnr = min(10, round(niter/4))
  , mh.diag = FALSE, part = NULL, print.level = 0
  , report.progress = ceiling(niter/10)
  , numderiv = 0, numderiv.method = c("Richardson", "simple")
  , numderiv.args = list()
  , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{Initial value for the MCMC chain.}
  \item{fghEval}{Log-density to be sampled from. A valid log-density can have one of 3 forms: 1) return log-density, but no gradient or Hessian, 2) return a list of \code{f} and \code{g} for log-density and its gradient vector, respectively, 3) return a list of \code{f}, \code{g}, and \code{h} for log-density, gradient vector, and Hessian matrix. Missing derivatives are computed numerically.}
  \item{niter}{Number of iterations to perform (in `nr' and `mcmc' mode combined).}
  \item{nnr}{Number of initial iterations to spend in `nr' mode.}
  \item{mh.diag}{Boolean flag, indicating whether detailed MH diagnostics such as components of acceptance test must be returned or not.}
  \item{part}{List describing partitioning of state space into subsets. Each element of the list must be an integer vector containing a set of indexes (between \code{1} and \code{length(x)} or \code{length(init)}) indicating which subset of all dimensions to jointly sample. These integer vectors must be mutually exclusive and collectively exhaustive, i.e. cover the entire state space and have no duplicates, in order for the partitioning to represent a valid Gibbs sampling approach. See \code{sns.make.part} and \code{sns.check.part}.}
  \item{print.level}{If greater than 0, print sampling progress report.}
  \item{report.progress}{Number of sampling iterations to wait before printing progress reports.}
  \item{numderiv}{Integer with value from the set \code{0,1,2}. If \code{0}, no numerical differentiation is performed, and thus \code{fghEval} is expected to supply \code{f}, \code{g} and \code{h}. If \code{1}, we expect \code{fghEval} to provide \code{f} amd \code{g}, and Hessian will be calculated numerically. If \code{2}, \code{fghEval} only returns log-density, and numerical differentiation is needed to calculate gradient and Hessian.}
  \item{numderiv.method}{Method used for numeric differentiation. This is passed to the \code{grad} and \code{hessian} functions in \pkg{numDeriv} package. See the package documentation for details.}
  \item{numderiv.args}{Arguments to the numeric differentiation method chosen in \code{numderiv.method}, passed to \code{grad} and \code{hessian} functions in \pkg{numDeriv}. See package documentation for details.}
  \item{\dots}{Other parameters to be passed to \code{fghEval}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\code{sns.run} returns an object of class \code{sns} with elements:
   \item{samplesMat}{A matrix object with \code{nsample} rows and \code{K} cols.}
   \item{acceptance}{Metropolis proposal percentage acceptance.}
   \item{burn.iters}{Number of burn-in ierations.}
   \item{sample.time}{Time in seconds spent in sampling.}
   \item{burnin.time}{Time in seconds spent in burn-in.}
}

\references{
Mahani A.S., Hasan A., Jiang M. &  Sharabiani M.T.A. (2016). Stochastic Newton Sampler: The R Package sns. Journal of Statistical Software, Code Snippets, 74(2), 1-33. doi:10.18637/jss.v074.c02
}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}
\note{
1. \code{sns.run} cannot be used if SNS is being run as part of a Gibbs cycle, such that the conditional distribution being sampled by SNS changes from one iteration to next. In such cases, \code{sns} must be used instead, inside an explicit Gibbs-cycle \code{for} loop.

2. See package vignette for more details on SNS theory, software, examples, and performance.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{sns}, \code{\link{summary.sns}}, \code{\link{plot.sns}}, \code{\link{predict.sns}}
}

\examples{
\dontrun{

# using RegressionFactory for generating log-likelihood and its derivatives
library(RegressionFactory)

loglike.poisson <- function(beta, X, y) {
  regfac.expand.1par(beta, X = X, y = y,
    fbase1 = fbase1.poisson.log)
}

# simulating data
K <- 5
N <- 1000
X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
beta <- runif(K, -0.5, +0.5)
y <- rpois(N, exp(X %*% beta))

beta.init <- rep(0.0, K)

# glm estimate (ML), for reference
beta.glm <- glm(y ~ X - 1, family = "poisson",
                start = beta.init)$coefficients

# sampling of likelihood
beta.smp <- sns.run(init = beta.init
  , fghEval = loglike.poisson, niter = 1000
  , nnr = 20, X = X, y = y)
smp.summ <- summary(beta.smp)

# compare mean of samples against ML estimate (from glm)
cbind(beta.glm, smp.summ$smp$mean)

# trying numerical differentiation
loglike.poisson.fonly <- function(beta, X, y) {
  regfac.expand.1par(beta, X = X, y = y, fgh = 0,
                     fbase1 = fbase1.poisson.log)
}
beta.smp <- sns.run(init = beta.init
  , fghEval = loglike.poisson.fonly, niter = 1000, nnr = 20
  , X = X, y = y, numderiv = 2)
smp.summ <- summary(beta.smp)
cbind(beta.glm, smp.summ$smp$mean)

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sampling}
\keyword{multivariate}
\keyword{mcmc}
\keyword{Metropolis}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
