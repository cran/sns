\name{sns}
\alias{sns}
\alias{sns.run}
\alias{print.sns}

\title{Stochastic Newton Sampler (SNS)}

\description{
SNS is a Metropolis-Hastings MCMC sampler with a multivariate Gaussian proposal function resulting from a local, second-order Taylor series expansion of log-density. The mean of the Gaussian proposal is identical to the full Newton-Raphson step from the current point. During burn-in, Newton-Raphson optimization can be performed to get close to the mode of the pdf which is unique due to convexity, resulting in faster convergence. For high dimensional densities, state space partitioning can be used to improve mixing.
}

\usage{
sns(x, fghEval, rnd=TRUE, gfit=NULL, mh.diag = FALSE, part = NULL, ...)
sns.run(init, fghEval, niter = 100, nnr = min(10, round(niter/4))
  , mh.diag = FALSE, part = NULL
  , print.level = 0, report.progress = ceiling(niter/10), ...)
\method{print}{sns}(x, ...)
}

\arguments{
    \item{x}{For \code{sns}, the current state vector. For \code{print.sns}, an object of class \code{sns}.}
    \item{fghEval}{function, should evaluate the log-density, its gradient and Hessian at any point. Must a return a list with labels: f - the log-probability density, g - gradient vector, h - Hessian matrix.}
    \item{rnd}{Runs 1 iteration of Newton-Raphson optimization method (non-stochastic or 'nr' mode) when \code{FALSE}. Runs Metropolis-Hastings (stochastic or 'mcmc' mode) for drawing a sample when \code{TRUE}.}
    \item{gfit}{Gaussian fit at point \code{init}. If \code{NULL} then \code{sns} will compute a Gaussian fit at \code{x}.}
    \item{mh.diag}{Boolean flag, indicating whether detailed MH diagnostics such as components of acceptance test must be returned or not.}
    \item{part}{List describing partitioning of state space into subsets. Each element of the list must be an integer vector containing a set of indexes (between \code{1} and \code{length(x)} or \code{length(init)}) indicating which subset of all dimensions to jointly sample. These integer vectors must be mutually exclusive and collectively exhaustive, i.e. cover the entire state space and have no duplicates, in order for the partitioning to represent a valid Gibbs sampling approach.}
    \item{init}{Starting state for \code{sns.run}.}
    \item{niter}{Number of iterations to perform (in 'nr' and 'mcmc' mode).}
    \item{nnr}{Number of Newton-Raphson (non-stochastic) iterations to perform at the beginning.}
    \item{print.level}{If greater than 0, print sampling progress report.}
    \item{report.progress}{Number of sampling iterations to wait before printing progress reports.}
    \item{...}{Arguments passed to \code{fghEval} in the case of \code{sns} and \code{sns.run}.}
}

\value{
   \code{sns.run} returns an object of class \code{sns} with elements:
   \item{samplesMat}{A matrix object with \code{nsample} rows and \code{K} cols.}
   \item{acceptance}{Metropolis proposal percentage acceptance.}
   \item{burn.iters}{Number of burn-in ierations.}
   \item{sample.time}{Time in seconds spent in sampling.}
   \item{burnin.time}{Time in seconds spent in burn-in.}

   \code{sns} returns the sample drawn as a vector, with attributes:
   \item{accept}{A boolean indicating whether the proposal was accepted.}
   \item{ll}{Value of the log-pdf at the sampled point.}
   \item{gfit}{List containing Gaussian fit to pdf at the sampled point.}
}

\note{
1. \code{sns.run} should be the function called by users when drawing multiple samples from the same density, while \code{sns} is a lower level function which can be directly called if required, for example when the sampled density is a conditional density and part of Gibbs sampling cycle, in which the conditional density is changing from one iteration to the next.

2. Proving log-concavity for arbitrary probability distributions is non-trvial. However distributions \emph{generated} by replacing parameters of a concave distribution with linear expressions are known to be log-concave. \emph{Mahani and Sharabiani (2013)} show how to automatically generate Hessian and gradient for the expanded parameter set distribution given those of the base distribution. The GLM expansion framework is available in the R package RegressionFactory.
 
3. Since SNS makes local Gaussian approximations to the density with the covariance being the Hessian, there is a strict requirement for the pdf to be log-concave.  
}

\seealso{
\code{\link{ess}}, \code{\link{summary.sns}}, \code{\link{plot.sns}}, \code{\link{predict.sns}}
}

\author{Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani}
\keyword{sampling, multivariate, mcmc, Metropolis}

\references{
 Mahani, Alireza S. and Sharabiani, Mansour T.A. (2013)
 \emph{Metropolis-Hastings Sampling Using Multivariate Gaussian Tangents}
 \url{http://arxiv.org/pdf/1308.0657v1.pdf} 
 
 Qi, Y. and Minka, T.P. (2002)
 \emph{Hessian-Based Markov Chain Monte Carlo Algorithms}
}

\examples{
\dontrun{

# use this library for constructing the
# log-likelihood and its gradient and Hessian
library(RegressionFactory)

# simulate data for logistic regression
N <- 1000 # number of observations
K <- 5 # number of attributes
X <- matrix(runif(N*K, min = -0.5, max = +0.5), ncol = K)
beta <- runif(K, min = -0.5, max = +0.5)
y <- 1*(runif(N) < 1/(1 + exp(- X \%*\% beta)))

# define log-likelihood using expander framework of
# RegressionFactory package
loglike <- function(beta, X, y) {
  regfac.expand.1par(beta, X, y, fbase1.binomial.logit, fgh = 2, n = 1)
}

# SNS sampling
nsmp <- 100
beta.init <- rep(0.0, K)
sns.out <- sns.run(beta.init, loglike, 200, X = X, y = y)
summary(sns.out)

}
}

