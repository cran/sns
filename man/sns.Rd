\name{sns}
\alias{sns}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Stochastic Newton Sampler (SNS)}

\description{
SNS is a Metropolis-Hastings MCMC sampler with a multivariate Gaussian proposal function resulting from a local, second-order Taylor series expansion of log-density. The mean of the Gaussian proposal is identical to the full Newton-Raphson step from the current point. During burn-in, Newton-Raphson optimization can be performed to get close to the mode of the pdf which is unique due to convexity, resulting in faster convergence. For high dimensional densities, state space partitioning can be used to improve mixing. Support for numerical differentiation is provided using \pkg{numDeriv} package. \code{sns} is the low-level function for drawing one sample from the distribution. For drawing multiple samples from a (fixed) distribution, consider using \code{sns.run}.
}
\usage{
sns(x, fghEval, rnd = TRUE, gfit = NULL, mh.diag = FALSE
  , part = NULL, numderiv = 0
  , numderiv.method = c("Richardson", "simple")
  , numderiv.args = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Current state vector.}
  \item{fghEval}{Log-density to be sampled from. A valid log-density can have one of 3 forms: 1) return log-density, but no gradient or Hessian, 2) return a list of \code{f} and \code{g} for log-density and its gradient vector, respectively, 3) return a list of \code{f}, \code{g}, and \code{h} for log-density, gradient vector, and Hessian matrix. Missing derivatives are computed numerically.}
  \item{rnd}{Runs 1 iteration of Newton-Raphson optimization method (non-stochastic or 'nr' mode) when \code{FALSE}. Runs Metropolis-Hastings (stochastic or 'mcmc' mode) for drawing a sample when \code{TRUE}.}
  \item{gfit}{Gaussian fit at point \code{init}. If \code{NULL} then \code{sns} will compute a Gaussian fit at \code{x}.}
  \item{mh.diag}{Boolean flag, indicating whether detailed MH diagnostics such as components of acceptance test must be returned or not.}
  \item{part}{List describing partitioning of state space into subsets. Each element of the list must be an integer vector containing a set of indexes (between \code{1} and \code{length(x)} or \code{length(init)}) indicating which subset of all dimensions to jointly sample. These integer vectors must be mutually exclusive and collectively exhaustive, i.e. cover the entire state space and have no duplicates, in order for the partitioning to represent a valid Gibbs sampling approach. See \code{sns.make.part} and \code{sns.check.part}.}
  \item{numderiv}{Integer with value from the set \code{0,1,2}. If \code{0}, no numerical differentiation is performed, and thus \code{fghEval} is expected to supply \code{f}, \code{g} and \code{h}. If \code{1}, we expect \code{fghEval} to provide \code{f} amd \code{g}, and Hessian will be calculated numerically. If \code{2}, \code{fghEval} only returns log-density, and numerical differentiation is needed to calculate gradient and Hessian.}
  \item{numderiv.method}{Method used for numeric differentiation. This is passed to the \code{grad} and \code{hessian} functions in \pkg{numDeriv} package. See the package documentation for details.}
  \item{numderiv.args}{Arguments to the numeric differentiation method chosen in \code{numderiv.method}, passed to \code{grad} and \code{hessian} functions in \pkg{numDeriv}. See package documentation for details.}
  \item{\dots}{Other arguments to be passed to \code{fghEval}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{sns} returns the sample drawn as a vector, with attributes:
  \item{accept}{A boolean indicating whether the proposed point was accepted.}
  \item{ll}{Value of the log-density at the sampled point.}
  \item{gfit}{List containing Gaussian fit to pdf at the sampled point.}
}

\references{
Mahani A.S., Hasan A., Jiang M. &  Sharabiani M.T.A. (2016). Stochastic Newton Sampler: The R Package sns. Journal of Statistical Software, Code Snippets, 74(2), 1-33. doi:10.18637/jss.v074.c02

Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains and their applications. Biometrika, 57(1), 97-109.

Qi, Y., & Minka, T. P. (2002). Hessian-based markov chain monte-carlo algorithms. 1st Cape Cod Workshop on Monte Carlo Methods.

}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}

\note{
1. Since SNS makes local Gaussian approximations to the density with the covariance matrix of the Gaussian proposal being the log-density Hessian, there is a strict requirement for the log-density to be concave.  

2. Proving log-concavity for arbitrary probability distributions is non-trvial. However, distributions \emph{generated} by replacing parameters of a concave distribution with linear expressions are known to be log-concave. This negative-definiteness invariance as well as expressions for full gradient and Hessian in terms of derivatives of low-dimensional base distributions are discussed in the vignette. The GLM expansion framework is available in the R package \pkg{RegressionFactory}.

3. See package vignette for more details on SNS theory, software, examples, and performance.
 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sns.run}}, \code{\link{sns.fghEval.numaug}}
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
y <- rpois(N, exp(X \%*\% beta))

beta.init <- rep(0.0, K)

# glm estimate, for reference
beta.glm <- glm(y ~ X - 1, family = "poisson",
                start = beta.init)$coefficients

# running SNS in non-stochastic mode
# this should produce results very close to glm
beta.sns <- beta.init
for (i in 1:20)
  beta.sns <- sns(beta.sns, loglike.poisson, X = X, y = y, rnd = F)

# comparison
all.equal(as.numeric(beta.glm), as.numeric(beta.sns))

# trying numerical differentiation
loglike.poisson.fonly <- function(beta, X, y) {
  regfac.expand.1par(beta, X = X, y = y, fgh = 0,
                     fbase1 = fbase1.poisson.log)
}

beta.sns.numderiv <- beta.init
for (i in 1:20)
  beta.sns.numderiv <- sns(beta.sns.numderiv, loglike.poisson.fonly
                  , X = X, y = y, rnd = F, numderiv = 2)
all.equal(as.numeric(beta.glm), as.numeric(beta.sns.numderiv))

# add numerical derivatives to fghEval outside sns
loglike.poisson.numaug <- sns.fghEval.numaug(loglike.poisson.fonly
  , numderiv = 2)

beta.sns.numaug <- beta.init
for (i in 1:20)
  # set numderiv to 0 to avoid repeating 
  # numerical augmentation inside sns
  beta.sns.numaug <- sns(beta.sns.numaug, loglike.poisson.numaug
                           , X = X, y = y, rnd = F, numderiv = 0)
all.equal(as.numeric(beta.glm), as.numeric(beta.sns.numaug))

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sampling, multivariate, mcmc, Metropolis}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
