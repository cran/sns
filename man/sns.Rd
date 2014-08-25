\name{sns}
\alias{sns}
\alias{sns.run}
\alias{print.sns}
\alias{summary.sns}

\title{Stochastic Newton Sampler (SNS)}

\description{
Multivariate, Markov Chain Monte Carlo sampler for twice differentiable, log-concave probability densities. SNS is a Metropolis-Hastings variant, with a special form of the proposal function. The proposal function is a local Gaussian approximation to the density obtained via Taylor expansion of the log-density around the current point. The mean of the Gaussian proposal is taken as the full Newton-Raphson step from the current point. During burn-in, Newton-Raphson optimization is performed to get close to the mode of the pdf which is unique due to convexity.
}

\usage{
sns(init, fghEval, rnd=TRUE, gfit=NULL, ...)
sns.run(K, nburnin, nsample, fghEval, start=NULL, print.level=0, report.progress=100, ...)
\method{print}{sns}(x, ...)
\method{summary}{sns}(object, show.means = FALSE, ...)
}

\arguments{
    \item{init}{Starting point for the Metropolis-Hastings procedure.}
    \item{fghEval}{function, should evaluate the log-density, its gradient and Hessian at any point. Must a return a list with labels: f - the log-probability density, g - gradient vector, h - Hessian matrix.}
    \item{rnd}{Runs 1 iteration of Newton's method (non-stochastic) when \code{FALSE}. Runs Metropolis-Hastings for draw a sample when \code{TRUE}. NOTE: Set to \code{FALSE} only during burn-in.}
    \item{gfit}{Gaussian fit at point \code{init}. If \code{NULL} then \code{sns} will compute a Gaussian fit at \code{init}.}
    \item{K}{Dimension of the space to draw samples from.}
    \item{nburnin}{Number of burn-in iterations (non-stochastic, Newton-Raphson).}
    \item{nsample}{Number of samples to draw (after burn-in).}
    \item{start}{Initial point for the Markov chain. Default: \code{rep(0.1, K)}.}
    \item{print.level}{If greater than 0, print sampling progress report.}
    \item{report.progress}{Number of sampling iterations between progress reports.}
    \item{...}{Extra args all passed to \code{fghEval} whenever it's called.}
    \item{x, object}{A \code{sns} class object.}
    \item{show.means}{Setting to \code{TRUE} enables summary.sns to print sample means and standard deviations for all variables.}
}

\value{
   \code{sns.run} returns an object of class \code{sns} with elements:
   \item{samplesMat}{A matrix object with \code{nsamples} rows and \code{K} cols.}
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
1. \code{sns.run} should be function called by users, while \code{sns} is a lower level function which can be directly called if required.

2. Proving log-concavity for arbitrary probability distributions is non-trvial. However distributions \emph{generated} by replacing parameters of a concave distribution with linear expressions are known to be log-concave. \emph{Mahani and Sharabiani (2013)} show how to automatically generate Hessian and gradient for the expanded parameter set distribution given those of the base distribution. For one paramter distributions, like GLMs, this expansion is used to implement function, gradient, Hessian evaluators for Bayesian Generalized Linear regression models in \code{\link{glmfgh}}.
 
3. Since SNS makes local Gaussian approximations to the log-pdf with the covariance being the Hessian, there is a strict requirement for the pdf to be log-concave.  
}

\seealso{\code{\link{ess}}, \code{\link{glmfgh}}}
\author{Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.S. Sharabiani}
\keyword{sampling, multivariate, mcmc, Metropolis}

\references{
 Mahani, Alireza S. and Sharabiani, Mansour T.S. (2013)
 \emph{Metropolis-Hastings Sampling Using Multivariate Gaussian Tangents}
 \url{http://arxiv.org/pdf/1308.0657v1.pdf} 
 
 Qi, Y. and Minka, T.P. (2002)
 \emph{Hessian-based markov Chain Monte Carlo algorithms}
}

\examples{
   library(sns)
   # Logistic log-likelihood, gradient, Hessian
   fghEval <- glmfgh(N=1000, K=5, glmtype="logistic")
   # Draw samples using SNS 
   sns.out <- sns.run(K=5, nburnin=20, nsample=1000, fghEval)
   summary(sns.out, show.means=TRUE)
}
