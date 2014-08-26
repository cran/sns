\name{bayesGLM}
\alias{glmfgh}

\title{Bayesian Regression of Generalized Linear Models}

\description{
  A closure which generates a function, gradient, Hessian evaluator using for log-likelihood functions of 4 types of GLMs. 
}

\usage{
glmfgh(N, K, glmtype = "logistic", X=NULL, y=NULL) 
}

\arguments{
    \item{N}{Number of observations.}
    \item{K}{Number of variables.}
    \item{glmtype}{Must be one of \code{logistic}, \code{poisson}, \code{geometric}, \code{exponential}. }
    \item{X}{Data matrix of the explanatory variable, must have \code{N} rows and \code{K} cols. If \code{NULL} then simulated data is used.}
    \item{y}{Vector having the dependent variable (response) corresponging to \code{X}}
}

\value{
  A function suitable for use in \code{sns.run} as a function, gradient, Hessian evaluator.
}

\note{
GLMs are distributions where a linear expression replaces the parameter in some `base' distributions. A general framework to compute the expanded Hessian from that of the base, following the `Log-concavity invariance theorem' (\code{Appendix A} of the preprint \emph{Mahani and Sharabiani (2013)}), is used to generate the Hessian and gradient.

\code{expand.1par.R} has the Hessian and gradient expander code for probability distributions derived from any one paramter base distributions. The derived distributions are guaranteed to be log-concave if the base distribution is log-concave as shown in the referred paper.  
}

\references{
 Mahani, Alireza S. and Sharabiani, Mansour T.A. (2013)
 \emph{Metropolis-Hastings Sampling Using Multivariate Gaussian Tangents}
 \url{http://arxiv.org/pdf/1308.0657v1.pdf}
}


\examples{
 \dontrun{
  # Evaluator for Logistic log-likelihood, gradient, Hessian
  fghEval <- glmfgh(N=100, K=2, glmtype="logistic")
 }
}
