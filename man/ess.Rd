\name{ess}
\alias{ess}

\title{Effective Sample Size Calculator}

\description{
Computes the effective sample size using the algorithm in Section 2.3 of the paper by Madeline Thompson. The algorithm is taken from earlier work on 'Initial Sequence Estimators' by multiple authors. The 'magic' argument is inspired by the Appendix in the No U-Turn Sampler paper by Hoffman & Gelman.
}

\usage{
ess(chain, mu=NULL, adj=TRUE, magic=0.05)
}

\arguments{
    \item{chain}{Matrix object with each sample (possibly multivariate) as a row.}
    \item{mu}{Vector of means with length = \code{ncol(chain)}.}
    \item{adj}{Set to \code{TRUE} to enable Initial Convex Sequence Estimator (see Section
2.3 of the paper by Thompson).}
    \item{magic}{cutoff used in Appendix of the NUTS paper by Gelman & Hoffmann. NOT used if \code{adj == TRUE}.}
}

\value{
  effective sample sizes for the time series in each column of \code{chain}.
}

\note{
Using the 'magic' argument (and choosing its value) is not mathematically rigorous in all situations. However, it's widely used in practice (e.g. in work of Hoffman and Gelman) and usually produces sensible results.
}

\references{
  Thompson, Madeline (2010) 
  \emph{A Comparison of Methods for Computing Autocorrelation Time}
  \url{http://arxiv.org/pdf/1011.0175v1.pdf}

  Hoffman, Matthew D. and  Gelman, Andrew (2014)
  \emph{The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo}
  \url{http://jmlr.org/papers/v15/hoffman14a.html} 
}

\examples{
  \dontrun{
   # Evaluator for Logistic log-likelihood, gradient, Hessian
   fghEval <- glmfgh(N=100, K=2, glmtype="logistic")
   # Draw samples using SNS
   sns.out <- runSNS(K, burnin.iters=20, nsamples=1000, fghEval)
   # Compute effective sample size
   apply(sns.out$samplesMat, 2, ess)
  }
}
