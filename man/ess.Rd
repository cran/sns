\name{ess}
\alias{ess}

\title{Effective Sample Size Calculator}

\description{
Computes the effective sample size of MCMC chains, using the algorithm in Section 2.3 of the paper by Madeline Thompson. The algorithm is taken from earlier work on `Initial Sequence Estimators' by multiple authors.
}

\usage{
ess(x, method = c("coda", "ise"))
}

\arguments{
    \item{x}{Matrix object with each sample (possibly multivariate) as a row. Effective sample size calculation is done independently for each column of \code{x}.}
    \item{method}{Method of calculating effective size. Current options are "coda" which calls \code{effectiveSize} function in \code{coda} package, and "ise" which uses the 'Initial Sequence Estimators' method described in Section 2.3 of Thompson (2010).}
}

\value{
  Vector with effective sample sizes for the time series in each column of \code{x}.
}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}

\references{
  Thompson, Madeleine (2010) 
  \emph{A Comparison of Methods for Computing Autocorrelation Time}
  \url{http://arxiv.org/pdf/1011.0175v1.pdf}
}

