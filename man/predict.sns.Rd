\name{predict.sns}
\alias{predict.sns}
\alias{summary.predict.sns}
\alias{print.summary.predict.sns}

\title{
Sample-based prediction using "sns" Objects
}

\description{
Method for sample-based prediction using the output of \code{\link{sns.run}}.
}

\usage{
\method{predict}{sns}(object, fpred
  , nburnin = max(nrow(object)/2, attr(object, "nnr"))
  , end = nrow(object), thin = 1, ...)
\method{summary}{predict.sns}(object
  , quantiles = c(0.025, 0.5, 0.975)
  , ess.method = c("coda", "ise"), ...)
\method{print}{summary.predict.sns}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class "sns" (output of \code{\link{sns.run}}) or "predict.sns" (output of \code{predict.sns}).}
  \item{fpred}{Prediction function, accepting a single value for the state vector and producing a vector of outputs.}
  \item{nburnin}{Number of burn-in iterations discarded for sample-based prediction.}
  \item{end}{Last iteration used in sample-based prediction.}
  \item{thin}{One out of \code{thin} iterations within the specified range are used for sample-based prediction.}
  \item{quantiles}{Values for which sample-based quantiles are calculated.}
  \item{ess.method}{Method used for calculating effective sample size. Default is to call \code{effectiveSize} from package \code{coda}.}
  \item{x}{An object of class "summary.predict.sns".}
  \item{...}{Arguments passed to/from other functions.}
}

\value{
\code{predict.sns} produces a matrix with number of rows equal to the length of prediction vector produces by \code{fpred}. Its numnber of columns is equal to the number of samples used within the user-specified range, and after thinning (if any). \code{summary.predict.sns} produces sample-based prediction mean, standard deviation, quantiles, and effective sample size.
}

\references{
Mahani A.S., Hasan A., Jiang M. &  Sharabiani M.T.A. (2016). Stochastic Newton Sampler: The R Package sns. Journal of Statistical Software, Code Snippets, 74(2), 1-33. doi:10.18637/jss.v074.c02
}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}

\note{
See package vignette for more details on SNS theory, software, examples, and performance.
}

\seealso{
\code{\link{sns.run}}
}

\examples{

\dontrun{

# using RegressionFactory for generating log-likelihood and derivatives
library("RegressionFactory")

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
beta.smp <- sns.run(beta.init, loglike.poisson,
  niter = 1000, nnr = 20, mh.diag = TRUE, X = X, y = y)

# prediction function for mean response
predmean.poisson <- function(beta, Xnew) exp(Xnew \%*\% beta)
ymean.new <- predict(beta.smp, predmean.poisson,
                     nburnin = 100, Xnew = X)
summary(ymean.new)

# (stochastic) prediction function for response
predsmp.poisson <- function(beta, Xnew)
  rpois(nrow(Xnew), exp(Xnew \%*\% beta))
ysmp.new <- predict(beta.smp, predsmp.poisson
                    , nburnin = 100, Xnew = X)
summary(ysmp.new)

}

}


