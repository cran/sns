\name{sns.check.logdensity}
\alias{sns.check.logdensity}
\alias{print.sns.check.logdensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility function for validating log-density
}
\description{
Utility function for validating log-density: 1) dimensional consistency of function argument, gradient and Hessian, 2) finiteness of function, gradient and Hessian, 3) closeness of analytical and numerical derivatives, and 4) negative definiteness of Hessian.
}
\usage{
sns.check.logdensity(x, fghEval
  , numderiv.method = c("Richardson", "complex")
  , numderiv.args = list()
  , blocks = append(list(1:length(x)), as.list(1:length(x)))
  , dx = rep(1, length(x)), nevals = 100, negdef.tol = 1e-08, ...)
\method{print}{sns.check.logdensity}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{For \code{sns.check.logdensity}, initial point, around which a random collection of points are generated to perform validation tests. For \code{print.sns.check.logdensity}, an object of class \code{sns.check.logdensity}, typically the output of \code{sns.check.logdensity} function.}
  \item{fghEval}{Log-density to be validated. A valid log-density can have one of 3 forms: 1) return log-density, but no gradient or Hessian, 2) return a list of \code{f} and \code{g} for log-density and its gradient vector, respectively, 3) return a list of \code{f}, \code{g}, and \code{h} for log-density, gradient vector, and Hessian matrix.}
  \item{numderiv.method}{Method used for numeric differentiation. This is passed to the \code{grad} and \code{hessian} functions in \pkg{numDeriv} package. See the package documentation for details.}
  \item{numderiv.args}{Arguments to the numeric differentiation method chosen in \code{numderiv.method}, passed to \code{grad} and \code{hessian} functions in \pkg{numDeriv}. See package documentation for details.}
  \item{blocks}{A list of state space subsets (identified by their positional indexes), for which negative-definiteness of Hessian blocks are to be tested. The default is to test for 1) entire state space, and 2) each dimension individually.}
  \item{dx}{A vector of same length as \code{x}. For \code{i}'th dimension, \code{nevals} values are sampled from a uniform distribution with min/max values equal to \code{x[i]-0.5*dx[i]} and \code{x[i]+0.5*dx[i]}, respectively. Vectors smaller than \code{length(x)} are extended as needed by recycling the provided values for \code{dx}.}
  \item{nevals}{Number of points in state space, for which validation tests will be performed.}
  \item{negdef.tol}{Lower bound for absolute value of (negative) eigenvalues of Hessian, evaluated at each of the \code{nevals} points in the state space. If one or more eigenvalues have absolute values smaller than \code{ngdef.tol}, log-density is declared non-log-concave at that point.}
  \item{\dots}{Other arguments to be passed to \code{fghEval}.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\code{sns.check.logdensity} returns a list of class \code{sns.check.logdensity}, with the following elements:
  \item{check.ld.struct}{Boolean flag, indicating whether log-density \code{fghEval} has one of the 3 forms of output, described above.}
  \item{numderiv}{Integer with values of \code{0,1,2}. A value of \code{0} means analytical gradient and Hessian have been provided, and thus there is no need for numerical differentiation. \code{1} means analytical gradient is provided, but Hessian must be calculated numerically. \code{2} means both gradient and Hessian must be numerically calculated. Users can pass this value to subsequent \code{sns} or \code{sns.run} calls.}
  \item{check.length.g}{Boolean flag, indicating whether length of gradient vector (element \code{g}) returned by \code{fghEval} equals \code{length(x)}.}
  \item{check.dim.h}{Boolean flag, indicating whether number of rows and columns of the Hessian matrix (element \code{h}) returned by \code{fghEval} equal \code{length(x)}.}
  \item{x.mat}{Collection of state space vectors (one per row), for which validation tests are performed. It has \code{nevals} rows and \code{length(x)} columns.}
  \item{t.evals}{Time spent on evaluating \code{fghEval} on \code{nevals} points chosen randomly in the neighborhood of \code{x}, as specified by \code{dx}. This includes log-density and, if provided, analytical evaluations of gradient and Hessian.}
  \item{t.num.evals}{Time spent on evaluating the numeric version of \code{fghEval}, in which gradient and Hessian are computed numerically, using \code{grad} and \code{hessian} functions in the \pkg{numDeriv} package. Comparison of this number with \code{t.evals} provides the user with insight into the relative speed of numerical differentiation compared to analytical versions.}
  \item{f.vec}{Vector of log-density values for state space vectors listed in \code{x.mat}.}
  \item{g.mat.num}{Collection of numerically-computed gradient vectors for state space values listed in \code{x.mat}, with the same dimension conventions.}
  \item{is.g.num.finite}{Boolean flag, indicating whether all numerically-computed gradient vectors have finite values.}
  \item{h.array.num}{Collection of numerically-computed Hessian matrices at points listed in \code{x.mat}. First dimension is of length \code{nevals}, and the remaining two dimensions equal \code{length(x)}.}
  \item{is.h.num.finite}{Boolean flag, indicating whether all numerically-computed Hessian matrices have finite values.}
  \item{g.mat}{Collection of analytically-computed gradient vectors for state space values listed in \code{x.mat}, with the same dimension conventions. This is only available if \code{fghEval} has a \code{g} field; otherwise \code{NA}.}
  \item{is.g.finite}{Boolean flag (if available), indicating whether all analytically-computed gradient vectors have finite values (if available).}
  \item{g.diff.max}{If available, maximum relative difference between analytical and numerical gradient vectors, over all \code{nevals} points in \code{x.mat}. Relative diference is defined as L2 norm of difference between the two gradient vectors, divided by the L2 norm of the analytical gradient vector.}
  \item{h.array}{If available, collection of analytically-computed Hessian matrices at points listed in \code{x.mat}. Dimensional conventions are the same as \code{h.array.num}.}
  \item{is.h.finite}{Boolean flag (if available), indicating whether all analytically-computed Hessian matrices have finite values.}
  \item{h.diff.max}{If available, maximum relative difference between analytical and numerical Hessian matrices, over all \code{nevals} points in \code{x.mat}. Relative difference is defined as the Frobenius norm of difference of analytical and numerical Hessian matrices, divided by the Frobenius norm of analytical Hessian.}
  \item{is.negdef.num}{Boolean flag, indicating whether numerical Hessian is negative-definite at all state space points indicated in \code{x.mat}.}
  \item{is.negdef}{Boolean flag, indicating whether analytical Hessian is negative-definite at all state space points indicated in \code{x.mat}.}
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}
\note{
1. Validation tests performed in \code{sns.check.logdensity} cannot prove that a log-density is twice-differentiable, or globally concave. However, when e.g. log-density Hessian is seen to be non-negative-definite at one of the points tested, we can definitively say that the Hessian is not globally negative-definite, and therefore \code{sns} should not be used for sampling from this distribution. Users must generally consider this function as a supplement to analytical work.

2. See package vignette for more details on SNS theory, software, examples, and performance.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

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

my.check <- sns.check.logdensity(beta.init, loglike.poisson
  , X = X, y = y, blocks = list(1:K))
my.check

# mistake in log-likelihood gradient
loglike.poisson.wrong <- function(beta, X, y) {
  ret <- regfac.expand.1par(beta, X = X, y = y,
                            fbase1 = fbase1.poisson.log)
  ret$g <- 1.2 * ret$g
  return (ret)
}
# maximum relative diff in gradient is now much larger
my.check.wrong <- sns.check.logdensity(beta.init
  , loglike.poisson.wrong, X = X, y = y, blocks = list(1:K))
my.check.wrong

# mistake in log-likelihood Hessian
loglike.poisson.wrong.2 <- function(beta, X, y) {
  ret <- regfac.expand.1par(beta, X = X, y = y,
                            fbase1 = fbase1.poisson.log)
  ret$h <- 1.2 * ret$h
  return (ret)
}
# maximum relative diff in Hessian is now much larger
my.check.wrong.2 <- sns.check.logdensity(beta.init
  , loglike.poisson.wrong.2, X = X, y = y, blocks = list(1:K))
my.check.wrong.2

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
