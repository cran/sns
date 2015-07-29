\name{sns.fghEval.numaug}
\alias{sns.fghEval.numaug}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility function for augmentation of a log-density function with numerical gradient and Hessian as needed
}
\description{
Augmenting a log-density with numerical gradient and Hessian, so it can be used by \code{sns} or \code{sns.run}. This augmentation will also be done inside the function, if the value of \code{numderiv} parameter passed to \code{sns} and \code{sns.run} is \code{1} or \code{2}. The advantage of using \code{sns.fghEval.numaug} outside these functions is efficiency, since the agumentation code will not have to be executed in every function call. Users must set \code{numderiv} to \code{0} when calling \code{sns} or \code{sns.run} if calling \code{sns.fghEval.numaug} first. See example.
}
\usage{
sns.fghEval.numaug(fghEval, numderiv = 0
  , numderiv.method = c("Richardson", "simple")
  , numderiv.args = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fghEval}{Log-density to be sampled from. A valid log-density can have one of 3 forms: 1) return log-density, but no gradient or Hessian, 2) return a list of \code{f} and \code{g} for log-density and its gradient vector, respectively, 3) return a list of \code{f}, \code{g}, and \code{h} for log-density, gradient vector, and Hessian matrix. Missing derivatives are computed numerically.}
  \item{numderiv}{This must be matched with \code{fghEval}: Integer with value from the set \code{0,1,2}. If \code{0}, no numerical differentiation is performed, and thus \code{fghEval} is expected to supply \code{f}, \code{g} and \code{h}. If \code{1}, we expect \code{fghEval} to provide \code{f} amd \code{g}, and Hessian will be calculated numerically. If \code{2}, \code{fghEval} only returns log-density, and numerical differentiation is needed to calculate gradient and Hessian.}
  \item{numderiv.method}{Method used for numeric differentiation. This is passed to the \code{grad} and \code{hessian} functions in \pkg{numDeriv} package. See the package documentation for details.}
  \item{numderiv.args}{Arguments to the numeric differentiation method chosen in \code{numderiv.method}, passed to \code{grad} and \code{hessian} functions in \pkg{numDeriv}. See package documentation for details.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A function, accepting same arguments as \code{fghEval}, but guaranteed to return the original log-density, plus gradient and Hessian (both of which could possibly by numerically calculated). If \code{numderiv=0}, \code{fghEval} is returned without change. The function will return log-density, gradient and Hessian as elements \code{f}, \code{g} and \code{h} of a list.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}

\note{
See package vignette for more details on SNS theory, software, examples, and performance.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sns}}, \code{\link{sns.run}}
}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
