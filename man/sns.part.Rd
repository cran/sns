\name{sns.make.part}
\alias{sns.make.part}
\alias{sns.check.part}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility Functions for Creating and Validating State Space Partitions
}
\description{
Utility functions for creating and validating state space partitions, to be used in SNS for improving the mixing of sampled chains for high-dimensional posteriors.
}
\usage{
sns.make.part(K, nsubset, method = "naive")
sns.check.part(part, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K}{Dimensionality of state space.}
  \item{nsubset}{Number of subsets to partition the state space dimensions into.}
  \item{method}{Method used for state space partitioning. Currently, only \code{naive} method is implemented, where coordinates are distributed evenly (or as evenly as possible) across subsets.}
  \item{part}{A list of length \code{nsubset}, with each element a vector of integer values, representing the coordinates belonging to a subset. This list is the output of \code{sns.make.part}.}
}

\value{
\code{sns.make.part} produces a list of integer vectors, each containing coordinates belonging to the same subset. \code{sns.check.part} produces a boolean flag, indicating whether or not the partition list is valid or not. The subset members must constitute a mutually-exclusive, collectively-exhaustive set relative to \code{1:K}.
}

\references{
Mahani A.S., Hasan A., Jiang M. &  Sharabiani M.T.A. (2016). Stochastic Newton Sampler: The R Package sns. Journal of Statistical Software, Code Snippets, 74(2), 1-33. doi:10.18637/jss.v074.c02
}

\author{
Alireza S. Mahani, Asad Hasan, Marshall Jiang, Mansour T.A. Sharabiani
}

\seealso{
\code{\link{sns}}, \code{\link{sns.run}}
}

\examples{
# creating a valid partition of a 6-dimensional state space
my.part.valid <- list(c(1,2,3), c(4,5,6))
is.valid.1 <- sns.check.part(my.part.valid, 6)
cat("is partition valid: ", is.valid.1, "\n")

# creating an invalid partition of a 6-dimensional state space
# (coordinate 4 is missing)
my.part.invalid <- list(c(1,2,3), c(5,6))
is.valid.2 <- sns.check.part(my.part.invalid, 6)
cat("is partition valid: ", is.valid.2, "\n")

}
