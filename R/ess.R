############################################################################# 
# Computes the effective sample size using the algorithm in Section 2.3 of
# the paper (http://arxiv.org/pdf/1011.0175v1.pdf) by Madeline Thompson.
# The algorithm is taken from earlier work on 'Initial Sequence Estimators'
# by multiple authors. 
# 
# The 'magic' argument is isnpired by the Appendix in the NUTS paper by
# Gelman & Hoffmann (2013).
# 
# Args: 
#   chain - matrix object with each sample (possibly multivariate) as a row
#   mu - vector of means (must of the same length as ncols(chain)
#   adj - Set to TRUE to enable Initial Convex Sequence Estimator (see Section
#         2.3 of the paper by Thompson, referred above).
#   magic - cutoff used in Appendix of NUTS paper by Gelman & Hoffmann (2013)
# Returns:
#   effective sample sizes for the time series in each col of 'chain'
############################################################################# 
ess <- function(chain, mu=NULL, adj=TRUE, magic=0.05) 
{
  dims <- dim(chain);
  M <- dims[1];
  K <- dims[2];

  # If M is of 1 dimensions
  if (is.null(M)) {
    M <- length(chain);
  }
  if (is.null(K) || is.na(K)) {
    K <- 1;
  }

  if (is.null(mu)) {
    if (K != 1) {
      mu = colMeans(chain); 
    }
    else {
      mu = mean(chain);
    }
  }

  x <- array();
  for (i in 1:K) {
    if(K != 1) { 
      sigma <- sqrt(rhof(chain[,i], M, 0, mu[i], 1)); # Center first term at 1
    }
    else {
      sigma <- sqrt(rhof(chain, M, 0, mu[i], 1)); # Center first term at 1
    }

    sum <- 0;
    prev <- 0; # storage for previous
    for (s in 1:(M-1)) {
      if (K != 1) {
        rho <- rhof(chain[,i], M, s, mu[i], sigma);
      }
      else {
        rho <- rhof(chain, M, s, mu[i], sigma);
      }

      # Break if less than magic number, or adjacent is negative
      if (!adj && rho < magic) {
        break;
      }
      else if (adj && (prev + rho) <= 0) {
        break;
      }

      else {
        #sum <- sum + (1 - s/M)*rho; # NUTS paper
        sum <- sum + rho; # Thompson version
        prev <- rho;
      }
    }
    x[i] <- (M/(1+2*sum));
  }

  return(x);
}

rhof <- function(chain, M, s, mu, sigma) {
  # Can optimize by saving the chain shifted by mu instead of doing it everytime...
  a <- sum((chain[1:(M-s)] - mu)*(chain[(s+1):M ] - mu)/(sigma^2*(M-s)));
  return(a);
}
