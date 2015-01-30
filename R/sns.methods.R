######################################
#                                    #
#  Methods for sns objects           #
#                                    #
######################################

# convenience function for partitioning the state space
sns.make.part <- function(K, nsubset, method = "naive") {
  if (method != "naive") stop("invalid method") #TODO: consider implementing more sophisticated partitioning methods
  if (nsubset > K) stop("number of partitions cannot exceed state space dimensionality")
  
  # deterimining number of coordinates per subset (nvec)
  nvec <- rep(0, nsubset)
  nleft <- K
  c <- 1
  while (nleft > 0) {
    nvec[c] <- nvec[c] + 1
    nleft <- nleft - 1
    c <- c %% nsubset + 1
  }

  # assigning coordinates to subsets
  ret <- list()
  c <- 0
  for (n in 1:nsubset) {
    ret[[n]] <- as.integer(c + 1:nvec[n])
    c <- c + nvec[n]
  }
  if (sns.check.part(ret, K)) return (ret)
  else stop("unexpectedly invalid state space partitioning")
}

# function for checking that state space partitioning is valid
# (mutually-exclusive and collectively-exhaustive)
sns.check.part <- function(part, K) {
  return (identical(as.integer(sort(unlist(part))), 1:K))
}

# predict methods
predict.sns <- function(object, fpred
  , nburnin = max(nrow(object)/2, attr(object, "nnr"))
  , end = nrow(object), thin = 1, ...) {

  niter <- nrow(object)
  nnr <- attr(object, "nnr")
  nmcmc <- niter - nnr
  if (nburnin < nnr) warning("it is strongly suggested that burnin period includes NR iterations (which are not valid MCMC iterations)")
  myseq <- seq(from = nburnin + 1, to = end, by = thin)

  pred <- apply(object[myseq, ], 1, fpred, ...)
  class(pred) <- "predict.sns"
  return (pred)
}
summary.predict.sns <- function(object, quantiles = c(0.025, 0.5, 0.975)
  , ess.method = c("coda", "ise"), ...) {
  smp.mean <- rowMeans(object)
  smp.sd <- apply(object, 1, sd)
  smp.ess <- ess(t(object), method = ess.method[1])
  smp.quantiles <- t(apply(object, 1, quantile, probs = quantiles))
  ret <- list(mean = smp.mean, sd = smp.sd, ess = smp.ess, quantiles = smp.quantiles, nseq = ncol(object))
  class(ret) <- "summary.predict.sns"
  return (ret)
}
print.summary.predict.sns <- function(x, ...) {
  cat("prediction sample statistics:\n")
  cat("\t(nominal sample size: ", x$nseq, ")\n", sep="")
  stats <- cbind(x$mean, x$sd, x$ess, x$quantiles)
  colnames(stats)[1:3] <- c("mean", "sd", "ess")
  rownames(stats) <- c(1:length(x$mean))
  printCoefmat(stats[1:min(length(x$mean), 6), ])
  if (length(x$mean) > 6) cat("...\n")
}

# print method 
print.sns <- function(x, ...) {
  cat("Stochastic Newton Sampler (SNS)\n")
  cat("state space dimensionality: ", ncol(x), "\n")
  if (!is.null(attr(x, "part"))) cat("state space partitioning: ", attr(x, "part"), " subsets\n")
  cat("total iterations: ", nrow(x), "\n")
  cat("\t(initial) NR iterations:", attr(x, "nnr"), "\n")
  cat("\t(final) MCMC iterations:", nrow(x) - attr(x, "nnr"), "\n")
}

# summary methods
# primary output:
# 1) acceptance rate
# 2) mean relative deviation (if available)
# 3) sample statistics (mean, sd, quantiles, ess, pval) (if available)
summary.sns <- function(object, quantiles = c(0.025, 0.5, 0.975)
  , pval.ref = 0.0, nburnin = max(nrow(object)/2, attr(object, "nnr"))
  , end = nrow(object), thin = 1, ess.method = c("coda", "ise"), ...) {
  K <- ncol(object)
  nnr <- attr(object, "nnr")
  if (nburnin < nnr) warning("it is strongly suggested that burnin period includes NR iterations (which are not valid MCMC iterations)")
  
  # number of subsets in state space partitioning
  npart <- max(1, length(attr(object, "part")))
    
  # average relative deviation of function value from quadratic approximation (post-burnin)
  if (!is.null(attr(object, "reldev"))) reldev.mean <- mean(attr(object, "reldev"), na.rm = TRUE)
  else reldev.mean <- NA
  
  nsmp <- end - nburnin
  if (nsmp > 0) {
    # average acceptance rate for MH transition proposals
    accept.rate <- sum(attr(object, "accept")[nburnin + 1:nsmp, ]) / length(attr(object, "accept")[nburnin + 1:nsmp, ])
    
    myseq <- seq(from = nburnin + 1, to = end, by = thin)
    nseq <- length(myseq)
    
    smp.mean <- colMeans(object[myseq, ])
    smp.sd <- apply(object[myseq, ], 2, sd)
    smp.ess <- ess(object[myseq, ], method = ess.method[1])
    smp.quantiles <- t(apply(object[myseq, ], 2, quantile, probs = quantiles))
    smp.pval <- apply(object[myseq, ], 2, sns.calc.pval, ref = pval.ref, na.rm = FALSE)
    
  } else {
    accept.rate <- NA
    nseq <- 0
    
    smp.mean <- NA
    smp.sd <- NA
    smp.ess <- NA
    smp.quantiles <- NA
    smp.pval <- NA
  }
  ret <- list(K = K, nnr = nnr, nburnin = nburnin, end = end, thin = thin
    , niter = nrow(object), nsmp = nsmp, nseq = nseq, npart = npart
    , accept.rate = accept.rate, reldev.mean = reldev.mean
    , pval.ref = pval.ref, ess.method = ess.method
    , smp = list(mean = smp.mean, sd = smp.sd, ess = smp.ess, quantiles = smp.quantiles, pval = smp.pval))
  class(ret) <- "summary.sns"
  return (ret)
}

print.summary.sns <- function(x, ...) {
  cat("Stochastic Newton Sampler (SNS)\n")
  cat("state space dimensionality: ", x$K, "\n")
  if (x$npart > 1) cat("state space partitioning: ", x$npart, " subsets\n")
  cat("total iterations: ", x$niter, "\n")
  cat("\tNR iterations: ", x$nnr, "\n")
  cat("\tburn-in iterations: ", x$nburnin, "\n")
  cat("\tend iteration: ", x$end, "\n")
  cat("\tthinning interval: ", x$thin, "\n")
  cat("\tsampling iterations (before thinning): ", x$nsmp, "\n")
  #cat("\tsampling iterations (after thinning): ", x$nseq, "\n")
  cat("acceptance rate: ", x$accept.rate, "\n")
  if (!is.na(x$reldev.mean)) cat("\tmean relative deviation from quadratic approx:", format(100*x$reldev.mean, digits=3), "% (post-burnin)\n")
  if (x$nsmp > 0) {
    cat("sample statistics:\n")
    cat("\t(nominal sample size: ", x$nseq, ")\n", sep="")
    stats <- cbind(x$smp$mean, x$smp$sd, x$smp$ess, x$smp$quantiles, x$smp$pval)
    colnames(stats)[c(1:3, 4 + ncol(x$smp$quantiles))] <- c("mean", "sd", "ess", "p-val")
    rownames(stats) <- c(1:x$K)
    printCoefmat(stats[1:min(x$K, 6), ], P.values = TRUE, has.Pvalue = TRUE)
    if (x$K > 6) cat("...\n")
    cat("summary of ess:\n")
    print(summary(x$smp$ess))
  }
}

# plot method
plot.sns <- function(x, nburnin = max(nrow(x)/2, attr(x, "nnr"))
  , select = if (length(x) <= 10) 1:5 else 1, ...) {
  init <- attr(x, "init")
  lp.init <- attr(x, "lp.init")
  lp <- attr(x, "lp")
  
  # in all cases, vertical line delineates transition from nr to mcmc mode
  K <- ncol(x)
  niter <- nrow(x)
  nnr <- attr(x, "nnr")
  if (nburnin < nnr) warning("it is strongly suggested that burnin period includes NR iterations (which are not valid MCMC iterations)")
  
  # log-probability trace plot
  if (1 %in% select) {
    plot(0:niter, c(lp.init, lp), type = "l"
      , xlab = "iter", ylab = "log-probability", main = "Log-Probability Trace Plot")
    if (nnr > 0 && nnr < niter) abline(v = nnr + 0.5, lty = 2, col = "red")
  }
  
  # state vector trace plots
  if (2 %in% select) {
    for (k in 1:K) {
      plot(0:niter, c(init[k], x[, k]), type = "l"
        , xlab = "iter", ylab = paste("x[", k, "]", sep = ""), main = "State Variable Trace Plot")
      if (nnr > 0 && nnr < niter) abline(v = nnr + 0.5, lty = 2, col = "red")
    }
  }
  
  if (nburnin < niter) {

    if (3 %in% select) {
      # effective sample size (horizontal line is maximum possible effective sample size)
      my.ess <- ess(x[(nburnin + 1):niter, ])
      plot(1:K, my.ess, xlab = "k", ylab = "effective sample size", ylim = c(0, niter - nburnin), main = "Effective Sample Size by Coordinate")
      abline(h = niter - nburnin, lty = 2, col = "red")
    }
    
    if (4 %in% select) {
      # state vector (univariate) histograms
      K <- ncol(x)
      for (k in 1:K) {
        hist(x[(nburnin + 1):niter, k], xlab = paste("x[", k, "]", sep = ""), main = "State Variable Histogram (post-burnin)")
        abline(v = mean(x[(nburnin + 1):niter, k]), lty = 2, col = "red")
      }
    }
  
    if (5 %in% select) {
      # state vector (univariate) autocorrelation plots
      K <- ncol(x)
      for (k in 1:K) {
        acf(x[(nburnin + 1):niter, k], xlab = paste("x[", k, "]", sep = ""), main = "State Variable Autocorrelation Plot (post-burnin)")
      }
    }

  }
}

sns.calc.pval <- function(x, ref=0.0, na.rm = FALSE) { # add flag for one-sided vs. two-sided
  if (na.rm) x <- x[!is.na(x)]
  bigger <- median(x)>ref
  if (sd(x)<.Machine$double.eps) {
    ret <- NA
  } else {
    ret <- max(1/length(x), 2*length(which(if (bigger) x<ref else x>ref))/length(x)) # TODO: justify minimum value
  }
  attr(ret, "bigger") <- bigger
  return (ret)
}

