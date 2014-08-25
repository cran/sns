######################################
#                                    #
#  Methods for sns objects           #
#                                    #
######################################

# print method 
print.sns <- function(x, ...)
{
  cat(paste0("Acceptance rate % = ", x$acceptance, "\n"))
  eff.samp <- round(ess(x$samplesMat)) 
  cat(paste0("Effective (independent) samples from ", nrow(x$samplesMat)),
      " draws:\n\tmin    - ", min(eff.samp), "\n\tmedian - ", median(eff.samp),
      "\n\tmax    - ", max(eff.samp), "\n")
}

# summary method 
summary.sns <- function(object, show.means = FALSE, ...)
{
  cat("-----------------------\n")
  cat("MCMC sampling using SNS\n")
  cat("-----------------------\n")
  print(object)
  cat(paste0("No. of burn-in iterations = ", object$burn.iters),"\nTimings:\n")
  cat(paste0("\tburn-in : ", round(object$burnin.time, 1), " sec\n"))
  cat(paste0("\tsampling: ", round(object$sample.time, 1), " sec\n"))
  if (show.means) {
    cat("\n")
    stats <- cbind(colMeans(object$samplesMat), 
                   apply(object$samplesMat, 2, sd))
    colnames(stats) <- c("sample mean", "std dev")
    rownames(stats) <- c(1:ncol(object$samplesMat))
    printCoefmat(stats)
  }
  cat("-----------------------\n")
}
