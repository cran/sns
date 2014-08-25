###############################################################################
#                        Package: sns                                         #
#                                                                             #
# Stochastic Newton Sampler (SNS)- implements the MH-MGT (Metropolis-Hastings #
# with Multivariate Gaussian Tangents) algorithm described in the preprint:   #
#                        http://arxiv.org/abs/1308.0657                       #
# Draws samples from twice differentiable, log-concave pdf.                   #
#                                                                             #
# Version: 1.0                                                                #
#                                                                             #
#               Scientific Computing Group, Sentrana Inc.                     #
###############################################################################

###############################################################################
# 	Core sampling function: draws samples from a multivariate pdf
# Args:
#   init - starting point for the Markov chain
#   f    - function, gradient, Hessian evaluator for the log-density. 
#          Must a return a list with labels:
#          f - the log-probability density
#          g - gradient vector 
#          h - Hessian matrix
#   rnd  - Runs 1 iteration of Newton's method (non-stochastic) when FALSE
#          Runs Metropolis-Hastings for draw a sample when TRUE
#          NOTE: Set to FALSE only during burn-in 
#   gfit - Gaussian fit at 'init'. If NULL, Gaussian fit at 'init' is computed
#   ...  - Extra args all passed to evaluator whenever it's called
#
# Output:
#   The sample, drawn from the pdf, as a vector, with attributes:
#     accept - TRUE/FALSE, specifying whether the Metropolis move was accepted
#     ll     - value of the function at the sampled point
#     gfit   - Gaussian fit at the sampled point
###############################################################################
sns <- function(init, fghEval, rnd=TRUE, gfit=NULL, ...)
{ 
  f <- fghEval
  x <- init

  fitGaussian <- function(x, f, ...) 
  {
    ret <- f(x,...)                # Evaluate the function at 'x'
    Sigma <- solve(-ret$h)           
    mu <- x + Sigma %*% ret$g

    return (list(mu=as.vector(mu), # Newton method solution
               Sigma=Sigma,        # Inverse Hessian or Covariance matrix
               iSigma=-ret$h,      # Inverse covariance or Hessian
               f=ret$f,            # function value 
               g=ret$g))           # gradient 
  }

  # rnd: if FALSE, perform Newton's optimization (non-stochastic)
  # Fit Gaussian at x
  if (is.null(gfit)) gfit <- fitGaussian(x = x, f = f, ...)
  mu     <- gfit$mu 
  Sigma  <- gfit$Sigma   # Covariance
  iSigma <- gfit$iSigma  # Inverse covariance

  K <- length(x);
      
  if (rnd) {
    # Draw sample from proposal distribution (Gaussian fit at x)
    x.prop <- as.vector(rmvnorm(n=1, mean=mu, sigma=Sigma))
  } else {
    # Run (non-stochastic) Newton optimization 
    rho <- 0.5; c <- 0.5;
    alphak <- 1; 
    d <- mu - x; # use newton's direction as step
    search_x <- as.vector(mu);
        
    fk <- gfit$f; # Values at the current point
    gk <- gfit$g; 
    fk1 <- f(search_x, ...)$f; # Function value at searching point
    ls_iter <- 1;
    # Linesearch by backtracking from full Newton step
    while (fk1 < fk + c*alphak*(t(gk)%*%d) && ls_iter < 20) { 
      alphak <- alphak*rho; # if so, then go half way
      search_x <- x + alphak*d;
      fk1 <- f(search_x, ...)$f;
      ls_iter <- ls_iter + 1;
    }
    x.prop <- as.vector(search_x);
  }

  log.q.prop <- dmvnorm(as.vector(x.prop), mu, Sigma, log=TRUE)
  
  # fit Gaussian at x.prop
  gfit.prop <- fitGaussian(x=x.prop,f=f,...)
  mu.prop <- gfit.prop$mu
  Sigma.prop <- gfit.prop$Sigma
  iSigma.prop <- gfit.prop$iSigma

  # create MH acceptance ratio
  log.q <- dmvnorm(as.vector(x), mu.prop, Sigma.prop, log=TRUE)
  
  log.p <- gfit$f
  log.p.prop <- gfit.prop$f
  log.ratio <- (log.p.prop-log.p) + (log.q-log.q.prop)
  ratio <- min(1,exp(log.ratio))
  
  # perform acceptance test
  if (ratio==1 || runif(1)<ratio || !rnd) {
	 gfit <- gfit.prop
	 x <- x.prop;
   attr(x,"sample") <- x.prop
	 attr(x,"accept") <- TRUE
	 attr(x,"ll") <- log.p.prop
  } else {
	 attr(x,"accept") <- FALSE
	 attr(x,"ll") <- log.p
  }
  attr(x,"gfit") <- gfit
  return (x)
}

###############################################################################
#                 Main user function
# Args:
#   K        - dimension of the space to draw samples from 
#   nburnin  - number of burn-in iteration (non-stochastic, Newton-Raphson)
#   nsample  - number of samples to draw (after burn-in)
#   fghEval  - function, gradient, Hessian evaluator for the log-density. 
#              Must a return a list with labels:
#                f - the log-probability density
#                g - gradient vector 
#                h - Hessian matrix
#   start    - initial point for the Markov chain. Default: rep(0.1, K)          
#   print.level - if non zero, prints sampling progress
#   report.progress - number of sampling iterations between printing progress 
#   ...  - Extra args all passed to evaluator whenever it's called
# 
# Output:
#   An object of class sns
#
# Note: 
#   The sampler is a Metropolis-Hastings Markov chain Monte Carlo variant, with
#   a special form of the proposal function. During burn-in, a non-stochastic
#   Newton-Raphson optimization is performed to get close to the pdf's mode.
#   
#   Currently restricted to log-concave, twice differentiable densities.   
###############################################################################
sns.run <- function(K, nburnin, nsample, fghEval, start=NULL, print.level=0, 
                    report.progress=100, ...)
{
  if (report.progress <= 0) {
      warning("Invalid value specifiec for 'report.progress', using default.")
      report.progress <- 100
  }
  if (is.null(start)) start <- rep(0.1, K) 
  if (!is.null(start) && length(start) != K)
      stop("Mismatch between args 'K' and 'start'")

  # Burn In iterations
  sample <- start
  t0 <- proc.time()
  for (i in 1:nburnin) {
      sample <- sns(sample, fghEval, rnd = FALSE)
  }
  t1 <- proc.time()
  burninTime <- as.numeric(t1 - t0)[3]
  if (print.level)
      cat(paste0("Finished ", nburnin, " burn-in iterations.\n"))

  # MCMC sampling
  acceptCnt <- 0
  chain <- matrix( , nrow=nsample, ncol=K)
  chain[1, ] <- attr(sample, "sample")
  t1 <- proc.time()
  for (i in 2:nsample) {
      sample <- sns(sample, fghEval)
      if (attr(sample, "accept")) acceptCnt <- acceptCnt + 1
      chain[i, ] <- attr(sample, "sample")
      if (print.level && (i %% report.progress == 0))
          cat(paste0("Finished  ", i, " sampling iterations out of ", nsample, ".\n"))
  }
  t2 <- proc.time()
  sampleTime <- as.numeric(t2 - t1)[3]
  acceptRate <- acceptCnt * 100 / nsample

  return(structure(list(
           samplesMat = chain,
           acceptance = acceptRate,
           burn.iters = nburnin,
           sample.time= sampleTime,
           burnin.time= burninTime),
           class = "sns")) 
}
