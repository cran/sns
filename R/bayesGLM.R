# A closure which generates a function, gradient, Hessian evaluator
# for log-likelihood functions of 4 types of GLMs.
#
# Args:
#   N - number of observations
#   K - number of variables
#   glmtype - must be one of the 4 strings
#   X - data matrix of the explanatory variable, must have N rows & K cols
#   y - vector having the dependent vaiable (reponse) corresponging to X
# Output:
#   Am evaluator which can be passes to sns()
glmfgh <- function(N, K, glmtype = "logistic", X=NULL, y=NULL) 
{
  # Check user supplied data for consistency
  if (!is.null(X)){
      stopifnot(!is.matrix(X) || ncol(X) != K || nrow(X) != N)
      stopifnot(!is.null(y) || length(y) != N)
  }

  # Set the link & distribution for specified GLM
  if (glmtype == "logistic")
      base.func <- ll_bern_logit
  else if (glmtype == "poisson")
      base.func <- ll_pois_log
  else if (glmtype == "exponential")
      base.func <- ll_exp_log
  else if (glmtype == "geometric")
      base.func <- ll_geom_logit
  else
      stop(paste("Unrecognized Arg 'glmtype'. Must be one of:",
                 "logistic, poisson, exponential, geomteric."))

  if (is.null(X) || is.null(y)) {
      # Generate simulated data
      stopifnot(N > K)
      data <- simGLM(N, K, glmtype = glmtype) 
      model.coef <- data$coef
      y <- data$df$y
      data$df$y <- NULL
      X <- data.matrix(data$df)
      data$df$y <- y 
  }

  # function, gradient, Hessian evaluator
  sns.fghEval <- function(beta, ...)
  {
      out <- ll_1par_expand(beta, X, y, base.func)
      return(list("f" = out$f, "g" = out$g, "h" = out$h))
  }
 
  return(sns.fghEval)
}

# Simulates data for a generalized linear model
# N - number of observations
# K - number of variables
# glmtype - must be one of: logistic, poisson, exponential, geomteric
simGLM <- function(N, K, glmtype="logistic")
{
  set.seed(exp(1));
  stopifnot(K >= 1)
  stopifnot(N > K)

  # Generate data matrix
  sim.df <-  matrix(c(rep(1, N), runif(N * (K - 1), min = -1, max = 1)), 
                    nrow = N, ncol = K)
  # Generate coefficient vec 
  coeff <- runif(K, -.5, .5)
 
  # Assign var names, first data var is the Intercept 
  colnames(sim.df) <- paste0("x", seq(1, K))

  # Compute linear predictors 
  eta <- sim.df %*% matrix(coeff, nrow = length(coeff), ncol = 1)
  
  # Combine into a data.frame
  sim.df <- data.frame(sim.df)
 
  # Compute response
  if (glmtype == "logistic")
      sim.df$y <- as.numeric(sapply(eta, function(x) 1/(1+exp(-x))) >= runif(N))
  else if (glmtype == "poisson")
      sim.df$y <- rpois(n=N, lambda = exp(eta))
  else if (glmtype == "exponential")
      sim.df$y <- rexp(n=N, rate = exp(eta))
  else if (glmtype == "geometric")
      sim.df$y <- rgeom(n=N, prob = sapply(eta, function(x) 1/(1+exp(-x))))
  else
      stop(paste("Unrecognized Arg 'glmtype'. Must be one of:", 
                 "logistic, poisson, exponential, geomteric."))

  # Return data and coefficients
  output <- list(df = sim.df, coef = coeff, glmtype = glmtype)
}

