## Hessian, gradient expansion functions for 1-parameter base distributions
#  Useful for GLM and GLM-like models
#  Reference: Theorem A.1 of http://arxiv.org/pdf/1308.0657v1.pdf

### One-parameter distributions
ll_1par_expand <- function(beta,X,y,ll_1par_base) {
  # obtain base distribution derivatives
  ret <- ll_1par_base(X%*%beta,y) 

  # expand base derivatives
  f <- sum(ret$f)
  g <- t(X)%*%ret$g
	xtw <- 0*X
	for (k in 1:ncol(X)) xtw[,k] <- X[,k]*ret$h
	h <- t(xtw)%*%X
  
  return (list(f=f,g=g,h=h))
}

## log-likelihood for base GLM distributions
# Bernoulli distribution with logit link function: y ~ dBern(1/(1+exp(-u))), y=0,1
ll_bern_logit <- function(u,y) {
  eu <- exp(u)
  f <- -log(1+1/eu)-(1-y)*u
  g <- 1/(1+eu)-(1-y)
  h <- -eu/(1+eu)^2
  return (list(f=f,g=g,h=h))
}
# Poisson distribution with log link function: y ~ dPois(exp(u)), y=0,1,2,...
ll_pois_log <- function(u,y) {
  eu <- exp(u)
  f <- y*u-eu - lfactorial(y)
  g <- y-eu
  h <- -eu
  return (list(f=f,g=g,h=h))
}
# exponential distribution with log link function: y ~ dExp(exp(u)), y in [0,+Inf] 
ll_exp_log <- function(u,y) {
  eu <- exp(u)
  f <- u-y*eu
  g <- 1-y*eu
  h <- -y*eu
  return (list(f=f,g=g,h=h))
}
# goemetric distribution with logit link function: y ~ dGeom(1/(1+exp(-u)))
ll_geom_logit <- function(u,y) {
  eu <- exp(u)
  f <- -(y*u+(1+y)*log(1+1/eu))
  g <- -y+(1+y)/(1+eu)
  h <- -(1+y)*eu/(1+eu)^2
  return (list(f=f,g=g,h=h))
}
