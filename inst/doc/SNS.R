### R code from vignette source 'SNS.Rnw'

###################################################
### code chunk number 1: SNS.Rnw:111-112
###################################################
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: SNS.Rnw:315-318
###################################################
library("sns")
library("mvtnorm")
my.seed <- 0


###################################################
### code chunk number 3: SNS.Rnw:321-328
###################################################
logdensity.mvg <- function(x, mu, isigsq) {
  f <- dmvnorm(x = as.numeric(x),
    mean = mu, sigma = solve(isigsq), log = TRUE)
  g <- - isigsq %*% (x - mu)
  h <- -isigsq
  return (list(f = f, g = g, h = h))
}


###################################################
### code chunk number 4: SNS.Rnw:331-340
###################################################
set.seed(my.seed)
K <- 3
mu <- runif(K, min = -0.5, max = +0.5)
isigsq <- matrix(runif(K*K, min = 0.1, max = 0.2), ncol = K)
isigsq <- 0.5*(isigsq + t(isigsq))
diag(isigsq) <- rep(0.5, K)
x.init <- rep(0.0, K)
x.smp <- sns.run(x.init, logdensity.mvg, niter = 500,
  mh.diag = TRUE, mu = mu, isigsq = isigsq)


###################################################
### code chunk number 5: SNS.Rnw:343-344
###################################################
summary(x.smp)


###################################################
### code chunk number 6: SNS.Rnw:354-361
###################################################
have_regression_factory <- require("RegressionFactory")
if (have_regression_factory) {
  loglike.poisson <- function(beta, X, y) {
    regfac.expand.1par(beta, X = X, y = y,
      fbase1 = fbase1.poisson.log)
  }
}


###################################################
### code chunk number 7: SNS.Rnw:363-367 (eval = FALSE)
###################################################
## loglike.poisson <- function(beta, X, y) {
##   regfac.expand.1par(beta, X = X, y = y,
##     fbase1 = fbase1.poisson.log)
## }


###################################################
### code chunk number 8: SNS.Rnw:370-376
###################################################
set.seed(my.seed)
K <- 5
N <- 1000
X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
beta <- runif(K, -0.5, +0.5)
y <- rpois(N, exp(X %*% beta))


###################################################
### code chunk number 9: SNS.Rnw:379-382
###################################################
beta.init <- rep(0.0, K)
beta.glm <- glm(y ~ X - 1, family = "poisson",
  start = beta.init)$coefficients


###################################################
### code chunk number 10: SNS.Rnw:385-391
###################################################
if (have_regression_factory) {
  beta.sns <- sns.run(beta.init, fghEval = loglike.poisson,
    niter = 20, nnr = 20, X = X, y = y)
  beta.nr <- beta.sns[20, ]
  cbind(beta.glm, beta.nr)
}


###################################################
### code chunk number 11: SNS.Rnw:393-397 (eval = FALSE)
###################################################
## beta.sns <- sns.run(beta.init, fghEval = loglike.poisson,
##   niter = 20, nnr = 20, X = X, y = y)
## beta.nr <- beta.sns[20, ]
## cbind(beta.glm, beta.nr)


###################################################
### code chunk number 12: SNS.Rnw:400-404
###################################################
if (have_regression_factory) {
  beta.smp <- sns.run(beta.init, loglike.poisson
    , niter = 200, nnr = 20, mh.diag = TRUE, X = X, y = y)
}


###################################################
### code chunk number 13: SNS.Rnw:406-408 (eval = FALSE)
###################################################
## beta.smp <- sns.run(beta.init, loglike.poisson
##   , niter = 200, nnr = 20, mh.diag = TRUE, X = X, y = y)


###################################################
### code chunk number 14: lp_plot
###################################################
if (have_regression_factory) {
  plot(beta.smp, select = 1)
} else {
  plot(1:10)
}


###################################################
### code chunk number 15: SNS.Rnw:418-419 (eval = FALSE)
###################################################
## plot(beta.smp, select = 1)


###################################################
### code chunk number 16: fig1
###################################################
if (have_regression_factory) {
  plot(beta.smp, select = 1)
} else {
  plot(1:10)
}


###################################################
### code chunk number 17: SNS.Rnw:429-432
###################################################
if (have_regression_factory) {
  summary(beta.smp)
}


###################################################
### code chunk number 18: SNS.Rnw:434-435 (eval = FALSE)
###################################################
## summary(beta.smp)


###################################################
### code chunk number 19: SNS.Rnw:440-445
###################################################
if (have_regression_factory) {
  beta.smp <- sns.run(beta.init, loglike.poisson,
    niter = 1000, nnr = 20, mh.diag = TRUE, X = X, y = y)
  predmean.poisson <- function(beta, Xnew) exp(Xnew %*% beta)
}


###################################################
### code chunk number 20: SNS.Rnw:447-450 (eval = FALSE)
###################################################
## beta.smp <- sns.run(beta.init, loglike.poisson,
##   niter = 1000, nnr = 20, mh.diag = TRUE, X = X, y = y)
## predmean.poisson <- function(beta, Xnew) exp(Xnew %*% beta)


###################################################
### code chunk number 21: SNS.Rnw:453-457
###################################################
if (have_regression_factory) {
  ymean.new <- predict(beta.smp, predmean.poisson,
    nburnin = 100, Xnew = X)
}


###################################################
### code chunk number 22: SNS.Rnw:459-461 (eval = FALSE)
###################################################
## ymean.new <- predict(beta.smp, predmean.poisson,
##   nburnin = 100, Xnew = X)


###################################################
### code chunk number 23: SNS.Rnw:466-472
###################################################
if (have_regression_factory) {
  predsmp.poisson <- function(beta, Xnew)
    rpois(nrow(Xnew), exp(Xnew %*% beta))
  ysmp.new <- predict(beta.smp, predsmp.poisson
    , nburnin = 100, Xnew = X)
}


###################################################
### code chunk number 24: SNS.Rnw:474-478 (eval = FALSE)
###################################################
## predsmp.poisson <- function(beta, Xnew)
##   rpois(nrow(Xnew), exp(Xnew %*% beta))
## ysmp.new <- predict(beta.smp, predsmp.poisson
##   , nburnin = 100, Xnew = X)


###################################################
### code chunk number 25: SNS.Rnw:481-484
###################################################
if (have_regression_factory) {
  summary(ymean.new)
}


###################################################
### code chunk number 26: SNS.Rnw:486-487 (eval = FALSE)
###################################################
## summary(ymean.new)


###################################################
### code chunk number 27: SNS.Rnw:493-496
###################################################
if (have_regression_factory) {
  summary(ysmp.new)
}


###################################################
### code chunk number 28: SNS.Rnw:498-499 (eval = FALSE)
###################################################
## summary(ysmp.new)


###################################################
### code chunk number 29: SNS.Rnw:511-522
###################################################
if (have_regression_factory) {
  set.seed(my.seed)
  K <- 100
  X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
  beta <- runif(K, -0.5, +0.5)
  y <- rpois(N, exp(X %*% beta))
  beta.init <- glm(y ~ X - 1, family = "poisson")$coefficients
  beta.smp <- sns.run(beta.init, loglike.poisson,
    niter = 100, nnr = 10, mh.diag = TRUE, X = X, y = y)
  summary(beta.smp)
}


###################################################
### code chunk number 30: SNS.Rnw:524-533 (eval = FALSE)
###################################################
## set.seed(my.seed)
## K <- 100
## X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
## beta <- runif(K, -0.5, +0.5)
## y <- rpois(N, exp(X %*% beta))
## beta.init <- glm(y ~ X - 1, family = "poisson")$coefficients
## beta.smp <- sns.run(beta.init, loglike.poisson,
##   niter = 100, nnr = 10, mh.diag = TRUE, X = X, y = y)
## summary(beta.smp)


###################################################
### code chunk number 31: SNS.Rnw:541-547
###################################################
if (have_regression_factory) {
  beta.smp.part <- sns.run(beta.init, loglike.poisson,
    niter = 100, nnr = 10, mh.diag = TRUE,
    part = sns.make.part(K, 10), X = X, y = y)
  summary(beta.smp.part)
}


###################################################
### code chunk number 32: SNS.Rnw:549-553 (eval = FALSE)
###################################################
## beta.smp.part <- sns.run(beta.init, loglike.poisson,
##   niter = 100, nnr = 10, mh.diag = TRUE,
##   part = sns.make.part(K, 10), X = X, y = y)
## summary(beta.smp.part)


###################################################
### code chunk number 33: ssp_plot
###################################################
if (have_regression_factory) {
  par(mfrow = c(1,2))
  plot(beta.smp, select = 1)
  plot(beta.smp.part, select = 1)
} else {
  plot(1:10)
}


###################################################
### code chunk number 34: SNS.Rnw:570-573 (eval = FALSE)
###################################################
## par(mfrow = c(1,2))
## plot(beta.smp, select = 1)
## plot(beta.smp.part, select = 1)


###################################################
### code chunk number 35: fig1
###################################################
if (have_regression_factory) {
  par(mfrow = c(1,2))
  plot(beta.smp, select = 1)
  plot(beta.smp.part, select = 1)
} else {
  plot(1:10)
}


###################################################
### code chunk number 36: SNS.Rnw:591-603
###################################################
loglike.linreg.het <- function(coeff, X, Z, y) {
  K1 <- ncol(X)
  K2 <- ncol(Z)
  beta <- coeff[1:K1]
  gamma <- coeff[K1 + 1:K2]
  
  mu <- X %*% beta
  sigma <- sqrt(exp(Z %*% gamma))
  f <- sum(dnorm(y, mu, sigma, log = TRUE))
  
  return (f)
}


###################################################
### code chunk number 37: SNS.Rnw:606-617
###################################################
set.seed(my.seed)
K1 <- 5
K2 <- 5
N <- 1000
X <- matrix(runif(N * K1, -0.5, +0.5), ncol = K1)
Z <- matrix(runif(N * K2, -0.5, +0.5), ncol = K2)
beta <- runif(K1, -0.5, +0.5)
gamma <- runif(K1, -0.5, +0.5)
mu <- X %*% beta
var <- exp(Z %*% gamma)
y <- rnorm(N, X %*% beta, sd = sqrt(var))


###################################################
### code chunk number 38: SNS.Rnw:620-625
###################################################
coeff.init <- rep(0.0, K1 + K2)
check.logdensity <- sns.check.logdensity(coeff.init, loglike.linreg.het
  , X = X, Z = Z, y = y, dx = 1, nevals = 10
  , blocks = list(1:(K1+K2), 1:K1, K1 + 1:K2))
check.logdensity


###################################################
### code chunk number 39: SNS.Rnw:649-665
###################################################
if (have_regression_factory) {
  loglike.linreg.het.beta <- function(beta, gamma, X, Z, y) {
    K1 <- length(beta)
    ret <- regfac.expand.2par(c(beta, gamma), X, Z, y
      , fbase2 = fbase2.gaussian.identity.log)
    return (list(f = ret$f, g = ret$g[1:K1], h = ret$h[1:K1, 1:K1]))
  }
  loglike.linreg.het.gamma <- function(gamma, beta, X, Z, y) {
    K1 <- length(beta)
    K2 <- length(gamma)
    ret <- regfac.expand.2par(c(beta, gamma), X, Z, y
      , fbase2 = fbase2.gaussian.identity.log)
    return (list(f = ret$f, g = ret$g[K1 + 1:K2]
             , h = ret$h[K1 + 1:K2, K1 + 1:K2]))
  }
}


###################################################
### code chunk number 40: SNS.Rnw:667-681 (eval = FALSE)
###################################################
## loglike.linreg.het.beta <- function(beta, gamma, X, Z, y) {
##   K1 <- length(beta)
##   ret <- regfac.expand.2par(c(beta, gamma), X, Z, y
##     , fbase2 = fbase2.gaussian.identity.log)
##   return (list(f = ret$f, g = ret$g[1:K1], h = ret$h[1:K1, 1:K1]))
## }
## loglike.linreg.het.gamma <- function(gamma, beta, X, Z, y) {
##   K1 <- length(beta)
##   K2 <- length(gamma)
##   ret <- regfac.expand.2par(c(beta, gamma), X, Z, y
##     , fbase2 = fbase2.gaussian.identity.log)
##   return (list(f = ret$f, g = ret$g[K1 + 1:K2]
##            , h = ret$h[K1 + 1:K2, K1 + 1:K2]))
## }


###################################################
### code chunk number 41: SNS.Rnw:684-703
###################################################
if (have_regression_factory) {
  nsmp <- 100
  beta.iter <- rep(0.0, K1)
  gamma.iter <- rep(0.0, K2)
  beta.smp <- array(NA, dim = c(nsmp, K1))
  gamma.smp <- array(NA, dim = c(nsmp, K1))
  for (n in 1:nsmp) {
    beta.iter <- sns(beta.iter, loglike.linreg.het.beta
      , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
    gamma.iter <- sns(gamma.iter, loglike.linreg.het.gamma
      , beta = beta.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
    beta.smp[n, ] <- beta.iter
    gamma.smp [n, ] <- gamma.iter
  }
  beta.est <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
  gamma.est <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
  print(cbind(beta, beta.est))
  print(cbind(gamma, gamma.est))
}


###################################################
### code chunk number 42: SNS.Rnw:705-722 (eval = FALSE)
###################################################
## nsmp <- 100
## beta.iter <- rep(0.0, K1)
## gamma.iter <- rep(0.0, K2)
## beta.smp <- array(NA, dim = c(nsmp, K1))
## gamma.smp <- array(NA, dim = c(nsmp, K1))
## for (n in 1:nsmp) {
##   beta.iter <- sns(beta.iter, loglike.linreg.het.beta
##     , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
##   gamma.iter <- sns(gamma.iter, loglike.linreg.het.gamma
##     , beta = beta.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
##   beta.smp[n, ] <- beta.iter
##   gamma.smp [n, ] <- gamma.iter
## }
## beta.est <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
## gamma.est <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
## cbind(beta, beta.est)
## cbind(gamma, gamma.est)


###################################################
### code chunk number 43: SNS.Rnw:725-746
###################################################
have_mfu <- require("MfUSampler")
if (have_mfu && have_regression_factory) {
  loglike.linreg.het.gamma.fonly <- function(gamma, beta, X, Z, y) {
    return (regfac.expand.2par(c(beta, gamma), X, Z, y
      , fbase2 = fbase2.gaussian.identity.log, fgh = 0))
  }
  beta.iter <- rep(0.0, K1)
  gamma.iter <- rep(0.0, K2)
  for (n in 1:nsmp) {
    beta.iter <- sns(beta.iter, loglike.linreg.het.beta
      , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
    gamma.iter <- MfU.Sample(gamma.iter, loglike.linreg.het.gamma.fonly
      , beta = beta.iter, X = X, Z = Z, y = y)
    beta.smp[n, ] <- beta.iter
    gamma.smp [n, ] <- gamma.iter
  }
  beta.est.mix <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
  gamma.est.mix <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
  print(cbind(beta, beta.est.mix))
  print(cbind(gamma, gamma.est.mix))
}


###################################################
### code chunk number 44: SNS.Rnw:748-766 (eval = FALSE)
###################################################
## loglike.linreg.het.gamma.fonly <- function(gamma, beta, X, Z, y) {
##   return (regfac.expand.2par(c(beta, gamma), X, Z, y
##     , fbase2 = fbase2.gaussian.identity.log, fgh = 0))
## }
## beta.iter <- rep(0.0, K1)
## gamma.iter <- rep(0.0, K2)
## for (n in 1:nsmp) {
##   beta.iter <- sns(beta.iter, loglike.linreg.het.beta
##     , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
##   gamma.iter <- MfU.Sample(gamma.iter, loglike.linreg.het.gamma.fonly
##     , beta = beta.iter, X = X, Z = Z, y = y)
##   beta.smp[n, ] <- beta.iter
##   gamma.smp [n, ] <- gamma.iter
## }
## beta.est.mix <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
## gamma.est.mix <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
## cbind(beta, beta.est.mix)
## cbind(gamma, gamma.est.mix)


###################################################
### code chunk number 45: SNS.Rnw:770-788
###################################################
if (have_regression_factory && have_mfu) {
  loglike.linreg.het.gamma.numaug <- 
    sns.fghEval.numaug(loglike.linreg.het.gamma.fonly, numderiv = 2)
  beta.iter <- rep(0.0, K1)
  gamma.iter <- rep(0.0, K2)
  for (n in 1:nsmp) {
    beta.iter <- sns(beta.iter, loglike.linreg.het.beta
      , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
    gamma.iter <- sns(gamma.iter, loglike.linreg.het.gamma
      , beta = beta.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
    beta.smp[n, ] <- beta.iter
    gamma.smp [n, ] <- gamma.iter
  }
  beta.est.num <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
  gamma.est.num <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
  print(cbind(beta, beta.est.num))
  print(cbind(gamma, gamma.est.num))
}


###################################################
### code chunk number 46: SNS.Rnw:790-806 (eval = FALSE)
###################################################
## loglike.linreg.het.gamma.numaug <- 
##   sns.fghEval.numaug(loglike.linreg.het.gamma.fonly, numderiv = 2)
## beta.iter <- rep(0.0, K1)
## gamma.iter <- rep(0.0, K2)
## for (n in 1:nsmp) {
##   beta.iter <- sns(beta.iter, loglike.linreg.het.beta
##     , gamma = gamma.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
##   gamma.iter <- sns(gamma.iter, loglike.linreg.het.gamma
##     , beta = beta.iter, X = X, Z = Z, y = y, rnd = nsmp>10)
##   beta.smp[n, ] <- beta.iter
##   gamma.smp [n, ] <- gamma.iter
## }
## beta.est.num <- colMeans(beta.smp[(nsmp/2+1):nsmp, ])
## gamma.est.num <- colMeans(gamma.smp[(nsmp/2+1):nsmp, ])
## cbind(beta, beta.est.num)
## cbind(gamma, gamma.est.num)


###################################################
### code chunk number 47: SNS.Rnw:875-876
###################################################
sessionInfo()


