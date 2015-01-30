### R code from vignette source 'SNS.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: SNS.Rnw:108-109
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: SNS.Rnw:214-216
###################################################
library(sns)
library(mvtnorm)


###################################################
### code chunk number 3: SNS.Rnw:221-228
###################################################
logdensity.mvg <- function(x, mu, isigsq) {
  f <- dmvnorm(x = as.numeric(x)
    , mean = mu, sigma = solve(isigsq), log = TRUE)
  g <- - isigsq %*% (x - mu)
  h <- -isigsq
  return (list(f = f, g = g, h = h))
}


###################################################
### code chunk number 4: SNS.Rnw:231-239
###################################################
K <- 3
mu <- runif(K, min = -0.5, max = +0.5)
isigsq <- matrix(runif(K*K, min = 0.1, max = 0.2), ncol = K)
isigsq <- 0.5*(isigsq + t(isigsq))
diag(isigsq) <- rep(0.5, K)
x.init <- rep(0.0, K)
x.smp <- sns.run(x.init, logdensity.mvg, niter = 500
  , mh.diag = TRUE, mu = mu, isigsq = isigsq)


###################################################
### code chunk number 5: SNS.Rnw:242-243
###################################################
summary(x.smp)


###################################################
### code chunk number 6: SNS.Rnw:253-258
###################################################
library(RegressionFactory)
loglike.poisson <- function(beta, X, y) {
  regfac.expand.1par(beta, X = X, y = y
    , fbase1 = fbase1.poisson.log)
}


###################################################
### code chunk number 7: SNS.Rnw:261-266
###################################################
K <- 5
N <- 1000
X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
beta <- runif(K, -0.5, +0.5)
y <- rpois(N, exp(X %*% beta))


###################################################
### code chunk number 8: SNS.Rnw:269-271
###################################################
beta.init <- rep(0.0, K)
beta.glm <- glm(y ~ X - 1, family = "poisson", start = beta.init)$coefficients


###################################################
### code chunk number 9: SNS.Rnw:274-278
###################################################
beta.sns <- sns.run(beta.init, fghEval = loglike.poisson
  , niter = 20, nnr = 20, X = X, y = y)
beta.nr <- beta.sns[20, ]
cbind(beta.glm, beta.nr)


###################################################
### code chunk number 10: SNS.Rnw:281-283
###################################################
beta.smp <- sns.run(beta.init, loglike.poisson
  , niter = 200, nnr = 20, mh.diag = TRUE, X = X, y = y)


###################################################
### code chunk number 11: lp_plot
###################################################
plot(beta.smp, select = 1)


###################################################
### code chunk number 12: fig1
###################################################
plot(beta.smp, select = 1)


###################################################
### code chunk number 13: SNS.Rnw:299-300
###################################################
summary(beta.smp)


###################################################
### code chunk number 14: SNS.Rnw:305-308 (eval = FALSE)
###################################################
## beta.smp <- sns.run(beta.init, loglike.poisson
##   , niter = 1000, nnr = 20, mh.diag = TRUE, X = X, y = y)
## predmean.poisson <- function(beta, Xnew) exp(Xnew %*% beta)


###################################################
### code chunk number 15: SNS.Rnw:311-313 (eval = FALSE)
###################################################
## ymean.new <- predict(beta.smp, predmean.poisson
##   , nburnin = 100, Xnew = X)


###################################################
### code chunk number 16: SNS.Rnw:318-322 (eval = FALSE)
###################################################
## predsmp.poisson <- function(beta, Xnew)
##   rpois(nrow(Xnew), exp(Xnew %*% beta))
## ysmp.new <- predict(beta.smp, predsmp.poisson
##   , nburnin = 100, Xnew = X)


###################################################
### code chunk number 17: SNS.Rnw:325-327 (eval = FALSE)
###################################################
## summary(ymean.new)
## summary(ysmp.new)


###################################################
### code chunk number 18: SNS.Rnw:329-332
###################################################
load("summs.pred")
print(summ.ymean)
print(summ.ysmp)


###################################################
### code chunk number 19: SNS.Rnw:339-347 (eval = FALSE)
###################################################
## K <- 100
## X <- matrix(runif(N * K, -0.5, +0.5), ncol = K)
## beta <- runif(K, -0.5, +0.5)
## y <- rpois(N, exp(X %*% beta))
## beta.init <- glm(y ~ X - 1, family = "poisson")$coefficients
## beta.smp <- sns.run(beta.init, loglike.poisson
##   , niter = 100, nnr = 10, mh.diag = TRUE, X = X, y = y)
## summary(beta.smp)


###################################################
### code chunk number 20: SNS.Rnw:349-351
###################################################
load("summs")
print(summ.nopart)


###################################################
### code chunk number 21: SNS.Rnw:355-359 (eval = FALSE)
###################################################
## beta.smp.part <- sns.run(beta.init, loglike.poisson
##   , niter = 100, nnr = 10, mh.diag = TRUE
##   , part = sns.make.part(K, 10), X = X, y = y)
## summary(beta.smp.part)


###################################################
### code chunk number 22: SNS.Rnw:361-362
###################################################
print(summ.part)


###################################################
### code chunk number 23: SNS.Rnw:366-369 (eval = FALSE)
###################################################
## par(mfrow = c(1,2))
## plot(beta.smp, select = 1)
## plot(beta.smp.part, select = 1)


