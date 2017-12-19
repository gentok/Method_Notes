# edvreg.R

# http://www.sscnet.ucla.edu/polisci/faculty/lewis/software/edvreg.R

# A function to fit the estimated dependent variable regression described in
# Jeffrey B. Lewis and Drew A. Linzer (2005), ``Estimating Regression Models in Which the Dependent Variable
# Is Based on Estimates.'' _Political_Analysis_ 13(4):345-364. Based on code written by Drew Linzer
# in March of 2005 ( MonteCarlo_EDV.R).
#
# Jeff Lewis (2013)
#

library(MASS)   # needed to use ginv(), the generalized inverse of a matrix
library(car)    # needed to access a fast wcrossprod function

# Find the trace of a square matrix
tr <- function(x) { sum(diag(x)) }

#
#  edvreg
#     mod: Regression formula
#
#     omegasq: Nx1 vector of estimated sampling/error variances of the dependent variable
#              (Note that unlike the STATA code, what is called for here is the sample variance not the
#              standard error)
#
#     proportional: "True" if omegasq is only proportional to the sampling/error variances of
#                   of the dependent variables
#
# Example call: res <- edvreg( y ~ x, omegasq = est.var.y )
#
edvreg <- function(mod,dt,omegasq,proportional=FALSE) {

  #**NEW**#
  vn <- all.vars(mod)
  dtx <- dt[,vn]
  omegasq <- omegasq[complete.cases(dtx)]
  dtx <- dtx[complete.cases(dtx),]
  attach(dtx)

  # Some housekeeping stolen from John Fox's corchrane-orcutt implementation
  X <- model.matrix(mod)
  xnames <- colnames(X)
  y <- model.response(model.frame(mod))

  # When omegasq is known up to a proportion (for example 1/N_i) ...
  if (proportional) {
    residsq <- (resid(lm(y~X-1)))^2
    steptwo <- lm(residsq~omegasq)
    zeroeps <- FALSE
    if (coef(steptwo)[1] < 0) {
      warning("Variance of epsilon estimated to be negative! Setting to 0...")
      steptwo <- lm(residsq~0+omegasq)
      zeroeps <- TRUE
    }
    predval <- fitted.values(steptwo)
    if (sum(predval<0) < 0) {
      warning("Negative omegaSq estimates encountered, setting to 0.0001")
      predval[predval<0] <- 0.0001 # error check against very rare possibility of negative weight
    }
    model.FGLS <- lm(y~X-1,weights=(1/predval))
    model.FGLS$vareps <- ifelse(zeroeps,0,steptwo$coefficients[1])
    model.FGLS$omegasqFactor <- ifelse(zeroeps,steptwo$coefficients[1],steptwo$coefficients[2])
    model.FGLS$predval <- predval
  }

  # When omegasq is known exactly...
  else {
    N <- nrow(X)
    sigmahatsq <- ( sum((residuals(lm(y~X-1)))^2) - sum(omegasq) +  tr(ginv(crossprod(X)) %*% wcrossprod(X,X,w=omegasq)))/(N-ncol(X)-1)
    if (sigmahatsq < 0) {
      warning(sprintf("SigmaHatSq estimated to be %7.4f, setting to 0...",sigmahatsq))
      sigmahatsq <- 0.0
    }
    w <- 1/(omegasq+sigmahatsq)
    model.FGLS <- lm(y~X-1, weights=w)
    model.FGLS$sigmahatsq <- sigmahatsq
  }
  names(model.FGLS$coefficients) <- xnames

  #**NEW**#
  detach(dtx)

  model.FGLS
}


#
#  A few extra functions to run Monte Carlo experiments on edvreg (not needed for real applications)...
#
edvreg.sim <- function(N=1500,theta=0.9,rho=-0.9,bet0=0,bet1=1,Cval=0.9) {
  omegasq <- rgamma(N,1/theta,1/theta)
  u <- rnorm(N,mean=0,sd=sqrt(omegasq))
  eps <- rnorm(N,mean=0,sd=1)
  x <- matrix((sqrt(1-(rho^2)))*rnorm(N,mean=0,sd=1) + rho*omegasq/sqrt(theta))
  ystar <- bet0 + bet1*x + sqrt(Cval)*u + sqrt(1-Cval)*eps
  omegasq <- omegasq*Cval
  list(ystar=ystar,x=x,omegasq=omegasq)
}

edvreg.test <- function() {
  sims <- 1000
  alpha <- matrix(NA,sims,6)
  beta  <- alpha
  vareps <- rep(NA,sims)
  omegasqFactor <- vareps

  for (i in 1:sims) {
    cat(sprintf("Iteration %5i...\n",i))
    dat <- edvreg.sim()
    omegasq <- dat$omegasq
    ystar <- dat$ystar
    x <- dat$x
    res0 <- summary(lm(ystar ~ x))
    res1 <- summary( edvreg(ystar ~ x, omegasq=omegasq) )
    res2a <- edvreg(ystar ~ x, omegasq=omegasq,proportional=TRUE)
    res2 <- summary(res2a)
    rr <- cbind(res0$coefficients[,1:2],
                res1$coefficients[,1:2],
                res2$coefficients[,1:2])
    alpha[i,] <- rr[1,]
    beta[i,]  <- rr[2,]
    omegasqFactor[i] <- res2a$omegasqFactor
    vareps[i] <- res2a$vareps
  }
  list(alpha=alpha,beta=beta,vareps=vareps,omegasqFactor=omegasqFactor)
}

edvreg.run.test <- function() {
  res <- edvreg.test()
  print( apply(res$beta[,c(1,3,5)],2,sd) )
  print( apply(res$beta[,c(2,4,6)],2,mean) )
  print( mean(res$vareps) )
  print( sd( res$vareps) )
  print( mean(res$omegasqFactor) )
  print( sd(res$omegasqFactor) )
  return( res )
}
