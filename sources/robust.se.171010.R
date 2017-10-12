## Robust and Cluster-Robust SE Estimation ##
## From "https://thetarzan.wordpress.com/2011/05/28/heteroskedasticity-robust-and-clustered-standard-errors-in-r/"

## Heteroskedasticity-robust standard error calculation.
summary.robust <- function(model) {
  s <- summary(model)
  X <- model.matrix(model)
  u2 <- residuals(model)^2
  XDX <- 0
  
  ## Here one needs to calculate X'DX. But due to the fact that
  ## D is huge (NxN), it is better to do it with a cycle.
  for(i in 1:nrow(X)) {
    XDX <- XDX + u2[i]*X[i,]%*%t(X[i,])
  }
  
  # inverse(X'X)
  XX1 <- solve(t(X)%*%X)
  
  # Variance calculation (Bread x meat x Bread)
  varcovar <- XX1 %*% XDX %*% XX1
  
  # degrees of freedom adjustment
  dfc <- sqrt(nrow(X))/sqrt(nrow(X)-ncol(X))
  
  # Standard errors of the coefficient estimates are the
  # square roots of the diagonal elements
  stdh <- dfc*sqrt(diag(varcovar))
  
  t <- model$coefficients/stdh
  p <- 2*pnorm(-abs(t))
  results <- cbind(model$coefficients, stdh, t, p)
  dimnames(results) <- dimnames(s$coefficients)
  results
}

vcov.robust<-function(model) {
    s <- summary(model)
    X <- model.matrix(model)
    u2 <- residuals(model)^2
    XDX <- 0
    
    ## Here one needs to calculate X'DX. But due to the fact that
    ## D is huge (NxN), it is better to do it with a cycle.
    for(i in 1:nrow(X)) {
      XDX <- XDX + u2[i]*X[i,]%*%t(X[i,])
    }
    
    # inverse(X'X)
    XX1 <- solve(t(X)%*%X)
    
    # Variance calculation (Bread x meat x Bread)
    varcovar <- XX1 %*% XDX %*% XX1
    
    # degrees of freedom adjustment
    dfc <- (nrow(X))/(nrow(X)-ncol(X))
    
    # Standard errors of the coefficient estimates are the
    # square roots of the diagonal elements
    vcov <- dfc*(diag(varcovar))
    
    return(vcov)
  }
  
#To use the function written above, simply replace  summary()  
##with summary.robust() to look at your regression results ?\ like this:

#require(foreign)
#mrdr = read.dta(file="/Users/kevingoulding/data/MURDER.dta")

# run regression
#reg4 = lm(cmrdrte ~ cexec + cunem, data = subset(mrdr,year == 93))

# see heteroskedasticity-robust standard errors
#summary.robust(reg4)


## Cluster Robust Standard Errrors ##

cl.test <- function(dat,fm, cluster){
  attach(dat, warn.conflicts = F)
  require(lmtest)
  require(sandwich)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  detach(dat)
  coeftest(fm, vcovCL) 
}  

cl.vcov <- function(dat,fm, cluster){
  attach(dat, warn.conflicts = F)
  library(sandwich)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  detach(dat)
  return(vcovCL)
}  

cl.summary <- function(dat,fm, cluster){
  attach(dat, warn.conflicts = F)
  require(lmtest)
  require(sandwich)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  detach(dat)
  clsummary <- list(coeftest(fm, vcovCL),vcovCL)
  names(clsummary) <- c("summary.table","cluster.vcov")
  return(clsummary)
}

## Multiple Clusters. Taken from:
## http://stackoverflow.com/questions/8389843/double-clustered-standard-errors-for-panel-data

clm.test <- function(dat,fm, cluster1, cluster2){
  attach(dat, warn.conflicts = F)
  library(sandwich);library(lmtest)
  cluster12 = paste(cluster1,cluster2, sep="")
  M1  <- length(unique(cluster1))
  M2  <- length(unique(cluster2))   
  M12 <- length(unique(cluster12))
  N   <- length(cluster1)          
  K   <- fm$rank             
  dfc1  <- (M1/(M1-1))*((N-1)/(N-K))  
  dfc2  <- (M2/(M2-1))*((N-1)/(N-K))  
  dfc12 <- (M12/(M12-1))*((N-1)/(N-K))  
  u1j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster1,  sum)) 
  u2j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster2,  sum)) 
  u12j  <- apply(estfun(fm), 2, function(x) tapply(x, cluster12, sum)) 
  vc1   <-  dfc1*sandwich(fm, meat=crossprod(u1j)/N )
  vc2   <-  dfc2*sandwich(fm, meat=crossprod(u2j)/N )
  vc12  <- dfc12*sandwich(fm, meat=crossprod(u12j)/N)
  vcovMCL <- vc1 + vc2 - vc12
  detach(dat)
  coeftest(fm, vcovMCL)}

clm.vcov <- function(dat,fm, cluster1, cluster2){
  attach(dat, warn.conflicts = F)
  library(sandwich);library(lmtest)
  cluster12 = paste(cluster1,cluster2, sep="")
  M1  <- length(unique(cluster1))
  M2  <- length(unique(cluster2))   
  M12 <- length(unique(cluster12))
  N   <- length(cluster1)          
  K   <- fm$rank             
  dfc1  <- (M1/(M1-1))*((N-1)/(N-K))  
  dfc2  <- (M2/(M2-1))*((N-1)/(N-K))  
  dfc12 <- (M12/(M12-1))*((N-1)/(N-K))  
  u1j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster1,  sum)) 
  u2j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster2,  sum)) 
  u12j  <- apply(estfun(fm), 2, function(x) tapply(x, cluster12, sum)) 
  vc1   <-  dfc1*sandwich(fm, meat=crossprod(u1j)/N )
  vc2   <-  dfc2*sandwich(fm, meat=crossprod(u2j)/N )
  vc12  <- dfc12*sandwich(fm, meat=crossprod(u12j)/N)
  vcovMCL <- vc1 + vc2 - vc12
  detach(dat)
  return(vcovMCL)
  }