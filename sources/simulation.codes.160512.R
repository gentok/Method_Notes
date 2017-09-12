## Simulation Program for lm model
simu.lm<-function(mod,predprof,cls=NULL,cirange=0.95,iterate.num=1000,seedval=2345){
  coef_lm<-mod$coefficients
  if(is.null(cls)) vcov_lm<-vcov(mod)
  else vcov_lm<-vcovHC(mod,cluster=cls,type="HC1")
  ndraws<-iterate.num; set.seed(seedval)
  betadraw_lm <- mvrnorm(ndraws, coef_lm, vcov_lm)
  predres<-matrix(NA,nrow=nrow(predprof),ncol=5)
  colnames(predres)<-c("Mean","Median","SE","lowCI","upCI")
  cidef<-qnorm(1-(1-cirange)/2)
  for(i in 1:nrow(predprof)){
    xfix<-as.vector(c(1,predprof[i,]))
    predstore<-betadraw_lm%*%xfix
    meanpred<-mean(predstore)
    medianpred<-median(predstore)
    sdpred<-sd(predstore)
    cipred<-c(meanpred-cidef*sdpred,meanpred+cidef*sdpred)
    predres[i,]<-c(meanpred,medianpred,sdpred,cipred)
  }
  predres<-as.data.frame(predres)
  return(predres)
}

## Simulation Program for glm logit model
simu.logit<-function(mod,predprof,cls=NULL,cirange=0.95,iterate.num=1000,seedval=2345){
  require(faraway)
  coef_logit<-mod$coefficients
  if(is.null(cls)) vcov_logit<-vcov(mod)
  else vcov_logit<-vcovHC(mod,cluster=cls,type="HC1")
  ndraws<-iterate.num; set.seed(seedval)
  betadraw_logit <- mvrnorm(ndraws, coef_logit, vcov_logit)
  predres<-matrix(NA,nrow=nrow(predprof),ncol=5)
  colnames(predres)<-c("Mean","Median","SE","lowCI","upCI")
  cidef<-qnorm(1-(1-cirange)/2)
  for(i in 1:nrow(predprof)){
    xfix<-as.vector(c(1,predprof[i,]))
    predstore<-ilogit(betadraw_logit%*%xfix)
    meanpred<-mean(predstore)
    medianpred<-median(predstore)
    sdpred<-sd(predstore)
    cipred<-c(meanpred-cidef*sdpred,meanpred+cidef*sdpred)
    predres[i,]<-c(meanpred,medianpred,sdpred,cipred)
  }
  predres<-as.data.frame(predres)
  return(predres)
}

# Simulation Program For glm Probit Model
simu.probit<-function(mod,predprof,cls=NULL,cirange=0.95,iterate.num=1000,seedval=2345){
  coef_probit<-mod$coefficients
  if(is.null(cls)) vcov_probit<-vcov(mod)
  else vcov_probit<-vcovHC(mod,cluster=cls,type="HC1")
  ndraws<-iterate.num; set.seed(seedval)
  betadraw_probit <- mvrnorm(ndraws, coef_probit, vcov_probit)
  predres<-matrix(NA,nrow=nrow(predprof),ncol=5)
  colnames(predres)<-c("Mean","Median","SE","lowCI","upCI")
  cidef<-qnorm(1-(1-cirange)/2)
  for(i in 1:nrow(predprof)){
    xfix<-as.vector(c(1,predprof[i,]))
    predstore<-pnorm(betadraw_probit%*%xfix)
    meanpred<-mean(predstore)
    medianpred<-median(predstore)
    sdpred<-sd(predstore)
    cipred<-c(meanpred-cidef*sdpred,meanpred+cidef*sdpred)
    predres[i,]<-c(meanpred,medianpred,sdpred,cipred)
  }
  predres<-as.data.frame(predres)
  return(predres)
}

## Simulation Program for glm poisson model
simu.poisson<-function(mod,predprof,cls=NULL,cirange=0.95,iterate.num=1000,seedval=2345){
  coef_poisson<-mod$coefficients
  if(is.null(cls)) vcov_poisson<-vcov(mod)
  else vcov_poisson<-vcovHC(mod,cluster=cls,type="HC1")
  ndraws<-iterate.num; set.seed(seedval)
  betadraw_poisson <- mvrnorm(ndraws, coef_poisson, vcov_poisson)
  predres<-matrix(NA,nrow=nrow(predprof),ncol=5)
  colnames(predres)<-c("Mean","Median","SE","lowCI","upCI")
  cidef<-qnorm(1-(1-cirange)/2)
  for(i in 1:nrow(predprof)){
    xfix<-as.vector(c(1,predprof[i,]))
    predstore<-exp(betadraw_poisson%*%xfix)
    meanpred<-mean(predstore)
    medianpred<-median(predstore)
    sdpred<-sd(predstore)
    cipred<-c(meanpred-cidef*sdpred,meanpred+cidef*sdpred)
    predres[i,]<-c(meanpred,medianpred,sdpred,cipred)
  }
  predres<-as.data.frame(predres)
  return(predres)
}

## Simulation Program for glm negbinom model
simu.negbinom<-function(mod,predprof,cls=NULL,cirange=0.95,iterate.num=1000,seedval=2345){
  coef_negbinom<-mod$coefficients
  if(is.null(cls)) vcov_negbinom<-vcov(mod)
  else vcov_negbinom<-vcovHC(mod,cluster=cls,type="HC1")
  ndraws<-iterate.num; set.seed(seedval)
  betadraw_negbinom <- mvrnorm(ndraws, coef_negbinom, vcov_negbinom)
  predres<-matrix(NA,nrow=nrow(predprof),ncol=5)
  colnames(predres)<-c("Mean","Median","SE","lowCI","upCI")
  cidef<-qnorm(1-(1-cirange)/2)
  for(i in 1:nrow(predprof)){
    xfix<-as.vector(c(1,predprof[i,]))
    predstore<-exp(betadraw_negbinom%*%xfix)
    meanpred<-mean(predstore)
    medianpred<-median(predstore)
    sdpred<-sd(predstore)
    cipred<-c(meanpred-cidef*sdpred,meanpred+cidef*sdpred)
    predres[i,]<-c(meanpred,medianpred,sdpred,cipred)
  }
  predres<-as.data.frame(predres)
  return(predres)
}
