################################################################################# 
## File Name: conjoint_functions.R                                             ##
## Date: 27 Sep 2017                                                           ##
## Author: Gento Kato                                                          ##
## Project: Innumeracy and Conjoint Policy Preference                          ##
## Purpose: Create Funtions to be used in the analysis                         ##
################################################################################# 

## Import Necessary Libraries

library(ggplot2)

#########################
## Mean Graph Function ##
#########################

meanplot <- 
function(plotdt,dv,policy,policyID,facet=F,facetvar=NULL,
         titletxt=NULL,boxplot=F){ #utility,policy,policyID,

plotdt$dv <- plotdt[,dv]
plotdt$policy <- plotdt[,policy]
plotdt$policyID <- plotdt[,policyID]
plotdt$facetvar <- plotdt[,facetvar]
  
gktheme0 <-
  theme(axis.text=element_text(size=10, colour="black"),
        axis.title.x=element_text(size=12,face="bold", vjust=-1.5),
        axis.title.y=element_text(size=12,face="bold", vjust=1.5),
        plot.title=element_text(size=12,face="bold", vjust=2,hjust=0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.background = element_rect(fill=NA, colour="black", size=0.5, linetype=1),
        legend.background = element_rect(fill=NA,colour=NA),
        legend.position = "none", #c(0.75,0.68),
        legend.key.width = unit(1.5, "cm"))

if (facet==F){
  plotstart = ggplot(plotdt, 
                     aes(y = dv,
                         x = reorder(policy, (length(policy)+1)-seq(1,length(policy),1)),
                         fill = as.factor(policyID)))
} else if (facet==T) {
  plotstart = ggplot(plotdt[!is.na(plotdt$facetvar),], 
                     aes(y = dv,
                         x = reorder(policy, (length(policy)+1)-seq(1,length(policy),1)),
                         fill = as.factor(policyID))) 
} 

if (boxplot==F){
  plotmed = plotstart + geom_bar(stat = "summary", fun.y = "mean")
} else if (boxplot==T) {
  plotmed = plotstart + geom_boxplot(lwd=0.5) # , outlier.shape=NA + scale_y_continuous(limits=c(0,4))
}
  
plotfin = plotmed + coord_flip() + gktheme0 + 
  labs(title = titletxt, x = 'Policies', y = 'Mean')

if (facet==T) {
  plotfin = plotfin + facet_grid(. ~ facetvar, scales="fixed")
}

return(plotfin)
}

#####################################
## Simulation Program for lm model ##
#####################################

simu.lm<-function(mod,predprof,cls=NULL,cirange=0.95,iterate.num=1000,seedval=2345){
  require(MASS)
  coef_lm<-mod$coefficients
  if(is.null(cls)) vcov_lm<-vcov(mod)
  else vcov_lm<-vcovHC(mod,cluster=cls,type="HC1")
  ndraws<-iterate.num; set.seed(seedval)
  betadraw_lm <- mvrnorm(ndraws, coef_lm, vcov_lm)
  predres<-matrix(NA,nrow=nrow(predprof),ncol=5)
  colnames(predres)<-c("Mean","Median","SE","lowCI","upCI")
  cidef <- c((1-cirange)/2,1-(1-cirange)/2)
  #cidef <- qnorm(1-(1-cirange)/2)
  for(i in 1:nrow(predprof)){
    xfix<-as.vector(c(1,as.matrix(predprof)[i,]))
    predstore<-betadraw_lm%*%xfix
    meanpred<-mean(predstore)
    medianpred<-median(predstore)
    sdpred<-sd(predstore)
    cipred <- quantile(predstore,probs=cidef)
    #cipred<-c(meanpred-cidef*sdpred,meanpred+cidef*sdpred)
    predres[i,]<-c(meanpred,medianpred,sdpred,cipred)
  }
  predres<-as.data.frame(predres)
  return(predres)
}

####################################
## Mean Plot for Predicted Result ##
####################################

meanplot_pred <- 
  function(plotdt,dv,lowCI,upCI,policy,policyID,corrected,facet=F,facetvar=NULL,titletxt=NULL){ 
    
    plotdt$dv <- plotdt[,dv]
    plotdt$upCI <- plotdt[,upCI]
    plotdt$lowCI <- plotdt[,lowCI]
    plotdt$policy <- plotdt[,policy]
    plotdt$policyID <- plotdt[,policyID]
    plotdt$facetvar <- plotdt[,facetvar]
    plotdt$corrected <- plotdt[,corrected]
    
    gktheme1 <-
      theme(axis.text=element_text(size=10, colour="black"),
            axis.title.x=element_text(size=12,face="bold", vjust=-1.5),
            axis.title.y=element_text(size=12,face="bold", vjust=1.5),
            plot.title=element_text(size=12,face="bold", vjust=2,hjust=0.5),
            plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
            panel.grid.major = element_line(colour = "grey90"),
            panel.grid.minor = element_line(colour = "grey90"),
            panel.background = element_rect(fill=NA, colour="black", size=0.5, linetype=1),
            legend.background = element_rect(fill=NA,colour=NA),
            legend.position = "none", #c(0.75,0.68),
            legend.key.width = unit(1.5, "cm"))
    
    if (facet==F){
      plotstart = ggplot(plotdt, 
                         aes(y = dv,
                             x = reorder(policy, (length(policy)+1)-seq(1,length(policy),1)),
                             fill = as.factor(policyID), alpha = as.factor(corrected), 
                             group = as.factor(corrected)))
    } else if (facet==T) {
      plotstart = ggplot(plotdt[!is.na(plotdt$facetvar),], 
                         aes(y = dv,
                             x = reorder(policy, (length(policy)+1)-seq(1,length(policy),1)),
                             fill = as.factor(policyID), alpha = as.factor(corrected), 
                             group = as.factor(corrected))) # 
    } 
    
    plotmed = plotstart + geom_bar(stat = "identity", position = "dodge") + 
      geom_errorbar(aes(ymin=lowCI,ymax=upCI),width=0.5, position = "dodge", 
                    colour="black", alpha = 1) + 
      scale_alpha_manual(values = c(1,0.5))

    plotfin = plotmed + coord_flip() + gktheme1 + 
      labs(title = titletxt, x = 'Policies', y = 'Mean')
    
    if (facet==T) {
      plotfin = plotfin + facet_grid(. ~ facetvar, scales="fixed")
    }
    
    return(plotfin)
  }

