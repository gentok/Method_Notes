#######################
## Plot Coefficients ##
#######################

## Modified from the code in:
## http://www.surefoss.org/dataanalysis/plotting-odds-ratios-aka-a-forrestplot-with-ggplot2/

## Need to import gktheme
gktheme <-
  theme(axis.text=element_text(size=10, colour="black"),
        axis.title.x=element_text(size=12,face="bold", vjust=-1.5),
        axis.title.y=element_text(size=12,face="bold", vjust=1.5),
        plot.title=element_text(size=12,face="bold", vjust=2,hjust=0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.background = element_rect(fill=NA, colour="black", size=0.5, linetype=1),
        legend.background = element_rect(fill=NA,colour=NA),
        legend.position = c(0.75,0.68),
        legend.key.width = unit(1.5, "cm"))

plot_coef<-function(x, title = NULL, orderval="original", intercept=TRUE, direct=FALSE, odds=FALSE,
                    custom.variable.names = NULL){
#' @param title plot title (string)
#' @param orderval Order of Coefficients in the plot. "original" (default), "coeforder" or "asis"
#' @param intercept Boulean. If TRUE (default), intercept included in the plot.
#' @param direct Boulean. If FALSE (default), coefficients imported from model result (where coef() and confint() is applicable.)
#' If TRUE, coefficients imported directly from coefficients table (rows=variables, columns=(coefficient,lower CI, upper CI))
#' @param odds Boulean. FALSE (default). If TRUE, the exponent of the coefficients will be exported.

## Import Coefficients
if (direct){
  if (odds){
    tmp <- exp(x)
  } else {
    tmp <- x
  }
} else {
  if (odds){
    tmp<-data.frame(cbind(exp(coef(x)), exp(confint(x))))
  } else {
    tmp<-data.frame(cbind(coef(x), confint(x)))
  }
}
coefs <- tmp
names(coefs)<-c('CF', 'lower', 'upper')

## 
if (is.null(custom.variable.names)==TRUE) {
  coefs$vars <- row.names(coefs)  
  } else {
  coefs$vars <- c("(Intercept)", custom.variable.names)
  }

## Include/Exclude Intercept
if (!intercept) {
  coefs<-coefs[-1,]
}

## Odds Ratio or Not
if (odds){
  ticks<-c(seq(.1, 1, by =.1), seq(0, 10, by =1), seq(10, 100, by =10))
}

## Start Plotting
if (orderval=="asis"){
  plotstart = ggplot(coefs, aes(y= CF, x = vars ))
} else if (orderval=="coeforder") {
  plotstart = ggplot(coefs, aes(y= CF, x = reorder(vars, CF) ))
} else if (orderval=="original") {
  plotstart = ggplot(coefs, aes(y= CF, x = reorder(vars, (length(vars)+1) - seq(1,length(vars),1)) ))
}

## Intermediate Plot
plotmid <- plotstart + geom_point(size=3) +
geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
coord_flip() + gktheme

## Final Plot
if (odds){
  plotfin <- plotmid + scale_y_log10(breaks=ticks, labels = ticks) +
  geom_hline(yintercept = 1, linetype=2) +
  labs(title = title, x = 'Variables', y = 'Odds Ratio')
} else {
  plotfin <- plotmid + geom_hline(yintercept = 0, linetype=2) +
  labs(title = title, x = 'Variables', y = 'Coefficient')
}

## Return the Plot
return(plotfin)

}
