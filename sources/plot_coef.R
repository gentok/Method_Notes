#####################
## Plot Odds Ratio ##
#####################

## Modified from the code in:
## http://www.surefoss.org/dataanalysis/plotting-odds-ratios-aka-a-forrestplot-with-ggplot2/

## Need to import gktheme

plot_coef<-function(x, title = NULL, orderval="original", intercept=TRUE, direct=FALSE, odds=FALSE){
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

## Include/Exclude Intercept
if (intercept) {
  coefs <- tmp
} else {
  coefs<-tmp[-1,]
}

names(coefs)<-c('CF', 'lower', 'upper')
coefs$vars<-row.names(coefs)

if (odds){
  ticks<-c(seq(.1, 1, by =.1), seq(0, 10, by =1), seq(10, 100, by =10))
}


## Start Plotting
if (orderval=="asis"){
  plotstart = ggplot(coefs, aes(y= CF, x = vars ))
} else if (orderval=="coeforder") {
  plotstart = ggplot(coefs, aes(y= CF, x = reorder(vars, CF) ))
} else if (orderval=="original") {
  plotstart = ggplot(coefs, aes(y= OR, x = reorder(vars, (length(vars)+1) - seq(1,length(vars),1)) ))
}

## Intermediate Plot
plotmid <- plotstart + geom_point() +
geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() + gktheme

## Final Plot
if (odds){
  plotfin <- plotmid + scale_y_log10(breaks=ticks, labels = ticks) +
  labs(title = title, x = 'Variables', y = 'Odds Ratio')
} else {
  plotfin <- plotmid + labs(title = title, x = 'Variables', y = 'Coefficient')
}

## Return the Plot
return(plotfin)

}
