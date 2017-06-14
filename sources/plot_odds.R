#####################
## Plot Odds Ratio ##
#####################

## Modified from the code in:
## http://www.surefoss.org/dataanalysis/plotting-odds-ratios-aka-a-forrestplot-with-ggplot2/

## Need to import gktheme

plot_odds<-function(x, title = NULL, orderval="original", intercept=TRUE, direct=FALSE){

if (direct){
  tmp <- exp(x)
} else {
  tmp<-data.frame(cbind(exp(coef(x)), exp(confint(x))))
}

odds<-tmp[-1,]
if (intercept) {
  odds <- tmp
}
names(odds)<-c('OR', 'lower', 'upper')
odds$vars<-row.names(odds)
ticks<-c(seq(.1, 1, by =.1), seq(0, 10, by =1), seq(10, 100, by =10))
if (orderval=="asis"){
  plotstart = ggplot(odds, aes(y= OR, x = vars ))
} else if (orderval=="ORorder") {
  plotstart = ggplot(odds, aes(y= OR, x = reorder(vars, OR) ))
} else if (orderval=="original") {
  plotstart = ggplot(odds, aes(y= OR, x = reorder(vars, (length(vars)+1) - seq(1,length(vars),1)) ))
}
plotstart + geom_point() +
geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
labs(title = title, x = 'Variables', y = 'Odds Ratio') + gktheme
}
