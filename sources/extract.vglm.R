# Modifying extract function of texreg to export table
require(texreg)
extract.vglm <- function (model, 
                          include.aic = TRUE,
                          include.bic = TRUE,
                          include.loglik = TRUE, 
                          include.df = FALSE, 
                          include.nobs = TRUE,
                          beside = TRUE,
                          resp.names = NA,
                          ...) 
{
  s <- summary(model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    gof <- c(gof, AIC(model))
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, BIC(model))
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, VGAM::logLik.vlm(model))
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.df == TRUE) {
    gof <- c(gof, df <- s@df[2])
    gof.names <- c(gof.names, "DF")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, nobs(s))
    gof.names <- c(gof.names, "Num.\\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  besidereq <- nrow(s@coef3) > 
    length(s@extra$colnames.y) - 1 + length(all.vars(s@terms$terms)[-1])
  
  if (beside == TRUE & besidereq==TRUE) {
    trlist <- list()
    
    respcol <- s@extra$colnames.y
    if (is.na(resp.names)) resp.names <- respcol
    if (length(resp.names)!=length(respcol)) {
      warning("resp.names length does not match with number of response categories")
      resp.names <- respcol
    }
    
    for (i in 1:(length(respcol)-1)) {
      names <- rownames(coef(s))
      resploc <- grep(paste0(":",respcol[i],"$"),names)
      names <- gsub(paste0(":",respcol[i],"$"),"",names[resploc])
      co <- s@coef3[resploc, 1]
      se <- s@coef3[resploc, 2]
      pval <- s@coef3[resploc, 4]
      if (i==1) {
        tr <- createTexreg(coef.names = names, coef = co, se = se, 
                           pvalues = pval, gof.names = gof.names, 
                           gof = gof, gof.decimal = gof.decimal,
                           model.name = paste(resp.names[i],resp.names[i+1],sep="|"))
      } else {
        tr <- createTexreg(coef.names = names, coef = co, se = se, 
                           pvalues = pval, gof.names = character(), 
                           gof = numeric(), gof.decimal = logical(),
                           model.name = paste(resp.names[i],resp.names[i+1],sep="|"))
      }
      trlist <- c(trlist, tr)
    }
    if (length(trlist) == 1) {
      return(trlist[[1]])
    }
    else {
      return(trlist)
    }
  }
  else {
    names <- rownames(coef(s))
    co <- s@coef3[, 1]
    se <- s@coef3[, 2]
    pval <- s@coef3[, 4]
    tr <- createTexreg(coef.names = names, coef = co, se = se, 
                       pvalues = pval, gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
    return(tr)
  }
  
}
setMethod("extract", signature = className("vglm"), definition = extract.vglm)
