# script of GSA methods
# Jeffrey A. Walker
# August 22, 2016

library(data.table)

ols.fit <- function(dt,xcols,ycols,zcols,scale_it=TRUE,boot=TRUE,perms=2000,all=FALSE){
  # returns ols estimate and bootstrap SE
  # for parametric SE use Obrien?
  # if all=TRUE then return matrix with rows=perms
  n <- nrow(dt)
  m <- length(ycols)
  p <- length(xcols)
  rows <- 1:n
  b_mat <- matrix(NA,nrow=perms,ncol=length(zcols))
  colnames(b_mat) <- zcols
  for(iter in 1:perms){
    Y <- data.matrix(dt[rows,.SD,.SDcols=ycols])
    if(scale_it==TRUE){
      Y <- scale(Y)
      X1 <- data.matrix(dt[rows, .SD, .SDcols=setdiff(xcols,zcols)])
      X2 <- scale(data.matrix(dt[rows, .SD, .SDcols=zcols]))
      X.dm <- cbind(rep(1,n),X1,X2)
    }else{
      X.dm <- cbind(rep(1,n),data.matrix(dt[rows, .SD, .SDcols=xcols]))  # design matrix
      }
    fit <- lm.fit(X.dm,Y)
    b_mat[iter,] <- apply(fit$coefficients[zcols,],1,mean)
    rows <- sample(1:n,replace=TRUE)
  }
  b_table <- data.table(
    Ind.Var = zcols,
    Estimate = b_mat[1,],
    SE = apply(b_mat,2,sd)
    )
  b_table[,t:=Estimate/SE]
  b_table[,prob:=2*pt(abs(t),df=n-p-1,lower.tail = FALSE)] # conservative df, upper end should be n*m-p-1
  if(all==TRUE){
    return(b_mat)
  }else{
    return(b_table)
  }
}

ols.fit.long <- function(dt,xcols,ycols,zcols,scale_it=TRUE,boot=TRUE,perms=2000){
  # same as ols.fit but using long instead of wide format. Ouch!
  # microbenchmark results using default settings
  #Unit: seconds
  #expr      min       lq     mean   median       uq      max neval
  #ols.fit(dt, xcols, ycols, zcols) 11.50812 11.50812 11.50812 11.50812 11.50812 11.50812     1
  #ols.fit.long(dt, xcols, ycols, zcols) 97.90578 97.90578 97.90578 97.90578 97.90578 97.90578     1
  
  n <- nrow(dt)
  m <- length(ycols)
  p <- length(xcols)
  rows <- 1:n
  b_mat <- matrix(NA,nrow=perms,ncol=length(zcols))
  colnames(b_mat) <- zcols
  for(iter in 1:perms){
    dts <- dt[rows]
    if(scale_it==TRUE){
      Y <- scale(data.matrix(dt[rows, .SD, .SDcols=ycols]))
      X1 <- data.matrix(dt[rows, .SD, .SDcols=setdiff(xcols,zcols)])
      X2 <- scale(data.matrix(dt[rows, .SD, .SDcols=zcols]))
      dts <- data.table(X1,X2,Y)
    }
    dts[,subject:=factor(.I)]
    dtlong <- melt(dts,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
    dtlong[,gene:=factor(gene)]
    dtlong <- orderBy(~subject + gene, dtlong)
    form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
    Y <- dtlong[,expression]
    X.mm <- model.matrix(form,data=dtlong)
    fit <- lm.fit(X.mm,Y)
    b_mat[iter,] <- coefficients(fit)[zcols]
    rows <- sample(1:n,replace=TRUE)
  }
  b_table <- data.table(
    Ind.Var = zcols,
    Estimate = b_mat[1,],
    SD = apply(b_mat,2,sd)
  )
  b_table[,t:=Estimate/SD]
  b_table[,prob:=2*pt(abs(t),df=n-p-1,lower.tail = FALSE)] # conservative df, upper end should be n*m-p-1
  return(b_table)
}


ols_estimates <- function(dt,xcols,ycols,zcols){
  # estimates computed both long and wide format to show equivalence
  # if boot==TRUE then return bootstrap se
  
  # wide format (multivariate)
  Y <- as.matrix(dt[,.SD,.SDcols=ycols])
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fit.mv <- lm(form, data=dt)
  bhat <- apply(coefficients(fit.mv)[zcols,],1,mean)
  
  # long format with multiple outcomes (gene expression levels) stacked into single column
  dt[,subject:=factor(.I)]
  dtlong <- melt(dt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  dtlong <- orderBy(~subject + gene, dtlong)
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  fit.long <- lm(form, data=dtlong)
  bhat.long <- coefficients(fit.long)[zcols]
  
  return(bhat)
}

obrien.fit <- function(dt,xcols,ycols,zcols){
  # O'Brien's 1984 OLS test for continuous Z
  # fast version of obrien using lm.fit
  # requires that I compute standard errors, t, and prob
  # dt is the data.table with the X specified by xcols and Y specified by ycols
  # zcols is the subarray of xcols to test. currently works with only 1 zcol to test
  # this would be easy to fix by setting whole thing in loop and then rbinding results
  i <- 1 # only 1 zcol
  n <- nrow(dt)
  X.dm <- cbind(rep(1,n),data.matrix(dt[, .SD, .SDcols=xcols])) # design matrix
  Xred <- cbind(rep(1,n),data.matrix(dt[, .SD, .SDcols=setdiff(xcols,zcols[i])]))
  Y <- data.matrix(dt[,.SD,.SDcols=ycols])
  m <- length(ycols)
  df <- nrow(dt) - length(xcols) - 1
  
  # Get residuals from X (so excluding zcols) to find R - the correlation among the outcomes not explained by zcols
  fit <- lm.fit(Xred,Y)
  R <- cor(fit$residuals)
  
  # fit full model
  XTXI <- solve(t(X.dm)%*%X.dm)
  fit <- lm.fit(X.dm,Y)
  b <- fit$coefficients[zcols[i],]
  e <- fit$residuals
  se <- sqrt(diag((t(e)%*%e)/df)*XTXI[zcols[i],zcols[i]])
  t_value <- b/se
  #coef_table <- data.table(Estimate=b,se=se,t=t_value)
  t_sum <- sum(t_value)
  
  obrien.b <- mean(b)
  obrien.t <- t_sum/sqrt(sum(R))
  obrien.sd <- obrien.b/obrien.t
  obrien.p <- 2*pt(abs(obrien.t),df=df,lower.tail = FALSE)
  obrien_table <- data.table(zcols=zcols,b=obrien.b,sd=obrien.sd,t=obrien.t,p=obrien.p)
  return(obrien_table)
}

obrien <- function(dt,xcols,ycols,zcols){
  # O'Brien's 1984 OLS test for continuous Z
  # dt is the data.table with the X specified by xcols and Y specified by ycols
  # zcols is the subarray of xcols to test. currently works with only 1 zcol to test
  i <- 1 # only 1 zcol
  Y <- data.matrix(dt[,.SD,.SDcols=ycols])
  m <- length(ycols)
  n <- nrow(dt)
  df <- nrow(dt) - length(xcols) - 1
  
  # coefficients (save t-value in addition to coefficients)
  coef <- matrix(0,nrow=m,ncol=4)
  colnames(coef) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  
  non_zcols <- which(xcols!=zcols)
  # Get residuals from X (so excluding zcols)
  form <- formula(paste('Y',paste(xcols[non_zcols],collapse='+'),sep='~'))
  R <- cor(residuals(lm(form, data=dt)))
  
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv <- lm(form, data=dt)
  sum_fit <- summary(fitmv)
  row <- which(row.names(sum_fit[[paste('Response ',ycols[1],sep='')]]$coefficients)==zcols)
  for(j in 1:m){ # save t-values instead of coefficients
    coef[j,] <- sum_fit[[paste('Response ',ycols[j],sep='')]]$coefficients[row, ]
  }
  obrien.b <- mean(coef[,'Estimate'])
  obrien.t <- sum(coef[,'t value'])/sqrt(sum(R))
  obrien.p <- 2*pt(abs(obrien.t),df=df,lower.tail = FALSE)
  obrien_table <- data.table(zcols=zcols,b=obrien.b,t=obrien.t,p=obrien.p)
  return(obrien_table)
}

permutation_t.fit <- function(dt,xcols,ycols,zcols,method='R2',perms=2000, write_it=FALSE,fn){
  # no testing if there are no nuisance covariates
  # methods include different test statistics
    # t - mean t-value over the m responses
    # R2 - Anderson's R^2 squared partial correlation coefficient
    # maxmean - Efron and Tabrishini's maxmean statistic
  # faster fit using lm.fit
  # permutation using Anderson and Robinson 2001
  # dt is a data.table with the X regressors and Y responses
  # xcols are the regressors
  # ycols are the responses
  # zcols is the xcols to return statistics
  # method=resid uses Anderson permutation of residuals
  # method=pred uses GlobalAncova permutation of predictors
  # fn is the file name to write to

  m <- length(ycols)
  Y <- data.matrix(dt[, .SD, .SDcols=ycols])
  X <- cbind(rep(1,nrow(dt)),data.matrix(dt[, .SD, .SDcols=xcols]))
  XTXI <- solve(t(X)%*%X)
  df <- nrow(dt) - length(xcols) - 1
  
  res_table <- data.table(matrix(0.0,nrow=perms,ncol=length(zcols)))
  setnames(res_table,zcols)
  for(i in 1:length(zcols)){
    # get residuals from covariates
    if(length(xcols)==length(zcols)){ # no nuissance covariate
      covs <- xcols
    }else{
      covs <- setdiff(xcols, zcols[i])
    }
    X.red <- cbind(rep(1,nrow(dt)),data.matrix(dt[, .SD, .SDcols=covs]))
    fit.obs <- lm.fit(X.red,Y)
    e <- fit.obs$residuals
    yhat <- fit.obs$fitted.values
    
    Z <- data.matrix(dt[,.SD,.SDcols=zcols[i]])
    fit.obs <- lm.fit(X.red,Z)
    Rz.x <- fit.obs$residuals # residuals of Z on X
    Rz.x.sqr <- Rz.x^2
    rows <- 1:nrow(dt) # observed on first iter and permuted after
    for(iter in 1:perms){
      Y.pi <- yhat + e[rows,] # permuted
      #form <- formula(paste('Y.pi',paste(c(covs,zcols[i]),collapse='+'),sep='~'))
      #fitmv.pi <- lm(form, data=dt)
      #fitmv_sm <- summary(fitmv.pi)
      #fitmv_sm[[paste('Response ',ycols[j],sep='')]]$coefficients[zcols[i], ]
      
      if(method=='t'){
        fit.obs <- lm.fit(X,Y.pi)
        b.pi <- fit.obs$coefficients[zcols[i],]
        e.pi <- fit.obs$residuals
        se.pi <- sqrt(diag((t(e.pi)%*%e.pi)/df)*XTXI[zcols[i],zcols[i]])
        t_sum <- sum(b.pi/se.pi)
        res_table[iter,(zcols[i]):=t_sum]
      }
      
      # anderson and robertson partial correlation test statistic
      if(method=='R2'){
        fit.obs <- lm.fit(X.red,Y.pi)
        e.pi <- fit.obs$residuals
        num <- (sum(e.pi*Rz.x))^2
        denom <- sum(e.pi^2)*sum(Rz.x.sqr)
        R.sqr <- num/denom
        res_table[iter,(zcols[i]):=R.sqr]
      }
      
      # Efron and tabrishini maxmean statistic
      if(method=='maxmean'){
        fit.obs <- lm.fit(X,Y.pi)
        b.pi <- fit.obs$coefficients[zcols[i],]
        b.posbar <- 0
        b.negbar <- 0
        b.pos <- b.pi[b.pi>0]
        if(length(b.pos)>0){b.posbar <- mean(b.pos)}
        b.neg <- b.pi[b.pi<0]
        if(length(b.neg)>0){b.negbar <- abs(mean(b.neg))}
        maxmean <- max(b.posbar,b.negbar)
        res_table[iter,(zcols[i]):=maxmean]
      }
      
      # permute rows
      rows <- sample(1:nrow(dt))
    }
  }
  
  prob <- apply(res_table,2,function(x) length(which(abs(x) >= abs(x[1])))/perms)
  return(prob)
}

permutation_F.fit <- function(dt,xcols,ycols,zcols,method='resid',perms=2000, write_it=FALSE,fn){
  # permutation_F using lm.fit
  # GlobalAncova Fga statistic but permutation following Anderson and Robinson 2001
  # dt is a data.table with the X regressors and Y responses
  # xcols are the regressors
  # ycols are the responses
  # method is not implemented here but is in the original permutation_F.
  # method=resid uses Anderson permutation of residuals
  # method=pred uses GlobalAncova permutation of predictors
  # fn is the file name to write to
  # notes: the expected association between permuted hedonic score and gene expression is zero so the expected delta is zero E(E(b.h)-E(b.e))=0-0.

  p <- length(ycols)
  Y <- data.matrix(dt[, .SD, .SDcols=ycols])
  
  Fga <- data.table(matrix(0.0,nrow=perms,ncol=length(zcols)))
  setnames(Fga,zcols)
  for(i in 1:length(zcols)){
    # get residuals from covariates
    covs <- setdiff(xcols, zcols[i])
    #form <- formula(paste('Y',paste(covs,collapse='+'),sep='~'))
    #fit.obs <- lm(form, data=dt)
    #e <- residuals(fit.obs)
    #yhat <- predict(fit.obs) # Yhat = aX
    X.full <- cbind(rep(1,nrow(dt)),data.matrix(dt[, .SD, .SDcols=c(covs,zcols[i])]))
    X.red <- cbind(rep(1,nrow(dt)),data.matrix(dt[, .SD, .SDcols=covs]))
    fit.obs <- lm.fit(X.red,Y)
    e <- fit.obs$residuals
    yhat <- fit.obs$fitted.values
    # e = e1 and yhat=yhat1 to xxx decimal place
    
    rows <- 1:nrow(dt) # observed on first iter and permuted after
    for(iter in 1:perms){
      Y.pi <- yhat + e[rows,] # permuted
      # full
      fitmv.pi <- lm.fit(X.full,Y.pi) # full model
      e.pi <- fitmv.pi$residuals
      rss.full <- sum(e.pi^2)
      # reduced
      fitmv.pi <- lm.fit(X.red,Y.pi) # reduced model
      e.pi <- fitmv.pi$residuals
      rss.red <- sum(e.pi^2)
      Fga[iter,(zcols[i]):=(rss.red-rss.full)/rss.full]
      
      # permute rows
      rows <- sample(1:nrow(dt))
    }
  }
  
  prob <- apply(Fga,2,function(x) length(which(x>=x[1]))/perms)
  return(prob)
}
