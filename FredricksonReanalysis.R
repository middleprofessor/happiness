# script to re-analyze Fredrickson et al 2013 and 2015
# Jeffrey A. Walker
# August 22, 2016
# clean version of Fredrickson_multivariate_lm.peerJ.rev1.R
# and an exact replica of Fredrickson_peerJ.rev2.working.R
# to post to github
# the scripts used for this manuscript are
# Fredrickson_peerj.rev2.R - scripts mostly specific to these datasets
# GSAmethods.R - generalized methods for any dataset
# GSA1.R - the Monte Carlo simulation of type I, II, S, M errors

library(data.table)
library(ggplot2)
library(reshape2)
library(GlobalAncova) #bioconductor 
library(globaltest) #bioconductor 
library(limma) #bioconductor
library(mvtnorm)
library(nlme)
library(geepack) #ditto
library(doBy)
library(showtext) # needed for eps fonts to add Arial to .eps
font.add('Arial',regular='Arial.ttf')
library(gridExtra)
library("grid") # needed for "unit" function in ggplot

# to install bioconductor packages use
# setRepositories()
# and choose "1 2" which is CRAN + BioC software

# for Monte Carlo simulation of error rates run function start_here() in the script GSA1.R

do_Fredrickson <- function(){
  run_ols <- FALSE
  run_gee <- FALSE
  run_gls <- TRUE
  run_cole15 <- FALSE
  run_gls_parametric <- FALSE
  run_gls_bootstrap <- FALSE
  run_gls_permutation <- FALSE

  dt2013 <- get_cole_data(fn='cole1_clean.txt',year=2013,scale_it=TRUE)
  dt2015 <- get_cole_data(fn='cole2_clean.txt',year=2015,scale_it=TRUE)
  dt2015B <- get_cole_2015B(scale_it=TRUE)
  # combine the data
  # illness in FRED2015 is averaged over the 13 categories so divide illness in FRED13 by 13 to be in same scale as in FRED15. Even with this the range of FRED15 is about 2X that of FRED13.
  #dt2013[,illness:=illness/13] # comment out to replicate FRED15
  # remove IL6 from dt2013
  redcols <- c(get_xcols(),c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes()))
  dtCombi <- rbind(data.table(dt2013[,.SD,.SDcols=redcols],study=2013),data.table(dt2015,study=2015)) # use Fill=TRUE and full dt2015 data to retain IL6 to replicate FRED15
  # dtCombi <- rbind(data.table(dt2013,study=2013),data.table(dt2015,study=2015),fill=TRUE) # use Fill=TRUE and full dt2015 data to retain IL6 to replicate FRED15
  dtCombi[,study:=factor(study)]
  
  ycols13 <- c(pro_inflam_genes(year=2013),antibody_genes(),ifn_genes())
  ycols15 <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  xcols <- get_xcols() # for all FRED datasets
  zcols <- c('zhedonia','zeudaimonia')
  
  #do GLS and save so I don't have to keep doing this
  if(run_gls_parametric==TRUE){
    saveRDS(gls_with_correlated_error(dt2013, xcols=xcols, ycols=ycols13, zcols=zcols, method='gls'), "FRED13.gls.rds")
    saveRDS(gls_with_correlated_error(dt2015, xcols=xcols, ycols=ycols15, zcols=zcols, method='gls'), "FRED15.gls.rds")
    saveRDS(gls_with_correlated_error(dtCombi, xcols=c(xcols,'study'), ycols=ycols15, zcols=zcols, method='gls'), "FRED.Combi.gls.rds")
    # to replicate FRED15
    # saveRDS(gls_with_correlated_error(dtCombi, xcols=c(xcols,'study'), ycols=ycols13, zcols=zcols, method='gls'), "FRED.Combi-rep.gls.rds")
    
    # redo 2013 and 2015 datasets without zhedonia to compare to COLE15
    xcolsb <- get_xcolsb() # for COLE15
    xcols1 <- setdiff(xcols,'zhedonia')
    saveRDS(gls_with_correlated_error(dt2013, xcols=xcols1, ycols=ycols13,zcols='zeudaimonia', method='gls'), "FRED13.no-hedonia.gls.rds")
    saveRDS(gls_with_correlated_error(dt2015, xcols=xcols1, ycols=ycols15,zcols='zeudaimonia', method='gls'), "FRED15.no-hedonia.gls.rds")
    # COLE15 dropping 'hispanic' and 'ln_hh_income' # see original for more variations
    saveRDS(gls_with_correlated_error(dt2015B[,.SD,.SDcols=setdiff(colnames(dt2015B),c('ln_hh_income','hispanic'))], xcols=setdiff(xcolsb,c('ln_hh_income','hispanic')), ycols=ycols13,zcols='zeudaimonia', method='gls'), "COLE15.gls.rds")
    
    # redo all three without smoke to see effect on standard error
    xcolsc <- setdiff(xcols,'smoke')
    saveRDS(gls_with_correlated_error(dt2013[,.SD,.SDcols=c(xcolsc,ycols13)], xcols=xcolsc, ycols=ycols13, zcols=zcols, method='gls'), "FRED13-nosmoke.gls.rds")
    saveRDS(gls_with_correlated_error(dt2015[,.SD,.SDcols=c(xcolsc,ycols15)], xcols=xcolsc, ycols=ycols15, zcols=zcols, method='gls'), "FRED15-nosmoke.gls.rds")
    saveRDS(gls_with_correlated_error(dtCombi[,.SD,.SDcols=c(xcolsc,'study',ycols15)], xcols=c(xcolsc,'study'), ycols=ycols15, zcols=zcols, method='gls'), "FRED.Combi-nosmoke.gls.rds")
    
  }
  
  if(run_gls_bootstrap == TRUE){
    xcolsc <- setdiff(xcols,'smoke')
    bootstrap_models(dt2013[,.SD,.SDcols=c(xcolsc,ycols13)], xcols=xcolsc, ycols=ycols13, zcols=zcols,which_file='FRED13',tests=c('gls'), niter=2)
    bootstrap_models(dt2015[,.SD,.SDcols=c(xcolsc,ycols15)], xcols=xcolsc, ycols=ycols15, zcols=zcols,which_file='FRED15',tests=c('gls'), niter=201)
    bootstrap_models(dtCombi[,.SD,.SDcols=c(xcolsc,'study',ycols15)], xcols=c(xcolsc,'study'), ycols=ycols15, zcols=zcols,which_file='FRED.Combi',tests=c('gls'), niter=201)    
  }
  
  if(run_gls_permutation==TRUE){
    permutation_gls(dt=dt2013,xcols=xcols,ycols=ycols13,zcols=zcols,method='gls',perms=200, write_it=TRUE,fn='perm_gls.FRED13')
    permutation_gls(dt=dt2015,xcols=xcols,ycols=ycols15,zcols=zcols,method='gls',perms=200, write_it=TRUE,fn='perm_gls.FRED15')
    permutation_gls(dt=dtCombi,xcols=c(xcols,'study'),ycols=ycols15,zcols=zcols,method='gls',perms=200, write_it=TRUE,fn='perm_gls.FRED.Combi')
  }
  
  which_file_list <- c('FRED13', 'FRED15','FRED.Combi')
  gls_table <- NULL
  ols_table <- NULL
  gee_table <- NULL
  gls_supp_table <- NULL # supplement
  
  for(which_file in which_file_list){
    covs <- xcols
    if(which_file=='FRED13'){dt <- copy(dt2013)}
    if(which_file=='FRED15'){dt <- copy(dt2015)}
    if(which_file=='FRED.Combi'){
      dt <- copy(dtCombi)
      covs <- c(xcols,'study')
      }
    if('IL6' %in% colnames(dt)){ycols <- ycols13}else{ycols <- ycols15}
    
    # some statistics on the gene expression levels
    R <- cor(as.matrix(dt[,.SD,.SDcols=ycols]))
    mean(abs(R[lower.tri(R)]))
    max(abs(R[lower.tri(R)]))
    
    # OLS table
    if(run_ols==TRUE){
      boot_mat <- ols.fit(dt,xcols=covs,ycols,zcols,perms=10000,scale_it=TRUE,boot=TRUE,all=TRUE)
      boot_t <- make_ols_table(boot_mat)
      # obrien_t <- rbind(obrien.fit(dt,xcols=covs,ycols,zcols='zhedonia'),obrien.fit(dt,xcols=covs,ycols,zcols='zeudaimonia'))
      # replace obrien_t with table computed including delta
      obrien_t <- obrien.fit.delta(dt,xcols=covs,ycols,zcols)
      perm_R2 <- permutation_t.fit(dt,xcols=covs,ycols,zcols,method='R2',perms=10000, write_it=FALSE,fn)
      perm_f <- permutation_F.fit(dt,xcols=covs,ycols,zcols,perms=10000)
      ga <- run_GlobalAncova(dt,xcols=covs,ycols,zcols,perms=10000)
      rot_z <- run_Roast(dt,xcols=covs,ycols,zcols=zcols,perms=10000)
      
      # create NA for delta entry
      perm_R2 <- c(perm_R2,NA)
      perm_f <- c(perm_f,NA)
      ga <- c(ga,NA)

      # note that this assumes p-values are in same order as there is no checking of row.names
      ols_table <- rbind(ols_table,data.table(
        Data=which_file,
        boot_t[,.(Type=Type,Estimate,SE_boot=SE)],
        SE_obrien=obrien_t[,SE],
        obrien=obrien_t[,prob],
        permR2=perm_R2,
        ga=ga,
        permf=perm_f,
        roast=rot_z
      ))
    }
    
    if(run_gee==TRUE){
      # gee table
      fit <- gls_with_correlated_error(dt,xcols=covs,ycols,zcols,method='gee')
      gee_res <- make_gls_table(fit,method='gee')
      gee_part <- data.table(Type=row.names(gee_res),Data=which_file,gee_res[,c('Estimate', 'Std.err', 'Pr(>|W|)')])
      gee_part <- setNames(gee_part,c('Type','Data','Estimate','SE','p'))
      gee_table <- rbind(gee_table,gee_part)
    }
    
    if(run_gls==TRUE){
      # gls table
      # gls coefficients and p-values
      fit <- readRDS(paste(which_file,'.gls.rds',sep=''))
      gls_res <- make_gls_table(fit,method='gls')
      gls_part <- data.table(Type=row.names(gls_res),Data=which_file,gls_res[,c('Value', 'Std.Error', 'p-value')])
      gls_part <- setNames(gls_part,c('Type','Data','Estimate','SE','p'))
      
      #compute bootstrap.gls stats
      fn <- paste(which_file, '.bootstrap.gls.list.v2.txt',sep='')
      gls_boot_res <- read_bootstrap.gls_list(fn)
      nrow(gls_boot_res)
      gls_boot_res[,delta:=b.zhedonia-b.zeudaimonia]
      boot.se <- apply(gls_boot_res[,.SD,.SDcols=c('b.zhedonia','b.zeudaimonia','delta')],2,sd)
      gls_part <- cbind(gls_part,SE_boot=boot.se)
      
      #compute permutation.gls stats
      fn <- paste(which_file, '.permutation.gls.list.v1.txt',sep='')
      gls_perm_res.v1 <- read_permutation.gls_list(fn)
      fn <- paste(which_file, '.permutation.gls.list.v2.txt',sep='')
      gls_perm_res.v2 <- read_permutation.gls_list(fn)
      gls_perm_res.v2 <- convert_gls_permutation_to_old_format(gls_perm_res.v2)
      gls_perm_res <- gls_perm_res.v1 # v2 only to confirm new code
      nrow(gls_perm_res)
      gls_perm_p <- permutation_gls_p_value(gls_perm_res, statistic='t')
      gls_part <- cbind(gls_part,perm_p=gls_perm_p)
      gls_table <- rbind(gls_table,gls_part)
      
      #supplement (boot) table
      fit <- readRDS(paste(which_file,'-nosmoke.gls.rds',sep=''))
      gls_res_no_smoke <- data.table(make_gls_table(fit,method='gls'))
      gls_no_smoke_part <- gls_part[,.(Type,Data,b=Estimate,SE)]
      gls_no_smoke_part <- cbind(gls_no_smoke_part,gls_res_no_smoke[,.(b.nosmoke=Value,SE.nosmoke=Std.Error)])
      gls_no_smoke_part <- cbind(gls_no_smoke_part,gls_part[,.(SE_boot)])
      gls_supp_table <- rbind(gls_supp_table,gls_no_smoke_part)
      
     }
  } # end which file
  

  # clean tables
  
  if(run_ols==TRUE){
    ols_table_full <- copy(ols_table)
    write.table(ols_table_full, 'ols_table_full.txt',quote=FALSE,row.names=FALSE,sep='\t')
    ols_table <- data.table(read.table('ols_table_full.txt', header=TRUE, sep='\t'))
    # make Type first column and drop SE_boot, which as the smoke problem
    ols_table <- ols_table[,.(Type,Data,Estimate,SE_obrien,obrien,permR2,ga,permf,roast)]
    ols_table[, Type:=factor(Type)]
    ols_table <- orderBy(~-Type,ols_table)
    ols_table[,Estimate:=round(Estimate,3)]
    #ols_table[,SE_boot:=round(SE_boot,3)]
    ols_table[,SE_obrien:=round(SE_obrien,3)]
    ols_table[,obrien:=round(obrien,2)]
    ols_table[,permR2:=round(permR2,2)]
    ols_table[,ga:=round(ga,2)]
    ols_table[,permf:=round(permf,2)]
    ols_table[,roast:=round(roast,2)]
    # ols_table <- ols_table[Type!='delta'] # delta now in its own table
    write.table(ols_table, 'ols_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
  }
  
  if(run_gls==TRUE){
    gls_full_table <- copy(gls_table)
    gls_table[, Type:=factor(Type)]
    gls_table <- orderBy(~-Type, gls_table)
    gls_table[,Estimate:=round(Estimate,3)]
    gls_table[,SE:=round(SE,3)]
    gls_table[,p:=round(p,3)]
    gls_table[,SE_boot:=round(SE_boot,3)]
    gls_table[,perm_p:=round(perm_p,2)]
    write.table(gls_table, 'gls_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
    
    #supplemental table
    gls_supp_table_full <- copy(gls_supp_table)
    gls_supp_table[, Type:=factor(Type)]
    gls_supp_table <- orderBy(~-Type,gls_supp_table)
    gls_supp_table[,b:=round(b,3)]
    gls_supp_table[,SE:=round(SE,3)]
    gls_supp_table[,b.nosmoke:=round(b.nosmoke,3)]
    gls_supp_table[,SE.nosmoke:=round(SE.nosmoke,3)]
    gls_supp_table[,SE_boot:=round(SE_boot,3)]
    write.table(gls_supp_table, 'gls_supp_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
    
  }
 
  if(run_gee==TRUE){
    gee_table_full <- copy(gee_table)
    gee_table[, Type:=factor(Type)]
    gee_table <- orderBy(~-Type,gee_table)
    gee_table[,Estimate:=round(Estimate,3)]
    gee_table[,SE:=round(SE,3)]
    gee_table[,p:=round(p,2)]
    write.table(gee_table, 'gee_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
  }
  
  # make Cole15 table (table 2 in manuscript)
  if(run_cole15==TRUE){
    fit <- readRDS(paste('FRED13.no-hedonia.gls.rds'))
    cole15_table <- data.table(Data='FRED13',t(summary(fit)$tTable['zeudaimonia',]))
    fit <- readRDS(paste('FRED15.no-hedonia.gls.rds'))
    cole15_table <- rbind(cole15_table,data.table(Data='FRED15',t(summary(fit)$tTable['zeudaimonia',])))
    fit <- readRDS(paste('Cole15.gls.rds'))
    cole15_table <- rbind(cole15_table,data.table(Data='COLE15',t(summary(fit)$tTable['zeudaimonia',])))
    setnames(cole15_table,old=colnames(cole15_table),new=c('Data','b.eudaimonia','SE','t','p'))
    cole15_table <- cole15_table[,.(Data,b.eudaimonia,SE,p)]
    cole15_table[,b.eudaimonia:=round(b.eudaimonia,3)]
    cole15_table[,SE:=round(SE,3)]
    cole15_table[,p:=round(p,3)]
    write.table(cole15_table, 'cole15_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
  }
  
}

get_cole_data <- function(fn,year,scale_it=TRUE){
  dt <- read_file(fn,year=year)
  dt[,zhedonia:=scale(zhedonia)]
  dt[,zeudaimonia:=scale(zeudaimonia)]
  dt <- contrast_coefficients(dt) # convert to CTRA response
  if(scale_it==TRUE){
    xcols <- get_xcols()
    ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
    X <- dt[,.SD,.SDcols=xcols]
    Y <- scale(dt[,.SD,.SDcols=ycols])
    dt <- cbind(X,Y)
  }
  return(dt)
}

get_cole_2015B <- function(scale_it=TRUE){
  fn <- 'cole3_clean.txt'
  year <- 2013 #IL6 is present
  dt2015B <- data.table(read.table(fn,header=TRUE,sep='\t'))
  dt2015B[, female:=factor(female)]
  dt2015B[, black:=factor(black)]
  dt2015B[, smoke:=factor(smoke)]
  dt2015B[, hispanic:=factor(hispanic)]
  dt2015B[, alcohol:=factor(alcohol)]
  xcolsb <- get_xcolsb()
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  dt2015B <- na.omit(dt2015B[,.SD,.SDcols=c(xcolsb,ycols)])
  dt2015B[, zeudaimonia:=scale(zeudaimonia)]
  # note there is no zhedonia
  #scale
  if(scale_it==TRUE){
    X <- dt2015B[,.SD,.SDcols=xcolsb]
    Y <- scale(dt2015B[,.SD,.SDcols=ycols])
    dt2015B <- cbind(X,Y)
  }
  return(dt2015B)
}

read_file <- function(fn, year=2013){
  dt <- data.table(read.table(fn,header=TRUE,sep='\t'))
  dt[, male:=factor(male)]
  dt[, white:=factor(white)]
  dt[, smoke:=factor(smoke)]
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  dt <- na.omit(dt[,.SD,.SDcols=c(xcols,ycols)])
  return(dt)
}

contrast_coefficients <- function(dt){
  # dt is a matrix niter * p matrix of beta coefficients or raw data
  # compute contrast coefficients AND compute mean of these
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  rev_ycols <- c(antibody_genes(),ifn_genes())
  dt[,(rev_ycols):=lapply(.SD,"*",-1),.SDcols=rev_ycols]
  return(dt)
}

get_xcols <- function(){
  # These are the regressors
  xcols <- c('male' , 'age' , 'white' , 'bmi' , 'alcohol', 'smoke' , 'illness' , 'cd3d' , 'cd3e' , 'cd4', 'cd8a' , 'fcgr3a', 'cd19' , 'ncam1' , 'cd14' , 'zhedonia' , 'zeudaimonia')
  return(xcols)
}

get_xcolsb <- function(){
  xcolsb <- c('age','female','black','smoke','hispanic','bmi','diabcvdcastr','ln_hh_income','alcohol','CD3D','CD3E','CD4','CD8A','FCGR3A','CD19','NCAM1','CD14','zeudaimonia')
  return(xcolsb)
}

# ycol functions

pro_inflam_genes <- function(year=2013){
  # from Frederickson 2013, note that 2015 does not include IL6
  # 19 proinflammatory genes, which are up-regulated on average in the CTRA
  # if year=2013 include IL6, otherwise exclude it
  pro_inflam <- c('IL1A', 'IL1B', 'IL6', 'IL8', 'TNF', 'PTGS1', 'PTGS2', 'FOS', 'FOSB', 'FOSL1', 'FOSL2', 'JUN', 'JUNB', 'JUND', 'NFKB1', 'NFKB2', 'REL', 'RELA', 'RELB')
  if(year==2015){pro_inflam <- setdiff(pro_inflam, 'IL6')}
  return(pro_inflam)
}

antibody_genes <- function(){
  # from Frederickson 2013, note that 2015 does not include IL6
  # three genes involved in antibody synthesis, which are down-regulated on average in the CTRA
  antibody <- c('IGJ', 'IGLL1', 'IGLL3')
  return(antibody)
}

ifn_genes <- function(){
  # from Frederickson 2013, note that 2015 does not include IL6
  # 31 genes involved in type I IFN responses, which are down-regulated on average in the CTRA
  ifn <- c('GBP1', 'IFI16', 'IFI27', 'IFI27L1', 'IFI27L2', 'IFI30', 'IFI35', 'IFI44', 'IFI44L', 'IFI6', 'IFIH1', 'IFIT1', 'IFIT2', 'IFIT3', 'IFIT5', 'IFIT1L', 'IFITM1', 'IFITM2', 'IFITM3', 'IFITM4P', 'IFITM5', 'IFNB1', 'IRF2', 'IRF7', 'IRF8', 'MX1', 'MX2', 'OAS1', 'OAS2', 'OAS3', 'OASL')
  return(ifn)
}

obrien.fit.delta <- function(dt,xcols,ycols,zcols){
  # This is obrien.fit from the GSA_methods.R scripts but I've added the computation of the SE and p-value for the difference in effect between two of the zcols (Hedonia and Eudaimonia)
  n <- nrow(dt)
  Y <- data.matrix(dt[,.SD,.SDcols=ycols])
  m <- length(ycols)
  p <- length(xcols)
  df <- n - p - 1
  b <- matrix(NA,nrow=2,ncol=m)
  se <- matrix(NA,nrow=2,ncol=m)
  t_value <- matrix(NA,nrow=2,ncol=m)
  sumR <- numeric(2)
  
  X.dm <- cbind(rep(1,n),data.matrix(dt[, .SD, .SDcols=xcols])) # design matrix
  XTXI <- solve(t(X.dm)%*%X.dm)
  fit <- lm.fit(X.dm,Y)
  b <- fit$coefficients[zcols,]
  e <- fit$residuals
  for(i in 1:length(zcols)){
    Xred <- cbind(rep(1,n),data.matrix(dt[, .SD, .SDcols=setdiff(xcols,zcols[i])]))
    # Get residuals from X (so excluding zcols) to find R - the correlation among the outcomes not explained by zcols
    fit <- lm.fit(Xred,Y)
    R <- cor(fit$residuals)
    sumR[i] <- sum(R)
    se[i,] <- sqrt(diag((t(e)%*%e)/df)*XTXI[zcols[i],zcols[i]])
    t_value[i,] <- b[i,]/se[i,]
  }
  
  #coef_table <- data.table(Estimate=b,se=se,t=t_value)
  obrien.b <- apply(b,1,mean)
  obrien.t <- apply(t_value,1,sum)/sqrt(sumR)
  obrien.sd <- obrien.b/obrien.t
  obrien.p <- 2*pt(abs(obrien.t),df=df,lower.tail = FALSE)
  
  delta.b <- b[1,] - b[2,]
  delta.sd <- sqrt(se[1,]^2 + se[2,]^2)
  delta.t <- delta.b/delta.sd
  obrien.t.delta <- sum(delta.t)/sqrt(sumR[1] + sumR[2])
  p.delta <- 2*pt(abs(obrien.t.delta),df=df,lower.tail = FALSE)
  obrien_table <- data.table(Type=c(zcols,'delta'),Estimate=c(obrien.b,mean(delta.b)),SE=c(obrien.sd,mean(delta.b)/obrien.t.delta),t=c(obrien.t,obrien.t.delta),prob=c(obrien.p,p.delta))
  return(obrien_table)
}



gls_with_correlated_error <- function(dt,xcols,ycols,zcols, method='gls'){
  # generalized from original function to allow analysis of COLE15
  # dt is the data in wide format
  # xcols are the predictors
  # ycols are the resposes
  # zcols are the focal predictors to return statistics
  dt[,subject:=factor(.I)]
  dtlong <- melt(dt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  dtlong <- orderBy(~subject + gene, dtlong)
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  if(method=='gls'){
    fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE),na.action=na.omit)
  }
  if(method=='lme'){
    fit1 <- lme(form, random = ~1|subject, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(maxIter=100, msMaxIter = 500, tolerance=1e-6, msVerbose = FALSE)) # tolerance=1e-6 default
  }
  if(method=='gee'){
    fit1 <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='exchangeable', std.err="san.se")
  }
  return(fit1)
}



bootstrap_models <- function(dt, xcols, ycols, zcols, which_file, tests=c('mv','gls','gee'), niter=200, write_it=TRUE){
  # bootstrap estimates using multivariate regression, glm, and gee models
  covs <- copy(xcols)
  if('study' %in% colnames(dt)){
    combi <- TRUE
  }else{
      combi <- FALSE
  }
  Y <- data.matrix(dt[, .SD, .SDcols=ycols])
  
  rows <- 1:nrow(dt) # use if combi==FALSE
  # use if combi ==TRUE
  if(combi==TRUE){
    rows1 <- which(dt[,study]==levels(dt[,study])[1])
    rows2 <- which(dt[,study]==levels(dt[,study])[2])
    s_rows1 <- copy(rows1)
    s_rows2 <- copy(rows2)
  }
  mv_matrix <- matrix(0,nrow=niter,ncol=2)
  colnames(mv_matrix) <- zcols
  gee_table <- data.table(NULL)
  gls_table <- data.table(NULL)
  code <- paste(sample(LETTERS,4,replace=TRUE),collapse='')
  gee.out <- paste(which_file,code,'bootstrap.gee.table','txt',sep='.')
  gls.out <- paste(which_file,code,'bootstrap.gls.table','txt',sep='.')
  samp <- 'obs'
  for(iter in 1:niter){
    if(combi==FALSE){
      X <- dt[rows, .SD, .SDcols=covs]
      Y <- scale(as.matrix(dt[rows, .SD, .SDcols=ycols]))
      dts <- cbind(X,Y)
      dts[,zhedonia:=scale(zhedonia)]
      dts[,zeudaimonia:=scale(zeudaimonia)]
    }else{
      dt1 <- dt[s_rows1,]
      dt1[,zhedonia:=scale(zhedonia)]
      dt1[,zeudaimonia:=scale(zeudaimonia)]
      X1 <- dt1[, .SD, .SDcols=covs]
      Y1 <- scale(as.matrix(dt1[, .SD, .SDcols=ycols]))
      dt1 <- cbind(X1,Y1)
      
      dt2 <- dt[s_rows2,]
      dt2[,zhedonia:=scale(zhedonia)]
      dt2[,zeudaimonia:=scale(zeudaimonia)]
      X2 <- dt2[, .SD, .SDcols=covs]
      Y2 <- scale(as.matrix(dt2[, .SD, .SDcols=ycols]))
      dt2 <- cbind(X2,Y2)
      dts <- rbind(dt1,dt2)
    }
    
    if('mv' %in% tests){ # not implemented in update
      Y.samp <- as.matrix(dts[,.SD,.SDcols=ycols])
      form <- formula(paste('Y.samp~',paste(xcols,collapse='+'),sep=''))
      fit <- lm(form,data=dt[rows,])
      mv_matrix[iter,] <- apply(coefficients(fit)[zcols,], 1, mean)
    }
    if('gls' %in% tests){
      fit.gls <- gls_with_correlated_error(dts,xcols=covs,ycols=ycols,zcols=zcols,method='gls')
      estimate <- summary(fit.gls)$tTable[zcols,'Value']
      tvalue <- summary(fit.gls)$tTable[zcols,'t-value']
      gls_table <- rbind(gls_table,data.table(samp=samp,b=t(estimate),t=t(tvalue)))
      if(write_it==TRUE){write.table(gls_table,gls.out,quote=FALSE,row.names=FALSE,sep='\t')}
    }
    if('gee' %in% tests){
      fit.geeglm <- gls_with_correlated_error(dts,xcols=covs,ycols=ycols,zcols=zcols,method='gee')
      # summary(lm(form,data=dtlong))$coefficients[zcols,]
      # summary(fit.geeglm)$coefficients[zcols,]
      estimate <- summary(fit.geeglm)$coefficients[zcols,'Estimate']
      names(estimate) <- zcols
      gee_table <- rbind(gee_table,data.table(samp=samp,t(estimate)))
      if(write_it==TRUE){write.table(gee_table,gee.out,quote=FALSE,row.names=FALSE,sep='\t')}
    }
    rows <- sample(1:nrow(dt),replace=TRUE)
    if(combi==TRUE){
      s_rows1 <- sample(rows1,replace=TRUE)
      s_rows2 <- sample(rows2,replace=TRUE)
    }
    samp <- 'resample'
  }
  return(NULL)
}

permutation_gls <- function(dt,xcols,ycols,zcols,method='gls',perms=200, write_it=FALSE,fn){
  # uses Anderson permutation but fits with GLS
  # dt is a data.table with the X regressors and Y responses
  # xcols are the regressors
  # ycols are the responses
  # zcols are the responses that we care about
  # method=GLS other methods not available in this slimmed down version
  # fn is the file name to write to
  # notes: the expected association between permuted hedonic score and gene expression is zero so the expected delta is zero E(E(b.h)-E(b.e))=0-0.
  
  code <- sample(LETTERS,4,replace=TRUE)
  fn_full <- paste(fn,'.',paste(code,collapse=''),'.txt',sep='')
  
  p <- length(ycols)
  Y <- data.matrix(dt[, .SD, .SDcols=ycols])
  b_table <- data.table(NULL)
  # get residuals from covariates
  covs <- setdiff(xcols, zcols)
  X.full <- cbind(rep(1,nrow(dt)),data.matrix(dt[, .SD, .SDcols=xcols]))
  X.red <- cbind(rep(1,nrow(dt)),data.matrix(dt[, .SD, .SDcols=covs]))
  fit.obs <- lm.fit(X.red,Y)
  e <- fit.obs$residuals
  yhat <- fit.obs$fitted.values
  
  rows <- 1:nrow(dt) # observed on first iter and permuted after
  samp <- 'obs'
  for(iter in 1:perms){
    Y.pi <- yhat + e[rows,] # permuted
    dts <- cbind(dt[,.SD,.SDcols=xcols],Y.pi)
    dts[,subject:=factor(.I)] # need to recalc since dt is re-created
    dtlong <- melt(dts,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
    dtlong[,gene:=factor(gene)]
    dtlong <- orderBy(~subject + gene, dtlong)
    
    form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
    if(method=='gls'){
      fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene))
    }
    
    # save stats to table
    # rbinding is slow but relative to the time it takes to do the GLS calculation this is trivial
    # also this works when there are two variables in zcols.
    estimate <- summary(fit1)$tTable[zcols,'Value']
    t_value = summary(fit1)$tTable[zcols,'t-value']
    b_table <- rbind(b_table,data.table(
      permutation=samp,
      Ind.Var=zcols,
      Estimate=estimate,
      t.value=t_value
    ))
    
    if(write_it==TRUE){
      write.table(b_table, fn_full,quote=FALSE,sep='\t',row.names=FALSE)
    }
    # permute rows
    rows <- sample(1:nrow(dt))
    samp <- 'perm'
  }
  
  return(b_table)
}

run_GlobalAncova <- function(dt,xcols,ycols,zcols,perms=10000){
  prob <- NULL
  Y <- as.matrix(dt[, .SD, .SDcols=ycols])
  Yt <- t(Y) # genes as rows

  form.full <- formula(paste('~',paste(xcols,collapse='+'),sep=''))
  model.dat <- dt[, .SD, .SDcols=xcols]
  for(i in 1:length(zcols)){
    prob[zcols[i]] <- GlobalAncova(Yt, form.full, model.dat=model.dat, test.terms=zcols[i], method='permutation',perm=perms)$test.result['p.perm',1]
  }
  return(prob)
}

run_Roast <- function(dt,xcols,ycols,zcols,perms=10000){
  # Roast
  # specific to the FRED datasets analyzing effects of
  # zcols = c('zhedonia', 'zeudaimonia')
  prob <- NULL
  Y <- as.matrix(dt[, .SD, .SDcols=ycols])
  Yt <- t(Y) # genes as rows
  form <- formula(paste('~',paste(xcols,collapse='+'),sep=''))
  design <- model.matrix(form, data=dt)
  colnames(design)[1] <- 'Intercept' # change '(Intercept)' to 'Intercept'
  for(i in 1:length(zcols)){
    Z <- which(colnames(design)==zcols[i])
    prob[zcols[i]] <- roast(y=Yt,design=design,contrast=Z, nrot=perms)$p['UpOrDown','P.Value']
  }
  cont.matrix <- makeContrasts(delta="zhedonia-zeudaimonia",levels=design)
  prob['delta'] <- roast(y=Yt,design=design,contrast=cont.matrix, nrot=perms)$p['UpOrDown','P.Value']
  return(prob)
}

make_ols_table <- function(b_mat){
  # returns ols estimate and bootstrap SE
  # for parametric SE use Obrien?
  # if all=TRUE then return matrix with rows=perms
  delta <- b_mat[,'zhedonia'] - b_mat[,'zeudaimonia']
  b_mat <- cbind(b_mat,delta)
  b_table <- data.table(
    Type = colnames(b_mat),
    Estimate = b_mat[1,],
    SE = apply(b_mat,2,sd)
  )
  b_table[,t:=Estimate/SE]
  b_table[,prob:=2*pt(abs(t),df=n-p-1,lower.tail = FALSE)] # conservative df, upper end should be n*m-p-1
  return(b_table)
}

make_gls_table <- function(fit, method='gls'){
  zcols=c('zhedonia','zeudaimonia')
  if(method=='gls'){gls_table <- summary(fit)$tTable[zcols,]}
  if(method=='gee'){gls_table <- summary(fit)$coefficients[zcols,]}
  # get delta
  coef_names <- names(fit$coefficients)
  hed_i <- which(coef_names=='zhedonia')
  eud_i <- which(coef_names=='zeudaimonia')
  p <- length(coef_names)
  lambda <- numeric(p)
  lambda[hed_i] <- 1
  lambda[eud_i] <- -1
  delta_row <- esticon(fit,lambda)[,c('Estimate','Std.Error','X2.value','Pr(>|X^2|)')]
  row.names(delta_row)[1] <- 'delta'
  colnames(delta_row) <- colnames(gls_table)
  gls_table <- rbind(gls_table, delta_row)
  
  return(gls_table)
}

read_bootstrap.gls_list <- function(fn){ # this is the old format
  the_list <- as.character(read.table(fn)[,1])
  dt <- data.table(NULL)
  for(i in 1:length(the_list)){
    file_name <- the_list[i]
    part_dt <- data.table(read.table(file_name,header=TRUE))
    if(i>1){ #exclude rows with 'obs'
      inc <- which(part_dt[,samp!='obs'])
      part_dt <- part_dt[inc]
    }
    dt <- rbind(dt,part_dt)
  }
  return(dt)
}

read_permutation.gls_list <- function(fn){
  the_list <- as.character(read.table(fn)[,1])
  dt <- data.table(NULL)
  for(i in 1:length(the_list)){
    file_name <- the_list[i]
    part_dt <- data.table(read.table(file_name,header=TRUE))
    if(i>1){ #exclude rows with 'obs'
      inc <- which(part_dt[,permutation!='obs'])
      part_dt <- part_dt[inc]
    }
    dt <- rbind(dt,part_dt)
  }
  return(dt)
}

convert_gls_permutation_to_old_format <- function(long){
  # long is the data.table in the new format
  zhedonia <- long[Ind.Var=='zhedonia']
  zeudaimonia <- long[Ind.Var=='zeudaimonia']
  dt <- cbind(zhedonia[,.(permutation,coeff.zhedonia=Estimate,t.zhedonia=t.value)],
              zeudaimonia[,.(coeff.zeudaimonia=Estimate,t.zeudaimonia=t.value)])
  return(dt)
}

permutation_gls_p_value <- function(res, statistic='t'){
  if(statistic=='t'){
    inc <- which(substr(colnames(res),1,1)=='t')
    # get p-value for delta
    t.value <- res[,.SD,.SDcols=colnames(res)[inc]]
    inc <- which(substr(colnames(res),1,1)=='c')
    coeff.value <- res[,.SD,.SDcols=colnames(res)[inc]]
    # delta.se <- sqrt(se[Type=='hedonic',se]^2 + se[Type=='eudaimonic',se]^2)
    se.value <- coeff.value/t.value
    delta.se <- sqrt(apply(se.value^2,1,sum))
    delta.t <- (coeff.value[,coeff.zhedonia] - coeff.value[,coeff.zeudaimonia])/delta.se
    t.value <- cbind(t.value, delta=delta.t)
    p.value <- apply(abs(t.value),2,function(x) length(which(x >= x[1]))/length(x))
  }
  if(statistic=='c'){
    inc <- which(substr(colnames(res),1,1)=='c')
    p.value <- apply(abs(res[,.SD,.SDcols=colnames(res)[inc]]),2,permutation.p.value)
    p.value <- c(p.value, delta=permutation.p.value(abs(res[,coeff.zhedonia]-res[,coeff.zeudaimonia])))
  }
  return(p.value)
}

permutation_gls_p_value.v2 <- function(res, zcols,statistic='t'){
  # new formatting of input file
  prob <- NULL
  if(statistic=='t'){
    for(iv in zcols){
      t <- res[Ind.Var==iv,t.value]
      prob[iv] <- length(which(abs(t) >= abs(t[1])))/length(t)*100
    }
    # delta
  }
  if(statistic=='c'){
  }
  return(prob)
}

figure_1 <- function(){
  # residual vs. fitted for GLS and GEE
  dt <- copy(dt2015)
  xcols <- get_xcols()
  ycols <- ycols15
  zcols <- c('zhedonia','zeudaimonia')
  fit.gls <- readRDS(paste('FRED15','.gls.rds',sep=''))
  fit.gee <- gls_with_correlated_error(dt,xcols=xcols,ycols,zcols,method='gee')
  # gls residuals vs. fitted
  qplot(x= fitted(fit.gls),y=residuals(fit.gls))
  qplot(x=fit.gee$fitted.values,y=fit.gee$residuals)
  
  dt.gls <- data.table(fitted=fitted(fit.gls),residuals=residuals(fit.gls))
  gg1 <- ggplot(data=dt.gls,aes(x=fitted,y=residuals))
  gg1 <- gg1 + geom_point()
  gg1 <- gg1 + labs(x='Fitted',y = 'Residuals')
  gg1 <- gg1 + ggtitle('A')
  gg1 <- gg1 + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(0.26,.16),plot.margin=unit(x=c(0,0.1,0,0),'cm'))
  gg1
  
  dt.gee <- data.table(fitted=fit.gee$fitted.values,residuals=fit.gee$residuals)
  setnames(dt.gee,old=colnames(dt.gee),new=colnames(dt.gls))
  gg2 <- ggplot(data=dt.gee,aes(x=fitted,y=residuals))
  gg2 <- gg2 + geom_point()
  gg2 <- gg2 + labs(x='Fitted',y = 'Residuals')
  gg2 <- gg2 + ggtitle('B')
  gg2 <- gg2 + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(0.26,.16),plot.margin=unit(x=c(0,0.1,0,0),'cm'))
  gg2
  
  fig_name <- paste('Fig_1_diagnostics.pdf',sep='')
  pdf(fig_name,paper='special',onefile=FALSE,width=6.5,height=3)
  #  postscript('fig_01.eps',horizontal=FALSE,onefile=FALSE,paper='special',height=3,width=6.5)
  showtext.begin()
  print(gg1)
  print(gg2)
  grid.arrange(gg1,gg2,ncol=2,nrow=1)
  showtext.end()
  dev.off()
  
}

figure_2 <- function(){
  # A scatterplot of regression coefficients for x=hedonia y=eudaimonia for GLS bootstrap
  zcols <- c('zhedonia','zeudaimonia')
  # bootstrap gls
  fn <- paste('FRED15.bootstrap.gls.list.v2.txt',sep='')
  gls_boot <- read_bootstrap.gls_list(fn)
  # replace obs with the .nosmoke results since these are what is bootstrapped
  obs <- summary(readRDS(paste('FRED15','-nosmoke.gls.rds',sep='')))$tTable[zcols,'Value']
  gls_boot[samp=='obs',b.zhedonia:=obs[['zhedonia']]]
  gls_boot[samp=='obs',b.zeudaimonia:=obs[['zeudaimonia']]]
  # limit to iter=200
  if(nrow(gls_boot)>200){gls_boot <- gls_boot[1:200,]}
  # get standard errors and CI
  apply(gls_boot[,.SD,.SDcols=c('b.zhedonia','b.zeudaimonia')],2,sd)
  apply(gls_boot[,.SD,.SDcols=c('b.zhedonia','b.zeudaimonia')],2,quantile,probs=c(0.025,0.975))
  gls_boot[,color:=factor(ifelse(samp=='resample',0,1))]
  gls_boot <- orderBy(~color,data=gls_boot)
  r <- cor(gls_boot[,b.zhedonia],gls_boot[,b.zeudaimonia])
  gg <- ggplot(data=gls_boot,aes(x=b.zhedonia,y=b.zeudaimonia,color=color))
  gg <- gg + geom_point(size=2)
  gg <- gg + scale_colour_manual(values=c("grey", "black"), #, "#E69F00,#56B4E9",
                                 labels=c("Resampled", "Observed"))
  gg <- gg + labs(x='Hedonia',y = 'Eudaimonia')
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(0.26,.16),plot.margin=unit(x=c(0,0.1,0,0),'cm'))
  gg
  ggExtra::ggMarginal(gg, type = "histogram")
  gg1 <- gg
  
  fig_name <- paste('Fig_2_gls_boot.pdf',sep='')
  pdf(fig_name,paper='special',onefile=FALSE,width=3,height=3)
  #  postscript('fig_01.eps',horizontal=FALSE,onefile=FALSE,paper='special',height=3,width=6.5)
  showtext.begin()
  print(gg1)
  ggExtra::ggMarginal(gg1, type = "histogram")
  showtext.end()
  dev.off()
  
}

# figure_3 is computed in the GSA1.R script

figure_4 <- function(){
  # bivariate plot of permutation gls effects
  fn <- paste('FRED15', '.permutation.gls.list.v1.txt',sep='')
  gls_perm_res <- read_permutation.gls_list(fn)
  gls_perm_res[,color:=factor(ifelse(permutation=='perm',0,1))]
  #gls_boot <- rbind(gls_boot,data.table(samp='FRED13',b.zhedonia=fred13['zhedonia'],b.zeudaimonia=fred13['zeudaimonia'],t.zhedonia=NA,t.zeudaimonia=NA,color=2))
  gls_perm_res <- orderBy(~color,data=gls_perm_res)
  
  # compare 95% tile
  apply(gls_mc[model=='gls',.SD,.SDcols=c('b','b.hed')],2,quantile,prob=c(0.025,0.975))
  apply(gls_perm_res[,.SD,.SDcols=c('coeff.zeudaimonia','coeff.zhedonia')],2,quantile,prob=c(0.025,0.975))
  
  
  gg <- ggplot(data=gls_perm_res,aes(x=coeff.zhedonia,y=coeff.zeudaimonia, color=color))
  gg <- gg + geom_point(size=2)
  gg <- gg + scale_colour_manual(values=c("grey", "black"), #, "#E69F00,#56B4E9",
                                 labels=c("Permuted", "Observed"))
  gg <- gg + labs(x='Hedonia',y = 'Eudaimonia')
  #gg <- gg + ggtitle('B')
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(0.26,.16),plot.margin=unit(x=c(0,0.1,0,0),'cm'))
  gg
  gg2 <- gg
  
  
  fig_name <- paste('Fig_4_gls_sim.pdf',sep='')
  pdf(fig_name,paper='special',onefile=FALSE,width=3,height=3)
  #  postscript('fig_01.eps',horizontal=FALSE,onefile=FALSE,paper='special',height=3,width=6.5)
  showtext.begin()
  print(gg2)
  ggExtra::ggMarginal(gg2, type = "histogram")
  showtext.end()
  dev.off()
  
}
