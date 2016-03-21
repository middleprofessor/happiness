# simulation of Fredrickson, Barbara L., et al. "A functional genomic perspective on human well-being." Proceedings of the National Academy of Sciences 110.33 (2013): 13684-13689.
# Jeffrey A. Walker
# December 6, 2016
# created git site on Feb 16 and copied script file. Original script as of Feb 16 is archived as Fredrickson_multivariate_lm.init.R in the folder 'init code'
 
library(data.table)
library(car)
library(mvtnorm)
library(lmPerm)
library(nlme)
library(lme4)
library(reshape2)
library(doBy)
# plot libraries
library(ggplot2)
library(showtext) # needed for eps fonts to add Arial to .eps
font.add('Arial',regular='Arial.ttf')
library(gridExtra)
library(pls)
library(MRCE) # multivariate regression with correlated error
library(gee) # generalized estimating equations
library(geepack) #ditto
library(GlobalAncova) # GlobalAncova
library(limma) # Roast
library(CCA) # canonical correlation Analysis

# lmperm downloaded from https://github.com/kabacoff/RiA2/tree/master/lmPerm
# using R-Studio chose Tools > Install Packages and the Install From > Package Archive pop-up menu.

# originally coefficients computed raw then reversed during analysis but to compute the permutation GLH, these need to be reversed, so...
# January 13 recomputed by reversing expression levels and not coefficients
do_Fredrickson <- function(){
  
  clean <- TRUE
  if(clean==FALSE){
    clean_cole1()
    clean_cole2()
  }
  
  # generate resampled files to compute results
  are_tests_run <- FALSE
  are_gls_permutations_run <- TRUE
  are_lme_permutations_run <- TRUE
  are_gls_parametric_tests_run <- TRUE
  are_gee_permutations_run <- TRUE
  are_type_1_simulations_run <- TRUE
  
  fn <- 'cole1_clean.txt'
  dt2013 <- read_file(fn,year=2013)
  dt2013[,zhedonia:=scale(zhedonia)]
  dt2013[,zeudaimonia:=scale(zeudaimonia)]
  dt2013 <- contrast_coefficients(dt2013) # convert to CTRA response
  
  fn <- 'cole2_clean.txt'
  dt2015 <- read_file(fn,year=2015)
  dt2015[,zhedonia:=scale(zhedonia)]
  dt2015[,zeudaimonia:=scale(zeudaimonia)]
  dt2015 <- contrast_coefficients(dt2015)  # convert to CTRA response
  
  fn <- 'cole1_clean.txt'
  dt1 <- read_file(fn,year=2015)
  fn <- 'cole2_clean.txt'
  dt2 <- read_file(fn,year=2015)
  # illness in FRED2015 is averaged over the 13 categories so divide illness in FRED13 by 13 to be in same scale as in FRED15. Even with this the range of FRED15 is about 2X that of FRED13.
  dt1[,illness:=(illness/13)]
  dtCombi <- rbind(dt1,dt2)
  # rescale zhedonia and zeudaimonia
  dtCombi[,zhedonia:=scale(zhedonia)]
  dtCombi[,zeudaimonia:=scale(zeudaimonia)]
  dtCombi <- contrast_coefficients(dtCombi)  # convert to CTRA response
  
  if(are_tests_run==FALSE){
    # note that in FRED13, smoke has only 6/77 scored as 1 and some bootstraps will entirely miss this. See note in bootstrap_t_test
    permutation_t_tests(dt2013,which_file='FRED13',niter=2000)
    permutation_t_tests(dt2015,which_file='FRED15',niter=2000)
    permutation_t_tests(dtCombi,which_file='FRED.Combi',niter=2000)
    
    boot_t_tests(dt2013,which_file='FRED13',niter=2000)
    boot_t_tests(dt2015,which_file='FRED15',niter=2000)
    boot_t_tests(dtCombi,which_file='FRED.Combi',niter=2000)
    
    bootstrap_models(dt2013,which_file='FRED13',tests=c('gee'), niter=2000)
    bootstrap_models(dt2015,which_file='FRED15',tests=c('gee'), niter=2000)
    bootstrap_models(dtCombi,which_file='FRED.Combi',tests=c('gee'), niter=2000)
 
    bootstrap_models(dt2013,which_file='FRED13',tests=c('gls'), niter=100)
    bootstrap_models(dt2015,which_file='FRED15',tests=c('gls'), niter=100)
    bootstrap_models(dtCombi,which_file='FRED.Combi',tests=c('gls'), niter=100)

    bootstrap_obrien(dt2013,fn='FRED13.obrien', niter=1000,initer=1000)
    bootstrap_obrien(dt2015,fn='FRED15.obrien', niter=1000,initer=1000)
    bootstrap_obrien(dtCombi,fn='FRED.Combi.obrien', niter=1000,initer=1000)
    
    permutation_obrien(dt2013,fn='FRED13.obrien', niter=2000,initer=1000)
    permutation_obrien(dt2015,fn='FRED15.obrien', niter=2000,initer=1000)
    permutation_obrien(dtCombi,fn='FRED.Combi.obrien', niter=2000,initer=1000)
  }
  if(are_gls_permutations_run==FALSE){
      do_gls_tests(dt2013,which_file='FRED13', method='gls', niter=101)
      do_gls_tests(dt2015,which_file='FRED15', method='gls', niter=101)
      do_gls_tests(dtCombi,which_file='FRED.Combi', method='gls', niter=101)
  }
  if(are_lme_permutations_run==FALSE){
    do_gls_tests(dt2013,which_file='FRED13', method='lme', niter=101) # doens't converge
    do_gls_tests(dt2015,which_file='FRED15', method='lme', niter=101)
    do_gls_tests(dtCombi,which_file='FRED.Combi', method='lme', niter=101)
  }
  if(are_gee_permutations_run==FALSE){
    do_gls_tests(dt2013,which_file='FRED13', method='gee', niter=2000)
    do_gls_tests(dt2015,which_file='FRED15', method='gee', niter=2000)
    do_gls_tests(dtCombi,which_file='FRED.Combi', method='gee', niter=2000)
  }
  if(are_gls_parametric_tests_run==FALSE){
    saveRDS(gls_with_correlated_error(dt2013, year=2013), "FRED13.gls.rds")
    saveRDS(gls_with_correlated_error(dt2015, year=2015), "FRED15.gls.rds")
    saveRDS(gls_with_correlated_error(dtCombi, year=2015), "FRED.Combi.gls.rds")
    saveRDS(gls_with_correlated_error(dt2013, year=2015), "FRED13.yr2015.gls.rds")
  }
  if(are_type_1_simulations_run==FALSE){
    simulate_type1_error(run_simulations=TRUE)
  }else{
    type_1_error_table <- simulate_type1_error(run_simulations=FALSE)
    }

  # compute all tables
  which_file_list <- c('FRED13', 'FRED15', 'FRED.Combi')
  hypotheses <- c('zhedonia','zeudaimonia','delta')
  tests <- c('bootstrap_t','perm_t','perm_gls','perm_gee')
  zcols <- c('zhedonia','zeudaimonia')
  
  naive_p_table <- data.table(NULL) # naive t-test p-values
  ols_table <- data.table(NULL) # average beta over m responses
  gls_table <- data.table(NULL) # common effect estimated by GLS
  gee_table <- data.table(NULL) # common effect estimated by GEE
  delta_table <- data.table(NULL)
    
  effect_matrix <- matrix(0,nrow=3,ncol=2) # matrix of backtransformed effects
  colnames(effect_matrix) <- c('zhedonia','zeudaimonia')
  row.names(effect_matrix) <- which_file_list
  resid_df <- numeric(3) # put the df of the t-tests for the gls here
  names(resid_df) <- which_file_list
  gls_perm_type_1_error_table <- data.table(NULL) # not the error of the perm test but the error of the parametric computed from the perm results
  
  effective_size_table <- numeric(3)
  names(effective_size_table) <- which_file_list

  for(which_file in which_file_list){
    if(which_file=='FRED13'){dt <- copy(dt2013)}
    if(which_file=='FRED15'){dt <- copy(dt2015)}
    if(which_file=='FRED.Combi'){dt <- copy(dtCombi)}
    if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
    ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
    R <- cor(as.matrix(dt[,.SD,.SDcols=ycols]))
    mean(abs(R[lower.tri(R)]))
    max(abs(R[lower.tri(R)]))
    
    # OLS table
    fn <- paste('coeff.',which_file,'.permutation.txt',sep='')
    ols_perm <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
    ols_perm[, beta:=apply(ols_perm[, .SD, .SDcols=ycols], 1, mean)]
    beta <- ols_perm[data=='obs',.(Type,beta)]
    fn <- paste(which_file,'.bootstrap.txt',sep='')
    boot_coeffs <- data.table(read.table(fn,header=TRUE))  # read in regression coefficients
    boot_coeffs[, beta:=apply(boot_coeffs[, .SD, .SDcols=ycols], 1, mean)]
    beta.boot <- boot_coeffs[data=='obs',.(Type,beta)]
    part_table <- get_ols_table(boot_coeffs,which_file)[,.(Type,Data,beta,SE)] # delete CIs
    ols_boot_p <- smart_t_stats(boot_coeffs)
    part_table[, boot_p:=ols_boot_p]
    ols_perm_p <- permutation_t(ols_perm,statistic='t')
    part_table[, perm_p:=ols_perm_p]
    # obrien permutation
    fn <- paste(which_file,'.obrien.perm.txt',sep='')
    obrien_res <- data.table(read.table(fn,header=TRUE)) # O'Brien T statistics
    obrien_p <- apply(abs(obrien_res),2,permutation.p.value)
    part_table[, obrien_p:=obrien_p]
    # obrien t
    #fn <- paste(which_file,'.obrien.boot.txt',sep='')
    #obrien_res <- data.table(read.table(fn,header=TRUE)) # O'Brien T statistics
    #obrien_se <- apply(obrien_res,2,quantile,prob=c(0.025,0.975))
    # roast
    rot_p <- roast_it(dt)
    part_table[, rot_p:=rot_p]
    
    ols_table <- rbind(ols_table,part_table)
    
    delta_table <- rbind(delta_table, data.table(Test='OLS', Data=which_file, Delta=part_table[Type=='delta', beta], boot_SE=part_table[Type=='delta', SE], boot_p=ols_boot_p['delta'], perm_p=ols_perm_p['delta'], obrien_p=obrien_p['Tdelta'], rot_p=rot_p['delta']))
     
    # effect size back in units of FRED13
    # compute back-transformed effect size
    effect_matrix[which_file,] <- effect_size(dt) # effect size for observed data
    # do naive t test first following FRED13
    naive_t_table <- naive_t_stats(dt) #(reported in FRED13: eudaimonic, P = 0.0045; hedonic, P = 0.0047)
    naive_p_table <- rbind(naive_p_table, data.table(data=which_file,naive_t_table[stat=='p.value',]))

    # gls coefficients and p-values
    fit <- readRDS(paste(which_file,'.gls.rds',sep=''))
    gls_res <- summary(fit)$tTable[zcols,]
     # get bootstrap CIs of gls
    gls_boot_se <- c(NA, NA, NA)
    names(gls_boot_se) <- hypotheses
    if(which_file == 'FRED.Combi'){
      gls_boot_res <- read_bootstrap.gls_list('FRED.Combi.bootstrap.gls.list.txt')
      gls_boot_res[,delta:=zhedonia-zeudaimonia]
      gls_boot_se <- apply(gls_boot_res[,.SD,.SDcols=hypotheses],2,sd)
 #     apply(gls_boot_res[,.SD,.SDcols=zcols],2,quantile,probs=c(0.025,0.975))
    }
    #compute permutation.gls stats
    fn <- paste(which_file, '.permutation.gls.list.txt',sep='')
    gls_perm_res <- permutation.gls.p.value(read_permutation.gls_list(fn), statistic='t')
    gls_table <- rbind(gls_table, data.table(
      Type=row.names(gls_res),
      Data=which_file,
      gls_res[,c('Value','Std.Error','p-value')],
      boot_SE=gls_boot_se[zcols],
      perm_p=gls_perm_res[c('t.zhedonia','t.zeudaimonia')]
      ))
    # uncomment to get CI on p-value for paper
    # permutation.gls.p.value(res, statistic='t',do_ci=TRUE) # get 95% CI on p-values

    # do gee
    gee_res <- gee_delta(gee_with_correlated_error(dt))
    gee_boot_se <- gee_bootstrap_se(which_file)
    gee_p <- gee_permutation_p(which_file)
    gee_table <- rbind(gee_table, data.table(
      Type=c(zcols,'delta'),
      Data=which_file,
      beta=gee_res[, 'Estimate'],
      se_robust=gee_res[, 'Std.err'],
      se_boot=gee_boot_se,
      wald_boot=gee_res[, 'Estimate']^2/gee_boot_se^2,
      boot_p=wald_test(gee_res[, 'Estimate'],gee_boot_se),
      perm_p=gee_p
      ))
  
  }
  
  # clean tables
  gee_table[, Type:=factor(Type)]
  gee_table <- round_dt(orderBy(~-Type, gee_table), 3)
  gee_table[,boot_p:=round(boot_p,2)]
  gee_table[,perm_p:=round(perm_p,2)]
  gee_table <- gee_table[Type!='delta'] # delta now in its own table
  
  ols_table[, Type:=factor(Type)]
  ols_table <- round_dt(orderBy(~-Type,ols_table),3)
  ols_table[,boot_p:=round(boot_p,2)]
  ols_table[,perm_p:=round(perm_p,2)]
  ols_table[,rot_p:=round(rot_p,2)]
  ols_table[,obrien_p:=round(obrien_p,2)]
  ols_table <- ols_table[Type!='delta'] # delta now in its own table
  
  gls_table[, Type:=factor(Type)]
  gls_table <- round_dt(orderBy(~-Type, gls_table), 3)
  gls_table[,perm_p:=round(perm_p,2)]
 
  delta_table[, Test:=factor(Test)]
  delta_table <- round_dt(orderBy(~-Test, delta_table), 3)
  delta_table[,perm_p:=round(perm_p,2)]
  delta_table[,boot_p:=round(boot_p,2)]
  delta_table[,rot_p:=round(rot_p,2)]
  delta_table[,obrien_p:=round(obrien_p,2)]
 
  effect_matrix <- data.table(data=row.names(effect_matrix),(round(effect_matrix,1)))
  naive_p_table <- round_dt(naive_p_table,3)
  naive_p_table <- naive_p_table[,.SD,.SDcols=c('data','zhedonia','zeudaimonia')]
  naive_effect_table <- merge(effect_matrix,naive_p_table,by='data')

  naive_effect_table
  gls_table
  ols_table
  gee_table
  delta_table
  
  write.table(naive_effect_table,'naive_effect_table.txt',sep='\t',quote=FALSE,row.names=FALSE)
  write.table(gls_table,'gls_table.txt',sep='\t',quote=FALSE,row.names=FALSE)
  write.table(ols_table, 'ols_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
  write.table(gee_table, 'gee_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
  write.table(delta_table, 'delta_table.txt',quote=FALSE,row.names=FALSE,sep='\t')
  # end function
  
}

figure_0 <- function(){
  # plot of residual vs fitted for GLS and GEE for combined
  dt <- copy(dtCombi)
  year <- 2015
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  dtlong <- get_dtlong(dt, year=2015, center=TRUE)
  dtlong <- orderBy(~subject + gene,dtlong)
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  #form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))

  fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject, waves=gene, corstr='exchangeable', std.err="san.se")
  summary(fit.geeglm)$coefficients[zcols,]
  qplot(fit.geeglm$fitted.values,fit.geeglm$residuals)
  # compare this to the GLS fit
  
  fit.gls.hetero <- readRDS(paste(which_file='FRED.Combi','.gls.rds',sep=''))
  qplot(fit.gls.hetero$fitted,fit.gls.hetero$residuals)
  
  # is the bias due to the heterogenous, what about exchangeable?
  fit.gls.homo <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  qplot(fit.gls.homo$fitted,fit.gls.homo$residuals)
  # yes
  
  # what about lme with hetergenous?
  fit.lme <- lme(form, random = ~1|subject, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(maxIter=100, msMaxIter = 500, tolerance=1e-6, msVerbose = FALSE)) # tolerance=1e-6 default
  saveRDS(fit.lme, "FRED.Combi.lme.rds")
  summary(fit.lme)$tTable[zcols,]
  qplot(fit.lme$fitted[,'fixed'],fit.lme$residuals[,'fixed'])
  qplot(fit.lme$fitted[,'subject'],fit.lme$residuals[,'subject'])
  
}

figure_1 <- function(){
  # A scatterplot of regression coefficients for x=hedonia y=eudaimonia for two of the randomly permuted runs using the FRED.Combi data
  which_file <- 'FRED13'
  fn <- paste(which_file,'.permutation.txt',sep='')
  part <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  dt <- data.table(data=which_file,part)
  which_file <- 'FRED15'
  fn <- paste(which_file,'.permutation.txt',sep='')
  part <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  dt <- rbind(dt,data.table(data=which_file,part),fill=TRUE)
  which_file <- 'FRED.Combi'
  fn <- paste(which_file,'.permutation.txt',sep='')
  part <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  dt <- rbind(dt,data.table(data=which_file,part),fill=TRUE)
  dt <- dt[,.(data,I,Type,IL1A)]
  
  dtcast <- merge(dt[Type=='hedonic',], dt[Type=='eudaimonic',],by=c('I','data'))
  setnames(dtcast,c('I','gene','Type.x','hedonic','Type.y','eudaimonic'))
  dtcast[I!=1,I:=0] # color observed differently from permuted
  dtcast <- rbind(dtcast[I==0,],dtcast[I==1,]) # make the observed data last
  dtcast[,I:=factor(as.character(I))]
  dtcast[,cor(hedonic,eudaimonic),by=I]
  gg <- ggplot(data=dtcast,aes(x=hedonic,y=eudaimonic,color=I))
  gg <- gg + geom_point()
  gg <- gg + scale_colour_manual(values=c("grey", "black"), #, "#E69F00,#56B4E9"
                                 labels=c("Permuted", "Observed"))
  #gg <- gg + scale_color_brewer(name=NULL,labels=c("Permuted", "Observed"), palette="Paired")
  gg <- gg + labs(x='Hedonia',y = 'Eudaimonia')
  gg <- gg + ggtitle('A')
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(.85,.9),plot.margin=unit(x=c(0,0.1,0,0),'cm'))
  gg
  gg1 <- gg
  
  fn <- paste('FRED13', '.permutation.gls.list.txt',sep='')
  res <- read_permutation.gls_list(fn)
  setnames(res,c('permutation','hedonic','eudaimonic','t.zhedonia','t.zeudaimonia'))
  dt <- copy(res)
  fn <- paste('FRED15', '.permutation.gls.list.txt',sep='')
  res <- read_permutation.gls_list(fn)
  setnames(res,c('permutation','hedonic','eudaimonic','t.zhedonia','t.zeudaimonia'))
  dt <- rbind(dt,res)
  fn <- paste('FRED.Combi', '.permutation.gls.list.txt',sep='')
  res <- read_permutation.gls_list(fn)
  setnames(res,c('permutation','hedonic','eudaimonic','t.zhedonia','t.zeudaimonia'))
  dt <- rbind(dt,res)
  dt[permutation=='perm',I:=0]
  dt[permutation=='obs',I:=1]
  dt[,I:=factor(as.character(I))] # make I factor
  dt <- rbind(dt[I=='0',], dt[I=='1']) # plot perm first
  gg <- ggplot(data=dt,aes(x=hedonic,y=eudaimonic,color=I))
  gg <- gg + geom_point()
  gg <- gg + scale_colour_manual(values=c("grey", "black"), #, "#E69F00,#56B4E9"
                                 labels=c("Permuted", "Observed"))
  gg <- gg + labs(x='Hedonia',y = 'Eudaimonia')
  gg <- gg + ggtitle('B')
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(.85,.9),plot.margin=unit(x=c(0,0.1,0,0),'cm'))
  gg
  gg2 <- gg
  
  
  fig_name <- paste('Fig_1.pdf',sep='')
  pdf(fig_name,paper='special',onefile=FALSE,width=6.5,height=3)
  #  postscript('fig_01.eps',horizontal=FALSE,onefile=FALSE,paper='special',height=3,width=6.5)
  showtext.begin()
  print(gg1)
  print(gg2)
  grid.arrange(gg1,gg2,ncol=2,nrow=1)
  showtext.end()
  dev.off()
  
}

effective_size <- function(N, R){
  # N is the number of subjects
  # R is the within subject error matrix
  rho <- mean(abs(R[lower.tri(R)]))
  n <- nrow(R)
  Nn <- N*n
  Ne <- Nn/(1 + rho*(n-1))
  return(Ne)
}

clean_cole1 <- function(){
  # cole1 downloaded from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45330
  # do not confuse coleI with brown. brown is a cleaned version of cole1.
  # cole1 is a 34591 x 80 data set of gene expression from the original PNAS paper. The paper analyzed only 53 immune genes
  # cole1 is a xxx x xxx data set of gene expression from the follow-up paper in PNAS
  cole1 <- read.table('files_needed_to_make_cole1/cole1_data_matrix.txt',header=TRUE,sep='\t',row.names='ID_REF')
  cole1 <- t(cole1)
  case_ID <- row.names(cole1)
  cole1 <- data.table(cole1)
  ycols <- c(pro_inflam_genes(year=2013),antibody_genes(),ifn_genes())
  dt <- cole1[,.SD, .SDcols=ycols]
  dt <- log2(dt)
  
  # compare to brown
  #y <- cole1[,IFIT1L]
  #target <- brown[,IFIT1L] # okay so Brown data are just cole1 log2(cole1)

  # read regresssors
  #NOTE recoded 'Missing' as NA and '.' as NA
  x <- read.table('files_needed_to_make_cole1/cole1_regressors.txt',sep='\t')
  col_labels <- x[,1]
  x <- t(x[,-1])
  colnames(x) <- col_labels
  x <- data.table(x)
  
  # convert to integer
  x[,male:=as.integer(male)]
  x[,age:=as.integer(age)]
  x[,white:=as.integer(white)]
  x[,alcohol:=as.integer(alcohol)]
  x[,smoke:=as.integer(smoke)]
  x[,illness:=as.integer(illness)]
  
  # save and reopen
  write.table(x,'x.temp',quote=FALSE,row.names=FALSE,sep='\t')
  x <- data.table(read.table('x.temp',header=TRUE,sep='\t'))
  
  # check that zhedonia are correct z-scores. No they are not. Cole1 was scaled prior to culling some data.
  check <- data.table(zhedonia=x[,zhedonia], zscore=scale(x[,hedonia]))
  cov(check,use='pairwise.complete.obs')
  # rescale hedonia and eudaimonia
  # x[, zhedonia:=scale(hedonia)]
  # x[, zeudaimonia:=scale(eudaimonia)]
  
  cole1 <- cbind(x,dt)
  
  # rows 74 and 79 are completely NA for the x data so delete these
  inc <- setdiff(1:nrow(cole1), c(74,79))
  cole1 <- cole1[inc,]
  inc <- which(is.na(cole1),arr.ind=TRUE)
  # row 60 is missing a value for "alcohol" 
  inc <- setdiff(1:nrow(cole1), c(60))
  cole1 <- cole1[inc,]
  
  write.table(cole1,'cole1_clean.txt',quote=FALSE,sep='\t',row.names=FALSE)
  
}

clean_cole2 <- function(){
  # cole2 downloaded from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55762
  # first create cole II data
  cole2 <- read.table('files_needed_to_make_cole2/cole2_data_matrix.txt',header=TRUE,sep='\t',row.names='ID_REF')
  cole2 <- t(cole2)
  case_ID <- row.names(cole2)
  cole2 <- data.table(cole2)
  
  ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  dt <- cole2[,.SD, .SDcols=ycols]
  dt <- log2(dt)
  
  x <- read.table('files_needed_to_make_cole2/cole2_regressor_clean.txt',sep='\t')
  col_labels <- x[,1]
  x <- t(x[,-1])
  colnames(x) <- col_labels
  x <- data.table(x)
  
  # need to make hedonic and eudaimonic z scores
  x[, zhedonia:=scale(hedonic)]
  x[, zeudaimonia:=scale(eudaimonic)]
  
  # is illness variance-scaled? - no but it is on a different scale from the 2013 data.
  mean(x[,illness])
  sd(x[,illness])
  
  cole2 <- cbind(x,dt)
  write.table(cole2,'cole2_clean.txt',quote=FALSE,sep='\t',row.names=FALSE)
}

correlation_in_Y <- function(){
  dt2013 <- read_file(fn,year=2013)
  ycols <- c(pro_inflam_genes(2013),antibody_genes(),ifn_genes())
  Y <- dt2013[,.SD,.SDcols=ycols]
  R <- cor(Y)
  mean(R[lower.tri(R)])
  max(R[lower.tri(R)])
  min(R[lower.tri(R)])
  
  Ycc <- contrast_coefficients(copy(Y))
  Rcc <- cor(Ycc)
  mean(R[lower.tri(Rcc)])
  max(R[lower.tri(Rcc)])
  min(R[lower.tri(Rcc)])
  
}

compare_cole1_brown <- function(){
  # The composite mean is different in brown and cole1, why?
  x.cole1 <- cole1[,.SD,.SDcols=xcols]
  x.brown <- brown[,.SD,.SDcols=xcols]
  
  # The difference in x seems to be how the 'White=4' was coded. If this is coded as NA (my original analysis) I get numbers close to Frederickson/Brown, If coded as 0 (which is the corrected version) I get my cole1 numbers, at least to the 3rd decimal place. Note that my version of the downloaded data matrix did NOT have any White=4 so this must have been corrected.
  
  # There is also the the alcohol='.' in my downloaded version. What is this in Brown?
  data.frame(cole1=cole1[,alcohol],brown=brown[,alcohol]) #NA in both
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

read_permutation.gls_list <- function(fn){
  the_list <- as.character(read.table(fn)[,1])
  res <- data.table(NULL)
  for(i in 1:length(the_list)){
    sfn <- the_list[i]
    res <- rbind(res,read.table(sfn,header=TRUE))
  }
  
  # check for multiple rows with perm=obs
  obs_rows <- which(res[,permutation]=='obs')
  exc <- setdiff(obs_rows,obs_rows[1])
  if(length(exc>=1)){res <- res[-exc,]}
  return(res)
}

read_bootstrap.gls_list <- function(fn){
  the_list <- as.character(read.table(fn)[,1])
  res <- data.table(NULL)
  for(i in 1:length(the_list)){
    sfn <- the_list[i]
    res <- rbind(res,read.table(sfn,header=TRUE))
  }
  
  # check for multiple rows with perm=obs
  obs_rows <- which(res[,samp]=='obs')
  exc <- setdiff(obs_rows,obs_rows[1])
  if(length(exc>=1)){res <- res[-exc,]}
  return(res)
}

gee_bootstrap_se <- function(which_file='FRED13'){
  if(which_file=='FRED13'){gee_fn <-'FRED13.AGSA.gee.table.txt'}
  if(which_file=='FRED15'){gee_fn <-'FRED15.LWDK.gee.table.txt'}
  if(which_file=='FRED.Combi'){gee_fn <-'FRED.Combi.OQKZ.gee.table.txt'}
  gee_boot <- data.table(read.table(gee_fn,header=TRUE,sep='\t'))
  gee_boot[, delta:=zhedonia-zeudaimonia]
  ses <- apply(gee_boot[, .SD, .SDcols=c('zhedonia','zeudaimonia','delta')],2,sd)
  # CIs <- gee_boot[,.(lwr=quantile(delta,0.025),upr=quantile(delta,0.975))]
  return(ses)
}

gee_permutation_p <- function(which_file){
  if(which_file=='FRED13'){gee_fn <-'FRED13.permutation.gee.ZHRW.gee.txt'}
  if(which_file=='FRED15'){gee_fn <-'FRED15.permutation.gee.ULKU.gee.txt'}
  if(which_file=='FRED.Combi'){gee_fn <-'FRED.Combi.permutation.gee.MEIF.gee.txt'}
  gee_perm <- data.table(read.table(gee_fn,header=TRUE,sep='\t'))
  gee_perm[, delta:=zhedonia-zeudaimonia]
  gee_p <- apply(abs(gee_perm[,.SD,.SDcols=c('zhedonia', 'zeudaimonia','delta')]), 2, permutation.p.value)
  return(gee_p)
}

round_dt <- function(df, digits) {
  df <- data.frame(df)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- (round(df[,nums], digits = digits))
  df <- data.table(df)
  (df)
}

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  (df)
}


permutation_t_tests <- function(dt,which_file,niter=2000){
  # dt is a data.table of the responses and regressors
  # which_file is the file FRED13 or FRED15 Fredrickson et. al.

  # get xcols and ycols
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  permutation_test_anderson(dt,xcols,ycols,niter=niter,write_it=TRUE,fn=paste(which_file,'.permutation.txt',sep=''))
}

boot_t_tests <- function(dt,which_file,niter=2000, exclude_smoke=FALSE){
  # dt is a data.table of the responses and regressors
  # which_file is the file FRED13 or FRED15 Fredrickson et. al.
  # exclude_smoke=TRUE, see below)
  
  # get xcols and ycols
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  # note that in FRED13, smoke has only 6/77 scored as 1 and some bootstraps will entirely miss this. I deal with this by excluding smoke from xcols in the bootstrap samples with fewer than 2 smoke=1
  if(exclude_smoke==TRUE){xcols <- setdiff(xcols,'smoke')}
  bootstrap_test(dt,xcols,ycols,niter=niter,write_it=TRUE, fn=paste(which_file,'.bootstrap.txt',sep=''))
  
}


do_gls_tests <- function(dt,which_file,method='lme',niter=101){
  # dt is a data.table of the responses and regressors
  # which_file is the file FRED13 or FRED15 Fredrickson et. al.
  # exclude_smoke=TRUE, see below)
  
  # get xcols and ycols
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  #gls permutation
  zcols <- c('zhedonia','zeudaimonia')
  permutation_gls(dt,xcols=xcols,ycols=ycols,zcols=zcols,method=method,niter=niter,do_obs=TRUE,write_it=TRUE,fn=paste(which_file,'.permutation.',method,sep=''))
  
  #gls permutation on individual measures of happiness
#  zcols <- c('zhedonia','zeudaimonia')
#  happy <- 'zeudaimonia' # limit analysis to zcols
#  red_xcols <- setdiff(xcols,setdiff(zcols,happy)) # remove the other happy from xcols
#  permutation_gls(dt,xcols=red_xcols,ycols=ycols,zcols=happy,niter=80,do_obs=FALSE,write_it=TRUE,fn=paste(paste(happy,which_file,sep=''),'.permutation.gls',sep=''))
  
}

get_xcols <- function(){
  # These are the regressors
  xcols <- c('male' , 'age' , 'white' , 'bmi' , 'alcohol', 'smoke' , 'illness' , 'cd3d' , 'cd3e' , 'cd4', 'cd8a' , 'fcgr3a', 'cd19' , 'ncam1' , 'cd14' , 'zhedonia' , 'zeudaimonia')
  return(xcols)
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


coeff_vector <- function(dt,xcols,ycols){
  # Frederickson et al computation of the regression coefficient of each gene expression on the set of regressors
  # dt is a data.table
  # xcols are the cols of dt that are the regressors
  # ycols are the cols of dt that are the univariate response
  p <- length(ycols)
  b.hed <- numeric(p)
  b.eud <- numeric(p)
  se.hed <- numeric(p) # doing this just to get the vector of se's for modeling
  for(j in 1:p){
    form <- formula(paste(ycols[j],paste(xcols,collapse='+'),sep='~'))
    fit <- lm(form, data=dt)
    b.hed[j] <- coefficients(fit)['zhedonia']
    b.eud[j] <- coefficients(fit)['zeudaimonia']
    se.hed[j] <- coefficients(summary(fit))['zhedonia','Std. Error']
  }
  
  # check ... check!
#  Y <- as.matrix(dt[,.SD,.SDcols=ycols])
#  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
#  fitmv <- lm(form, data=dt)
#  b.hedm <- coefficients(fitmv)['zhedonia',]
#  b.eudm <- coefficients(fitmv)['zeudaimonia',]
#  data.table(uni=b.hed,multi=b.hedm)
  
  coeff_table <- data.frame(b.hed=b.hed,b.eud=b.eud)
  rownames(coeff_table) <- ycols
  return(coeff_table)
}

get_ols_table <- function(cc, which_file){
  # cc input are contrast coefficients
  # returns mean beta and SE and CIs
  if('IL6' %in% colnames(cc)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  b_table <- copy(cc)
  b_table[, beta:=apply(.SD,1,mean),.SDcols=ycols]
  means <- cbind(b_table[Type=='hedonic',.(data,hedonic=beta)], b_table[Type=='eudaimonic',.(eudaimonic=beta)])
  means[, delta:=hedonic-eudaimonic]
  beta_bar <- apply(means[data=='obs', .SD, .SDcols=c('hedonic','eudaimonic','delta')],2,mean) # using apply because it returns the correct format
  beta_se <- apply(means[, .SD, .SDcols=c('hedonic','eudaimonic','delta')],2,sd)
  beta_ci <- apply(means[, .SD, .SDcols=c('hedonic','eudaimonic','delta')],2,quantile, probs=c(0.025,0.975))
  ols_table <- data.table(
    Type=c('zhedonia','zeudaimonia','delta'),
    Data=which_file,
    beta=beta_bar,
    SE=beta_se,
    lwr=beta_ci['2.5%',],
    upr=beta_ci['97.5%',]
  )
   return(ols_table)
}

permutation_test <- function(dt,res,xcols,ycols,niter=1000, write_it=FALSE,fn){
  # permutation test of delta - the mean difference between hed and eud coeffs
  # dt is original data
  # res is the resulting coefficients of the regression
  # xcols are the regressors
  # ycols are the responses
  # fn is the file name to write to
  # notes: the expected association between permuted hedonic score and gene expression is zero so the expected delta is zero E(E(b.h)-E(b.e))=0-0
  # *** note - following Anderson 2001 this doesn't maintain exchangeability because all of Y is permuted instead of just error from X. So do permutation_test_anderson
  
  p <- length(ycols)
  hed_matrix <- matrix(0,nrow=niter,ncol=p) # matrix of coefficients
  eud_matrix <- matrix(0,nrow=niter,ncol=p) # matrix of coefficients
  colnames(hed_matrix) <- ycols
  colnames(eud_matrix) <- ycols
  hed_matrix[1,] <- res$b.hed
  eud_matrix[1,] <- res$b.eud
  
  for(iter in 2:niter){
    perm <- sample(1:nrow(dt))
    X <- dt[perm,.SD, .SDcols=xcols]
    Y <- dt[,.SD, .SDcols=ycols]
    permuted_dt <- cbind(X,Y)
    sim.res <- coeff_vector(permuted_dt,xcols,ycols)
    hed_matrix[iter,] <- sim.res$b.hed
    eud_matrix[iter,] <- sim.res$b.eud
  }
  b_matrix <- data.table(I=rep(1:niter,2), Type=rep(c('hedonic','eudaimonic'),each=niter), rbind(hed_matrix,eud_matrix))
  
  if(write_it==TRUE){write.table(b_matrix,fn,quote=FALSE,sep='\t',row.names=FALSE)}
  return(NULL)
}

permutation_test_anderson <- function(dt,xcols,ycols,niter=1000, write_it=FALSE,fn, multi=FALSE){
  # permutation test based on Anderson and Robinson 2001
  # dt is a data.table with the X regressors and Y responses
  # xcols are the regressors
  # ycols are the responses
  # fn is the file name to write to
  # multi: if TRUE then do GLH tests
  # notes: the expected association between permuted hedonic score and gene expression is zero so the expected delta is zero E(E(b.h)-E(b.e))=0-0.
  
  # returns two different results
  # bmatrix is a matrix of the 53 coefficients for each permuation (observed in row I=1)

  p <- length(ycols)
  Y <- scale(as.matrix(dt[, .SD, .SDcols=ycols]))
  dt[,zhedonia:=scale(zhedonia)]
  dt[,zeudaimonia:=scale(zeudaimonia)]
  
  # coefficients (save t-value in addition to coefficients)
  coeffs.hed <- matrix(0,nrow=niter,ncol=p)
  coeffs.eud <- matrix(0,nrow=niter,ncol=p)
  colnames(coeffs.hed) <- ycols
  colnames(coeffs.eud) <- ycols
  t.hed <- matrix(0,nrow=niter,ncol=p)
  t.eud <- matrix(0,nrow=niter,ncol=p)
  colnames(t.hed) <- ycols
  colnames(t.eud) <- ycols

  zcols <- c('zhedonia','zeudaimonia')
  xcols2 <- setdiff(xcols, zcols)
  # Get residuals from X (so excluding zhedonia and zeudaimonia)
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit.obs <- lm(form, data=dt)
  yhat <- predict(fit.obs) # Yhat = aX
  rows <- 1:nrow(dt) # observed on first iter and permuted after
  for(iter in 1:niter){
    e <- residuals(fit.obs)[rows,] # permuted Yhat - aX
    # yhat is that predicted by xcols2, e contains the residual or unpredicted or what is left to be predicted by zhedonia and zeudaimonia. So permute e and there should be no expected correlation between happiness and ypi.
    Ypi <- yhat + e # if iter=1 then ypi=y, otherwise permuted
    form <- formula(paste('Ypi',paste(xcols,collapse='+'),sep='~'))
    fitmv.pi <- lm(form, data=dt)
    # check when iter = 1 the obs data are analyzed correctly
    #form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
    #data.table(
    #  sep=coefficients(fitmv.pi)['zhedonia',],
    #  tog=coefficients(lm(form, data=dt))['zhedonia',])
    
    ss <- summary(fitmv.pi)
    for(j in 1:p){ # save t-values instead of coefficients
      coeffs.hed[iter,j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zhedonia', "Estimate"]
      coeffs.eud[iter,j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zeudaimonia', "Estimate"]
      t.hed[iter,j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zhedonia', "t value"]
      t.eud[iter,j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zeudaimonia', "t value"]
    }

    # permute rows
    rows <- sample(1:nrow(dt))
  }
  
  # matrix of regression coefficients for each gene
  b_table <- data.table(I=rep(1:niter,2), data=rep(c('obs',rep('perm',niter-1)),2), Type=rep(c('hedonic','eudaimonic'),each=niter), rbind(coeffs.hed,coeffs.eud))
  # matrix of regression t stats for each gene
  t_table <- data.table(I=rep(1:niter,2), data=rep(c('obs',rep('perm',niter-1)),2), Type=rep(c('hedonic','eudaimonic'),each=niter), rbind(t.hed,t.eud))
   
  fn1 <- paste('coeff.',fn,sep='')
  if(write_it==TRUE){write.table(b_table,fn1,quote=FALSE,sep='\t',row.names=FALSE)}
  fn2 <- paste('t.',fn,sep='')
  if(write_it==TRUE){write.table(t_table,fn2,quote=FALSE,sep='\t',row.names=FALSE)}
  

  return(NULL)
}

permutation_gls <- function(dt, xcols,ycols, zcols=c('zhedonia','zeudaimonia'), method='gls', niter=100, do_obs=TRUE, write_it=FALSE, fn){
  # hybrid of permutation_anderson and gls_with_correlated_error
  # basically use the permuted data to get a null distribution of b_hedonia and b_zeudaimonia
  
  #fn <- 'cole2_clean.txt'
  #dt <- read_file(fn,year=2015)
  #xcols <- get_xcols()
  #ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  
  # create matrix with scaled y
  #y.scale <- scale(dt[,.SD,.SDcols=ycols])
  Y <- scale(dt[,.SD,.SDcols=ycols])
  dt <- cbind(dt[,.SD,.SDcols=xcols],Y)
  dt[,subject:=factor(.I)]
  
  p <- length(ycols)
  xcols2 <- setdiff(xcols, zcols)
  # Get residuals from X (so excluding zhedonia and zeudaimonia)
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit.obs <- lm(form, data=dt)
  yhat <- predict(fit.obs) # Yhat = aX
  rows <- 1:nrow(dt) # observed on first iter and permuted after

  # coefficients (save t-value in addition to coefficients)
  b_matrix <- matrix(0,nrow=niter,ncol=3*length(zcols))
  #colnames(b_matrix) <- c('coeff.zhedonia','coeff.zeudaimonia','t.zhedonia','t.zeudaimonia')
  colnames(b_matrix) <- c(paste('coeff',zcols,sep='.'),paste('t',zcols,sep='.'),paste('p',zcols,sep='.'))
  
  #gee table added
  gee_perm_table <- data.table(NULL)
  samp <- 'obs'

  code <- sample(LETTERS,4,replace=TRUE)
  fn_full <- paste(fn,'.',paste(code,collapse=''),'.txt',sep='')
  fn_temp <- paste(fn,'.',paste(code,collapse=''),'.temp.txt',sep='')
  fn_gee <- paste(fn,'.',paste(code,collapse=''),'.gee.txt',sep='')
  
  for(iter in 1:niter){
    e <- residuals(fit.obs)[rows,] # permuted Yhat - aX
    Ypi <- scale(yhat + e)
    dt <- cbind(dt[,.SD,.SDcols=xcols],Ypi)
    dt[,subject:=factor(.I)] # need to recalc since dt is re-created
    # check apply(dt[,.SD,.SDcols=ycols],2,sd)
    
    # wide to long
    dtlong <- melt(dt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
    dtlong[,gene:=factor(gene)]
    dtlong <- orderBy(~subject + gene, dtlong)
    
    # this replicates the analysis of Fredrickson et al.
    form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
    if(method=='gls'){
      fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
    }
    if(method=='lme'){
      fit1 <- lme(form, random = ~1|subject, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(maxIter=100, msMaxIter = 500, tolerance=1e-6, msVerbose = FALSE)) # tolerance=1e-6 default
      #summary(fit1)$tTable[zcols,]
    }
    if(method=='lmer'){
      form3 <- formula(paste('expression~',paste(xcols,collapse='+'),'+(1|subject) + (1|gene)',sep=''))
      fit1 <- lmer(form3, data=dtlong, REML=FALSE) # tolerance=1e-6 default
      #summary(fit1)$tTable[zcols,]
    }
    if(method=='gee'){
      fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='exchangeable', std.err="san.se")
      estimate <- summary(fit.geeglm)$coefficients[zcols,'Estimate']
      names(estimate) <- zcols
      gee_perm_table <- rbind(gee_perm_table,data.table(samp=samp,t(estimate)))
    }
    b_matrix[iter,] <- c(
      summary(fit1)$tTable[zcols,'Value'],
      summary(fit1)$tTable[zcols,'t-value'],
      summary(fit1)$tTable[zcols,'p-value']
    )

    # permute rows for next iteration
    rows <- sample(1:nrow(dt))
    samp <- 'perm'
    
    # partial writing so I can check results
    if(method!='gee'){write.table(b_matrix,fn_temp,quote=FALSE,row.names = FALSE)}
  }
  
  b_matrix <- data.table(permutation=c('obs',rep('perm',niter-1)),b_matrix)
  # re-write with permutation column
  if(method!='gee'){write.table(b_matrix,fn_full,quote=FALSE,row.names = FALSE)}
  if(method=='gee'){write.table(gee_perm_table,fn_gee,quote=FALSE,row.names=FALSE,sep='\t')}
  return(NULL)
}

bootstrap_obrien <- function(dt,fn, niter=200,initer=1000){
  #initer is the inner bootstrap iteration to get the error on the correlation among the beta coefficients
  # iter is the outer iteration for bootstrap
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  xcols2 <- setdiff(xcols, zcols)
  N <- nrow(dt)
  p <- length(ycols)
  J <- matrix(1,nrow=p,ncol=1)
  tJ <- t(J)
  Bhed <- matrix(0,nrow=initer,ncol=p)
  Beud <- matrix(0,nrow=initer,ncol=p)
  Bdelta <- matrix(0,nrow=initer,ncol=p)
  Ttable <- data.table(NULL)
  Y <- dt[,.SD, .SDcols=ycols] # Y is rescaled in the loop
  rows <- 1:N # observed on first iter and permuted after
  for(iter in 1:niter){
    # don't need to rescale these at this point because they will be rescaled in the inner loop
    Ysamp <- Y[rows,]
    dt.samp <- dt[rows,]
    # beta for each iteration is the first row in the inner loop
    irows <- 1:nrow(dt) # observed on first iter and permuted after
    # now resample dt and Ypi many time
    for(iiter in 1:initer){
      Y.insamp <- scale(Ysamp[irows,])
      dt.insamp <- dt.samp[irows,]
      dt.insamp[, zhedonia:=scale(zhedonia)]
      dt.insamp[, zeudaimonia:=scale(zeudaimonia)]
      smoke <- dt.insamp[,smoke]
      if(length(which(smoke==1))<2){txcols <- setdiff(xcols,'smoke')}else{txcols<- xcols}
      form <- formula(paste('Y.insamp',paste(txcols,collapse='+'),sep='~'))
      fitmv.pi <- lm(form, data=dt.insamp)
      # these need to be the t-value not the coefficient
      Bhed[iiter,] <- coefficients(fitmv.pi)['zhedonia',]
      Beud[iiter,] <- coefficients(fitmv.pi)['zeudaimonia',]
      Bdelta[iiter,] <- coefficients(fitmv.pi)['zhedonia',] - coefficients(fitmv.pi)['zeudaimonia',]
      irows <- sample(1:N, replace=TRUE)
    }
    R <- cor(Bhed)
    num <- tJ%*%Bhed[1,]
    denom <- c(sqrt(tJ%*%R%*%J))
    Thed <- c(num/denom)
    R <- cor(Beud)
    num <- tJ%*%Beud[1,]
    denom <- c(sqrt(tJ%*%R%*%J))
    Teud <- c(num/denom)
    R <- cor(Bdelta)
    num <- tJ%*%Bdelta[1,]
    denom <- c(sqrt(tJ%*%R%*%J))
    Tdelta <- c(num/denom)
    Ttable <- rbind(Ttable,data.table(Thed=Thed,Teud=Teud,Tdelta=Tdelta))
    rows <- sample(1:N, replace=TRUE)
  }
  code <- sample(LETTERS,4,replace=TRUE)
  code <- 'boot'
  fn_full <- paste(fn,'.',paste(code,collapse=''),'.txt',sep='')
  write.table(Ttable,fn_full,quote=FALSE,sep='\t',row.names=FALSE)
}

permutation_obrien <- function(dt,fn, niter=200,initer=1000){
  # returns the null distribution for O'Brien's T
  #initer is the inner bootstrap iteration to get the error on the correlation among the beta coefficients
  # iter is the outer iteration for permutation
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  xcols2 <- setdiff(xcols, zcols)
  N <- nrow(dt)
  p <- length(ycols)
  J <- matrix(1,nrow=p,ncol=1)
  tJ <- t(J)
  Bhed <- matrix(0,nrow=initer,ncol=p)
  Beud <- matrix(0,nrow=initer,ncol=p)
  Bdelta <- matrix(0,nrow=initer,ncol=p)
  t.hed <- numeric(p)
  t.eud <- numeric(p)
  t.delta <- numeric(p)
  Ttable <- data.table(NULL)
  Y <- scale(dt[,.SD, .SDcols=ycols])
  dts <- copy(dt)
  dts[, zhedonia:=scale(zhedonia)]
  dts[, zeudaimonia:=scale(zeudaimonia)]
  # predicted and residuals for all X other than zcols
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit.obs <- lm(form, data=dts)
  yhat <- predict(fit.obs) # Yhat = aX
  rows <- 1:N # observed on first iter and permuted after
  for(iter in 1:niter){
    e <- residuals(fit.obs)[rows,] # permuted Yhat - aX
    # yhat is that predicted by xcols2, e contains the residual or unpredicted or what is left to be predicted by zhedonia and zeudaimonia. So permute e and there should be no expected correlation between happiness and ypi.
    Ypi <- yhat + e # if iter=1 then ypi=y, otherwise permuted
    irows <- 1:nrow(dt) # observed on first iter and permuted after
    # now resample dt and Ypi many time
    for(iiter in 1:initer){
      Yp.samp <- scale(Ypi[irows,])
      dt.samp <- dts[irows,]
      dt.samp[, zhedonia:=scale(zhedonia)]
      dt.samp[, zeudaimonia:=scale(zeudaimonia)]
      smoke <- dt.samp[,smoke]
      if(length(which(smoke==1))<2){txcols <- setdiff(xcols,'smoke')}else{txcols<- xcols}
      form <- formula(paste('Yp.samp',paste(txcols,collapse='+'),sep='~'))
      fitmv.pi <- lm(form, data=dt.samp)
      Bhed[iiter,] <- coefficients(fitmv.pi)['zhedonia',]
      Beud[iiter,] <- coefficients(fitmv.pi)['zeudaimonia',]
      Bdelta[iiter,] <- coefficients(fitmv.pi)['zhedonia',] - coefficients(fitmv.pi)['zeudaimonia',]
      if(iiter==1){
        ss <- summary(fitmv.pi)
        for(j in 1:p){ # save t-values instead of coefficients
          t.hed[j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zhedonia', "t value"]
          t.eud[j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zeudaimonia', "t value"]
          t.delta[j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zhedonia', "t value"] - ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zeudaimonia', "t value"]
        }
      }
      irows <- sample(1:N, replace=TRUE)
    }
    R <- cor(Bhed)
    num <- sum(t.hed)
    denom <- sqrt(sum(R)) # c(sqrt(tJ%*%R%*%J))
    Thed <- c(num/denom)
    R <- cor(Beud)
    num <- sum(t.eud)
    denom <- sqrt(sum(R)) # c(sqrt(tJ%*%R%*%J))
    Teud <- c(num/denom)
    R <- cor(Bdelta)
    num <- sum(t.delta)
    denom <- sqrt(sum(R)) # c(sqrt(tJ%*%R%*%J))
    Tdelta <- c(num/denom)
    Ttable <- rbind(Ttable,data.table(Thed=Thed,Teud=Teud,Tdelta=Tdelta))
    rows <- sample(1:N, replace=TRUE)
  }
  code <- sample(LETTERS,4,replace=TRUE)
  code <- 'perm'
  fn_full <- paste(fn,'.',paste(code,collapse=''),'.txt',sep='')
  write.table(Ttable,fn_full,quote=FALSE,sep='\t',row.names=FALSE)
}


roast_it <- function(dt){
  # Roast
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  Y <- scale(dt[, .SD, .SDcols=ycols])
  dt[,zhedonia:=scale(zhedonia)]
  dt[,zeudaimonia:=scale(zeudaimonia)]
  Yt <- t(Y)
  form <- formula(paste('~',paste(xcols,collapse='+'),sep=''))
  design <- model.matrix(form, data=dt)
  cont.matrix <- makeContrasts(delta="zhedonia-zeudaimonia",levels=design)
  hedonia <- which(colnames(design)=='zhedonia')
  eudaimonia <- which(colnames(design)=='zeudaimonia')
  phed <- roast(y=Yt,design=design,contrast=hedonia, nrot=1999)
  peud <- roast(y=Yt,design=design,contrast=eudaimonia, nrot=1999)
  pdelta <- roast(y=Yt,design=design,contrast=cont.matrix, nrot=1999)  
  
  
  p.roast <- c(zhedonia=phed$p.value['UpOrDown','P.Value'], zeudaimonia=peud$p.value['UpOrDown','P.Value'], delta=pdelta$p.value['UpOrDown','P.Value'])
  return(p.roast)
  
}

convert_gls_temp_file <- function(){
  # some of the gls permutation runs failed to converge and the function ended before completion so this function is converting the temp file to the final file
  fn <- 'FRED13.permutation.gls.KJBH.temp.txt'
  fn_out <- 'FRED13.permutation.gls.KJBH.txt'
  do_obs <- TRUE
  b_matrix <- read.table(fn, header=TRUE)
  niter <- nrow(b_matrix)
  if(do_obs==TRUE){
    b_matrix <- data.table(permutation=c('obs',rep('perm',niter-1)),b_matrix)
  }else{
    b_matrix <- data.table(permutation=rep('perm',niter),b_matrix)
  }
  # re-write with permutation column
  write.table(b_matrix,fn_out,quote=FALSE,row.names = FALSE)
}

bootstrap_test <- function(dt,xcols,ycols,niter=1000, write_it=FALSE, fn){
# resamples dt and computes the standardized beta of zhedonia and zeudaimonia
# first row is observed data
# dt is the data.table of X regressors and Y responses
# res is the resulting coefficients of the regression
# xcols are the regressors
# ycols are the responses
# fn is the file name to write to
  p <- length(ycols)
  hed_matrix <- matrix(0,nrow=niter,ncol=p) # matrix of coefficients
  eud_matrix <- matrix(0,nrow=niter,ncol=p) # matrix of coefficients
  colnames(hed_matrix) <- ycols
  colnames(eud_matrix) <- ycols

  rows <- 1:nrow(dt) # observed on first iter and permuted after
  for(iter in 1:niter){
    Y <- scale(as.matrix(dt[rows, .SD, .SDcols=ycols]))
    dts <- dt[rows,]
    dts[,zhedonia:=scale(zhedonia)]
    dts[,zeudaimonia:=scale(zeudaimonia)]
    # if too few cases with smoke=1 then drop it
    smoke <- dts[,smoke]
    if(length(which(smoke==1))<2){txcols <- setdiff(xcols,'smoke')}else{txcols<- xcols}
    form <- formula(paste('Y',paste(txcols,collapse='+'),sep='~'))
    fitmv <- lm(form, data=dts[,.SD,.SDcols=txcols])
    hed_matrix[iter,] <- coefficients(fitmv)['zhedonia',]
    eud_matrix[iter,] <- coefficients(fitmv)['zeudaimonia',]
    rows <- sample(1:nrow(dt),replace=TRUE)
  }
  b_matrix <- data.table(I=rep(1:niter,2), data=rep(c('obs',rep('resamp',niter-1)),2), Type=rep(c('hedonic','eudaimonic'),each=niter), rbind(hed_matrix,eud_matrix))
  
  if(write_it==TRUE){write.table(b_matrix,fn,quote=FALSE,sep='\t',row.names=FALSE)}
  return(NULL)
}

wald_test <- function(x,se){
  pchisq(x^2/se^2, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
}

gee_delta <- function(gee_res){
  # given the two estimates in est and the SEs of the estimates in se,
  # this computes the SE, Wald, and p-value of the difference in the estimates
  delta <- gee_res["zhedonia",'Estimate'] - gee_res["zeudaimonia",'Estimate']
  delta.se <- sqrt(sum(gee_res[,'Std.err']^2))
  delta.wald <- delta^2/delta.se^2
  delta.p <- pchisq(delta.wald, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  gee_res <- rbind(gee_res,c(delta,delta.se,delta.wald,delta.p))
  row.names(gee_res)[3] <- 'delta'
  return(gee_res)
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

glh.p.values <- function(glh_mod){
  tol <- 1e-12
  p <- nrow(glh_mod$SSPH)
  lambda.h <- eigen(glh_mod$SSPH)$values
  lambda.h <- lambda.h[lambda.h>tol]
  h <- length(lambda.h)
  lambda.e <- eigen(glh_mod$SSPE)$values
  lambda.e <- lambda.e[lambda.e>tol]
  e <- length(lambda.e)
  lambda.he <- eigen(glh_mod$SSPE+glh_mod$SSPE)$values
  lambda.he <- lambda.he[lambda.he>tol]
  he <- length(lambda.he)
  wilks <- det(glh_mod$SSPE)/det(glh_mod$SSPE + glh_mod$SSPH) 
  pillai <- sum(diag(glh_mod$SSPH%*%solve(glh_mod$SSPE + glh_mod$SSPH)))
  s <- min(p,h)
  m <- (abs(p-h)-1)/2
  n <- (e-p-1)/2
  ndf <- s*(2*m+s+1)
  ddf <- s*(2*n+s+1)
  pillaiF <- (2*n+s+1)*pillai/((2*m+s+1)*(s-pillai))
  pillaiP <- pf(pillaiF,ndf,ddf,lower.tail=FALSE)
  return(pillaiP)
}

t.test.p.value <- function(x){
  # usage:  t <- apply(res, 2, t.test.p.value)
  return(t.test(x)$p.value)
}

t.test.statistic <- function(x){
  # usage:  t <- apply(res, 2, t.test.statistic)
  return(t.test(x)$statistic)
}

permutation.p.value <- function(x){
  # obs statistic must be in first cell
  # for two sided, x must be abs(x)
  # usage:  t <- apply(res, 2, permutation.p.value)
  return(length(which(x >= x[1]))/length(x))
}

permutation.gls.p.value <- function(res, statistic='t',do_ci=FALSE){
  # res has columns for coefficient and t stat
  # read permuation.gls moved permutation=obs to 1st row and deleted additional permuation=obs
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
    p.value <- apply(abs(t.value),2,permutation.p.value)
    
    # bootstrap CIs of p-value
    if(do_ci==TRUE){
      t.obs <- t.value[1,]
      niter <- 2000
      bootp <- matrix(0,nrow=niter,ncol=3)
      for(iter in 1:niter){
        rows <- sample(2:nrow(res),(nrow(res)-1),replace=TRUE)
        inc <- which(substr(colnames(res),1,1)=='t')
        t.value <- res[rows,.SD,.SDcols=colnames(res)[inc]]
        inc <- which(substr(colnames(res),1,1)=='c')
        coeff.value <- res[rows,.SD,.SDcols=colnames(res)[inc]]
        se.value <- coeff.value/t.value
        delta.se <- sqrt(apply(se.value^2,1,sum))
        delta.t <- (coeff.value[,coeff.zhedonia] - coeff.value[,coeff.zeudaimonia])/delta.se
        t.value <- rbind(t.obs,cbind(t.value, delta=delta.t))
        bootp[iter,] <- apply(abs(t.value),2,permutation.p.value)
      }
      ci <- apply(bootp,2,quantile, probs=c(0.025,0.975))
      p.value <- rbind(p.value,ci)
      p.value <- data.table(stat=row.names(p.value),p.value)
    }
  }
  if(statistic=='c'){
    inc <- which(substr(colnames(res),1,1)=='c')
    p.value <- apply(abs(res[,.SD,.SDcols=colnames(res)[inc]]),2,permutation.p.value)
    p.value <- c(p.value, delta=permutation.p.value(abs(res[,coeff.zhedonia]-res[,coeff.zeudaimonia])))
  }
  return(p.value)
}


type.1.error <- function(x, alpha=0.05){
  # obs statistic must be in first cell
  # usage:  t <- apply(res, 2, type.1.error)
  return(length(which(x <= alpha))/length(x))
}

effect_size <- function(dt){
  # fit the UNstandardized Y
  # change in mean expression give 4sd change in hedonic score
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  Y <- as.matrix(dt[, .SD, .SDcols=ycols])
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv <- lm(form, data=dt[,.SD,.SDcols=xcols])
  coeffs <- data.table(zhedonia=coefficients(fitmv)['zhedonia',], zeudaimonia=coefficients(fitmv)['zeudaimonia',])
  obs2 <- (2^(4*apply(coeffs,2,mean)) - 1)*100
  return(obs2)
}

GLH <- function(dt, explore_it=FALSE){
  # General Linear Hypothesis to estimate p-value associated with
  # zhedonia=0
  # zeudaimonia=0
  # zhedonia = zeudaimonia
  # dt is the XY matrix of responses and regressors
  # returns stats for zhedonia, zeudaimonia, delta
 
  #
  if(explore_it == TRUE){
    dt <- copy(cole1)
  }
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  Y <- as.matrix(dt[,.SD,.SDcols=ycols])
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv <- lm(form, data=dt)
  GLH1 <- linearHypothesis(fitmv, "zhedonia = 0", test='Pillai') #
  GLH2 <- linearHypothesis(fitmv, "zeudaimonia = 0", test='Pillai')
  GLH3 <- linearHypothesis(fitmv, "zhedonia = zeudaimonia", test='Pillai')

  if(explore_it==TRUE){
    Xobs <- dt[,.SD,.SDcols=xcols]
    # convert factors to numeric
    Xobs[,male:=as.numeric(as.character(male))]
    Xobs[,white:=as.numeric(as.character(white))]
    Xobs[,smoke:=as.numeric(as.character(smoke))]
    Xobs <- as.matrix(Xobs)
    fit <- eigen(cov(Y))
    pc1 <- Y%*%fit$vectors[,1]
    summary(lm(pc1~Xobs))
    fitmvr <- lm(form, data=dt)
    GLH1r <- linearHypothesis(fitmvr, "zhedonia = 0")
    GLH2r <- linearHypothesis(fitmvr, "zeudaimonia = 0")
    GLH3r <- linearHypothesis(fitmvr, "zhedonia = zeudaimonia")
    # machs nichts
    
    #x on y - overall test makes no sense because of covariates but GLH would
    xcols1 <- c(ycols,setdiff(xcols,'zhedonia'))
    form <- formula(paste('zhedonia',paste(xcols,collapse='+'),sep='~'))
    fitmv <- lm(form, data=dt)
    summary(fitmv)
    
    # Y on X v. X on Y
    summary(lm(zhedonia~IL1A + IL1B, data=dt))
    Y <- as.matrix(dt[,.SD,.SDcols=c('IL1A', 'IL1B')])
    anova(lm(Y~zhedonia,data=dt))
    # so what if we took residuals of covariates then did manova?
    xcols1 <- setdiff(xcols,'zhedonia')
    Y <- as.matrix(dt[,.SD,.SDcols=ycols])
    form <- formula(paste('Y',paste(xcols1,collapse='+'),sep='~'))
    fitmv <- lm(form, data=dt)
    Yfit <- predict(fitmv)
    anova(lm(Yfit~zhedonia,data=dt))
    
  }
  
  return(list(GLH1=GLH1, GLH2=GLH2, GLH3=GLH3))
}

gls_type_I_error <- function(df,res){
  p.zhed <- 2*pt(abs(res[permutation=='perm',t.zhedonia]),df=df,lower.tail = FALSE)
  p.zeud <- 2*pt(abs(res[permutation=='perm',t.zeudaimonia]),df=df,lower.tail = FALSE)
  t1e.zhed <- type.1.error(p.zhed)
  t1e.zeud <- type.1.error(p.zeud)
  return(c(t1e.zhed=t1e.zhed, t1e.zeud=t1e.zeud))
}

naive_t_stats <- function(dt){
  # compute t-statistic of regulation as in Frederickson 2013
  # fit the UNstandardized Y
  # change in mean expression give 4sd change in hedonic score
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  Y <- as.matrix(dt[, .SD, .SDcols=ycols])
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv <- lm(form, data=dt[,.SD,.SDcols=xcols])
  coeffs <- data.table(zhedonia=coefficients(fitmv)['zhedonia',], zeudaimonia=coefficients(fitmv)['zeudaimonia',])
  # backtransform coeffs 
  beta <- apply(coeffs,2,mean)
  p <- nrow(coeffs)
  iters <- 200 # following the original
  yboot <- matrix(0,nrow=iters,ncol=2)
  rows <- 1:p
  for(iter in 1:iters){
 #   yboot[iter,] <- apply(coeffs[rows,],2,mean)
    yboot[iter,] <- (2^(4*apply(coeffs[rows,],2,mean)) - 1)*100
    rows <- sample(1:p,replace=TRUE)
  }
  colnames(yboot) <- colnames(coeffs)
  yboot <- data.table(yboot)
  yboot[,delta:=zhedonia-zeudaimonia]
  estimate <- as.numeric(yboot[1,])
  se <- apply(yboot,2,sd)
  se['delta'] <- sqrt(se['zhedonia']^2 + se['zeudaimonia']^2) # same as before
  t <- estimate/se
  p.value <- 2*pt(abs(t),df=(p-1),lower.tail = FALSE)
  t.table <- data.table(stat=c('estimate','se','t','p.value'),rbind(estimate,se,t,p.value))
  
  return(t.table)
  
}

smart_t_stats <- function(cc){
  # cc must be the contrast coefficients and not raw regression coefficients
  if('IL6' %in% colnames(cc)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  cc[,ME:=apply(.SD,1,mean),.SDcols=ycols, by=Type]
  # create new table with zhedonia and zeudaimonia mean coefficients in columns
  beta <- data.table(zhedonia=cc[Type=='hedonic',ME], zeudaimonia=cc[Type=='eudaimonic',ME])
  beta[, delta:=zhedonia-zeudaimonia]
  means <- unlist(beta[1,])
  se <- apply(beta,2,sd)
  t <- abs(means/se)
  m <- length(ycols)
  prob <- 2*pt(t,df=(m-1),lower.tail = FALSE)
  return(prob)
}


permutation_t <- function(cc,statistic='t'){
  # cc must be the contrast coefficients and not raw regression coefficients
  if('IL6' %in% colnames(cc)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  p <- length(ycols)
  cc[,ME:=apply(.SD,1,mean),.SDcols=ycols, by=Type]
  cc[,SD:=apply(.SD,1,sd),.SDcols=ycols, by=Type]
  cc[,SE:=SD/sqrt(p), by=Type]
  cc[,t:=ME/SE, by=Type]
  
  delta <- cc[Type=='hedonic',ME] - cc[Type=='eudaimonic',ME]
  
  if(statistic=='mean'){
    theta <- data.table(
      zhedonia=cc[Type=='hedonic',ME],
      zeudaimonia=cc[Type=='eudaimonic',ME],
      delta=delta
    )
  }
  if(statistic=='t'){
    theta <- data.table(
      zhedonia=cc[Type=='hedonic',t],
      zeudaimonia=cc[Type=='eudaimonic',t],
      delta=delta/sqrt(cc[Type=='hedonic',SD]^2/p + cc[Type=='eudaimonic',SD]^2/p)
      )
  }
  
  ols.perm.p <- apply(abs(theta),2,permutation.p.value)
  return(ols.perm.p)
}

permutation_lambda <- function(res){
  # res is the table of Wilk's lambda or Pillai's trace for each iteration
  niter <- nrow(res)
  p_matrix <- data.table(p.hed=NA,p.eud=NA,p.delta=NA)
  
  # is hedonic up regulated? This is two-sided despite the question
  zhedonia <- length(which(res$zhedonia.0 >= res$zhedonia.0[1]))/niter # permutation p
  # is eudamonic down regulated?
  zeudaimonia <- length(which(res$zeudaimonia.0 >= res$zeudaimonia.0[1]))/niter # permutation p
  # is there a difference in regulation?
  delta <- length(which(res$delta.0 >= res$delta.0[1]))/niter # permutation p
  res <- c(zhedonia,zeudaimonia,delta)
  names(res) <- c('zhedonia','zeudaimonia','delta')
  return(res)
}


bootstrap_stats <- function(cc){
  # cc must be the contrast coefficients and not raw regression coefficients
  if('IL6' %in% colnames(cc)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  cc[,ME:=apply(.SD,1,mean),.SDcols=ycols, by=Type]
  
  b.hed <- cc[Type=='hedonic', ME]
  b.eud <- cc[Type=='eudaimonic', ME]
  
  # is hedonic up regulated? This is two-sided despite the question
  hed.CI <- quantile(b.hed,c(0.0275,0.975))
  eud.CI <- quantile(b.eud,c(0.0275,0.975))
  delta.CI <- quantile((b.hed-b.eud),c(0.0275,0.975))
  ci_table <- data.table(label=c('hedonic','eudaimonic','delta'),rbind(hed.CI,eud.CI,delta.CI))
  return(ci_table)
}

correlation_p.value <- function(dt1, dt2){
  # dt must be the contrast coefficients and not raw regression coefficients
  ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  r <- numeric(nrow(dt1))
  for(i in 1:nrow(dt1)){
    y1 <- as.numeric(dt1[i,.SD, .SDcols=ycols])
    y2 <- as.numeric(dt2[i,.SD, .SDcols=ycols])
    y1s <- y1/sqrt(sum(y1^2))
    y2s <- y2/sqrt(sum(y2^2))
    # sum(y1s^2) # check
    # sum(y2s^2) # check
    r[i] <- as.numeric(t(y1s)%*%y2s)
  }
  niter <- nrow(dt1)/2
  res.table <- data.table(I=dt1[,I],Type=dt1[,Type], r=r)
  
  
  # check below
  #r.hed <- as.numeric(res.table[Type=='hedonic',r])
  #hed.p.value <- length(which(abs(r.hed) >= abs(r.hed[1])))/niter
  #r.eud <- as.numeric(res.table[Type=='eudaimonic',r])
  #eud.p.value <- length(which(abs(r.eud) >= abs(r.eud[1])))/niter
  return(res.table[, .(p.value=length(which(abs(r) >= abs(r[I==1])))/niter), by=Type])
}

explore_random_response <- function(){
  cole1 <- data.table(read.table('cole1_clean.txt',header=TRUE,sep='\t'))
  cole1[, male:=factor(male)]
  cole1[, white:=factor(white)]
  cole1[, smoke:=factor(smoke)]
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year=2013),antibody_genes(),ifn_genes())
  Y <- as.matrix(contrast_coefficients(cole1[,.SD,.SDcols=ycols]))
  
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv.obs <- lm(form, data=cole1)
  res.obs <- t(coefficients(fitmv.obs)[c('zhedonia','zeudaimonia'),])
  p.obs <- apply(res.obs, 2, t.test.p.value)
  t.obs <- apply(res.obs, 2, t.test.statistic)
  mu.obs <- apply(res.obs,2,mean)
  sd.obs <- apply(res.obs,2,sd)
  
  S <- cov(Y)
  niter <- 1000
  
  # simulation 1 - simply sampling coefficients with same sd as observed.
  p1 <- matrix(0,nrow=niter,ncol=2)
  t1 <- matrix(0,nrow=niter,ncol=2)
  mu1 <- matrix(0,nrow=niter,ncol=2)
  colnames(t) <- names(t.obs)
  colnames(mu) <- names(t.obs)
  t1[1,] <- t.obs
  mu1[1,] <- mu.obs
  p1[1,] <- p.obs
  res <- t(coefficients(fitmv.obs)[c('zhedonia','zeudaimonia'),]) # just creating this
  for(iter in 2:niter){
    res[,1] <- rnorm(n, sd=sd(res.obs[,1]))
    res[,2] <- rnorm(n, sd=sd(res.obs[,2]))
    p1[iter,] <- apply(res, 2, t.test.p.value)
    t1[iter,] <- apply(res, 2, t.test.statistic)
    mu1[iter,] <- apply(res, 2, mean)
  }
  apply(res.obs,2,sd)/sqrt(ncol(Y)) # se of mean coefficients from observed data
  apply(mu1,2,sd) # sd of mean from observed data, if the coefficients weren't estimated these should be the same as above
  apply(abs(t1),2,permutation.p.value)
  apply(abs(mu1),2,permutation.p.value)
  # result sim 1: 0.084 0.005 - these are the ~ naive t-test p.values
  # result sim 2: 0.4059      0.0876, quite common to get mean larger for hedonic
  apply(p1[2:niter,],2,type.1.error)
  
  
  # simulation 2 - no effect but coeff are estimated with Ssim=Diag(Sobs)
  niter <- 10^3
  p2 <- matrix(0,nrow=niter,ncol=2)
  t2 <- matrix(0,nrow=niter,ncol=2)
  mu2 <- matrix(0,nrow=niter,ncol=2)
  sd2 <- matrix(0,nrow=niter,ncol=2)
  colnames(t2) <- names(t.obs)
  colnames(mu2) <- names(t.obs)
  t2[1,] <- t.obs
  mu2[1,] <- mu.obs
  p2[1,] <- p.obs
  sd2[1,] <- sd.obs
  for(iter in 2:niter){
    Yr <- rmvnorm(nrow(Y), sigma=diag(diag(S)))
    colnames(Yr) <- colnames(Y)
    form <- formula(paste('Yr',paste(xcols,collapse='+'),sep='~'))
    fitmv <- lm(form, data=cole1)
    res <- t(coefficients(fitmv)[c('zhedonia','zeudaimonia'),])
    p2[iter,] <- apply(res, 2, t.test.p.value)
    t2[iter,] <- apply(res, 2, t.test.statistic)
    mu2[iter,] <- apply(res, 2, mean)
    sd2[iter,] <- apply(res, 2, sd)
  }
  apply(res.obs,2,sd)/sqrt(ncol(Y)) # se of mean coefficients from observed data
  apply(mu2[2:niter,],2,sd) # sd of mean from simulated data, if the coefficients weren't estimated these should be the same as above
  sd.obs
  apply(sd2[2:niter,],2,mean)
  apply(sd2[2:niter,],2,mean)/sqrt(ncol(Y))
  # result: true se or mean is alomst 2x estimated.
  apply(abs(t2),2,permutation.p.value)
  apply(abs(mu2),2,permutation.p.value)
  # result: 0.4059      0.0876, quite common to get mean larger for hedonic
  apply(p2[2:niter,],2,type.1.error)
  cbind(mu2,p1)[1:15,]
  
  # some exploration
  n <- p
  b <- rnorm(n, sd=sd(res.obs[,1])) # randomly sample regression coefficients"
  b.est <- b + rnorm(n, sd=mean(se.hed)) # add estimation error, se.hed grabbed from coeff_vector()
  sd(b.est) # = sqrt(2)*sd
  sd(res.obs[,1])/sqrt(p) # expected se of mean of coefficients
  sd(b.est)/sqrt(p) # se of mean given estimation error
  apply(mu,2,sd)[1] # compare
  
  # now repeat but simulate using full S and not just diag.
  p3 <- matrix(0,nrow=niter,ncol=2)
  t3 <- matrix(0,nrow=niter,ncol=2)
  mu3 <- matrix(0,nrow=niter,ncol=2)
  sd3 <- matrix(0,nrow=niter,ncol=2)
  colnames(t3) <- names(t.obs)
  colnames(mu3) <- names(t.obs)
  p3[1,] <- p.obs
  t3[1,] <- t.obs
  mu3[1,] <- mu.obs
  sd3[1,] <- sd.obs
  for(iter in 2:niter){
    Yr <- rmvnorm(nrow(Y), sigma=S)
    colnames(Yr) <- colnames(Y)
    form <- formula(paste('Yr',paste(xcols,collapse='+'),sep='~'))
    fitmv <- lm(form, data=cole1)
    res <- t(coefficients(fitmv)[c('zhedonia','zeudaimonia'),])
    p3[iter,] <- apply(res, 2, t.test.p.value)
    t3[iter,] <- apply(res, 2, t.test.statistic)
    mu3[iter,] <- apply(res, 2, mean)
    sd3[iter,] <- apply(res, 2, sd)
  }
  sd.obs/sqrt(ncol(Y)) # se of mean coefficients from observed data
  apply(mu3[2:niter,],2,sd) # sd of mean from observed data, if the coefficients weren't estimated these should be the same as above
  sd.obs
  apply(sd3[2:niter,],2,mean)
  apply(sd3[2:niter,],2,mean)/sqrt(ncol(Y))
  # result actual se about 6x estimated se
  apply(abs(t3),2,permutation.p.value)
  apply(abs(mu3),2,permutation.p.value)
  # result sim 1: 0.084     0.005 - these are the ~ naive t-test p.values
  # result sim 2: 0.4059    0.0876, quite common to get mean larger for hedonic
  # result sim 3: 0.784     0.509
  apply(p3[2:niter,],2,type.1.error)
  
  sd.obs
  apply(sd2[2:niter,],2,mean)
  apply(sd3[2:niter,],2,mean)
  apply(sd2[2:niter,],2,min)
  apply(sd3[2:niter,],2,min)

  sd.obs/sqrt(ncol(Y))
  apply(sd2[2:niter,],2,mean)/sqrt(ncol(Y))
  apply(sd3[2:niter,],2,mean)/sqrt(ncol(Y))
  
  hist(mu2[1,])
  hist(mu3[1,])
  
  # *** bias in mu + for hed and - for eud
  apply(mu2,2,mean)
  apply(mu3,2,mean)
  apply(mu3,2,sd)/sqrt(niter)*2
  # TO DO - cannot have 
  
  # does this bias depend on covariates
  p4 <- matrix(0,nrow=niter,ncol=2)
  t4 <- matrix(0,nrow=niter,ncol=2)
  mu4 <- matrix(0,nrow=niter,ncol=2)
  sd4 <- matrix(0,nrow=niter,ncol=2)
  colnames(t4) <- names(t.obs)
  colnames(mu4) <- names(t.obs)
  p4[1,] <- p.obs
  t4[1,] <- t.obs
  mu4[1,] <- mu.obs
  sd4[1,] <- sd.obs
  for(iter in 2:niter){
    Yr <- rmvnorm(nrow(Y), sigma=S)
    colnames(Yr) <- colnames(Y)
    form <- formula('Yr ~ zhedonia + zeudaimonia')
    fitmv <- lm(form, data=cole1)
    res <- t(coefficients(fitmv)[c('zhedonia','zeudaimonia'),])
    p4[iter,] <- apply(res, 2, t.test.p.value)
    t4[iter,] <- apply(res, 2, t.test.statistic)
    mu4[iter,] <- apply(res, 2, mean)
    sd4[iter,] <- apply(res, 2, sd)
  }
  apply(mu2,2,mean)
  apply(mu3,2,mean)
  apply(mu4,2,mean) # still + and - bias but not as big
  
  # Random hed and eud but with same correlation structure as obs.
  psy <- cole1[,.SD, .SDcols=c('zhedonia','zeudaimonia')]
  Spsy <- cov(psy)
  #Spsy[2,1] <- -Spsy[2,1] # this just reverses sign of cor(mu5)
  #Spsy[1,2] <- -Spsy[1,2]
  covcols <- setdif(xcols,c('zhedonia','zeudaimonia'))
  pseudocole1 <- copy(cole1)
  p5 <- matrix(0,nrow=niter,ncol=2)
  t5 <- matrix(0,nrow=niter,ncol=2)
  mu5 <- matrix(0,nrow=niter,ncol=2)
  sd5 <- matrix(0,nrow=niter,ncol=2)
  colnames(t5) <- names(t.obs)
  colnames(mu5) <- names(t.obs)
  p5[1,] <- p.obs
  t5[1,] <- t.obs
  mu5[1,] <- mu.obs
  sd5[1,] <- sd.obs
  for(iter in 2:niter){
    Xpsy <- rmvnorm(nrow(Y), sigma=Spsy)
    colnames(Xpsy) <- c('zhedonia','zeudaimonia')
    pseudocole1[, zhedonia:=Xpsy[,'zhedonia']]
    pseudocole1[, zeudaimonia:=Xpsy[,'zeudaimonia']]
    form <- formula(paste('Yr',paste(xcols,collapse='+'),sep='~'))
    fitmv <- lm(form, data=pseudocole1)
    res <- t(coefficients(fitmv)[c('zhedonia','zeudaimonia'),])
    p5[iter,] <- apply(res, 2, t.test.p.value)
    t5[iter,] <- apply(res, 2, t.test.statistic)
    mu5[iter,] <- apply(res, 2, mean)
    sd5[iter,] <- apply(res, 2, sd)
  }
  apply(mu3,2,mean)
  apply(mu5,2,mean) 
  cor(Xpsy)[2,1]
  cor(mu5)[2,1]
  qplot(x=zhedonia, y=zeudaimonia, data=data.table(mu5))
  # **** fascinating. Basically it's 50-50, if zhedonia is + or - but if + then zeudaimonia is negative. This is what has been inchoate in my mind since the start, that the negative correlation is a geometric artifact by having both high correlation Y in the model.Sort of like Pearson-Aitchison but the explanation isn't immediately obvious. If one partial regression somewhat like fitting the regression on residuals. Obviously there will be a small correlation with the residuals since much of the correlation is shared and it has to be in the opposite direction. so...
  
  # One happiness score without other in the model
  psy <- cole1[,.SD, .SDcols=c('zhedonia','zeudaimonia')]
  Spsy <- cov(psy)
  #Spsy[2,1] <- -Spsy[2,1] # this just reverses sign of cor(mu5)
  #Spsy[1,2] <- -Spsy[1,2]
  covcols <- setdif(xcols,c('zhedonia','zeudaimonia'))
  pseudocole1 <- copy(cole1)
  p6 <- matrix(0,nrow=niter,ncol=2)
  t6 <- matrix(0,nrow=niter,ncol=2)
  mu6 <- matrix(0,nrow=niter,ncol=2)
  sd6 <- matrix(0,nrow=niter,ncol=2)
  colnames(t6) <- names(t.obs)
  colnames(mu6) <- names(t.obs)
  p6[1,] <- p.obs
  t6[1,] <- t.obs
  mu6[1,] <- mu.obs
  sd6[1,] <- sd.obs
  for(iter in 2:niter){
    Xpsy <- rmvnorm(nrow(Y), sigma=Spsy)
    colnames(Xpsy) <- c('zhedonia','zeudaimonia')
    pseudocole1[, zhedonia:=Xpsy[,'zhedonia']]
    pseudocole1[, zeudaimonia:=Xpsy[,'zeudaimonia']]
    xcols1 <- setdiff(xcols,'zeudaimonia')
    form <- formula(paste('Yr',paste(xcols1,collapse='+'),sep='~'))
    fitmv <- lm(form, data=pseudocole1)
    res[,'zhedonia'] <- t(coefficients(fitmv)[c('zhedonia'),])
    xcols2 <- setdiff(xcols,'zhedonia')
    form <- formula(paste('Yr',paste(xcols2,collapse='+'),sep='~'))
    fitmv <- lm(form, data=pseudocole1)
    res[,'zeudaimonia'] <- t(coefficients(fitmv)[c('zeudaimonia'),])
    p6[iter,] <- apply(res, 2, t.test.p.value)
    t6[iter,] <- apply(res, 2, t.test.statistic)
    mu6[iter,] <- apply(res, 2, mean)
    sd6[iter,] <- apply(res, 2, sd)
  }
  qplot(x=zhedonia, y=zeudaimonia, data=data.table(mu6)) 
  
  # obvious!
  
  # simulate
  
  r <- 0.8
  b <- sqrt(r)
  n <- nrow(cole1)
  p <- length(ycols)
  z <- rnorm(n)
  x1 <- b*z + sqrt(1-r)*rnorm(n)
  x2 <- b*z + sqrt(1-r)*rnorm(n)
  Y <- matrix(rnorm(n*p),nrow=n)
  fitmv <- lm(Y~x1+x2)
  coeffs <- t(coefficients(fitmv)[c('x1','x2'),])
  qplot(x=x1, y=x2, data=data.table(coeffs)) # bam!
  apply(coeffs,2,mean)
  # *** because of negative correlation of coefficients, the mean of the coefficeints will usually be of opposite sign
  
  coeffs2 <- copy(coeffs)
  Yres <- residuals(lm(Y~x2))
  coeffs2[,'x1' ] <- t(coefficients(lm(Yres~x1))['x1',])
  Yres <- residuals(lm(Y~x1))
  coeffs2[,'x2'] <- t(coefficients(lm(Yres~x2))['x2',])
  cbind(coeffs,coeffs2)
  
  # visualize
  j <- 20
  Yres <- residuals(lm(Y~x2))
  coeffs[j,]
  
  qplot(x=x2,y=Y[,j])
  cor(x2,Y[,j])
  qplot(x=x1,y=Yres[,j])
  cor(x1,Yres[,j])
}

gls_with_correlated_error <- function(dt,year=2015, fit_lme=FALSE){
  #if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  dtlong <- get_dtlong(dt,year=year,center_Y=TRUE) # returns contrast coefficients
  dtlong <- orderBy(~subject + gene, dtlong)
  
  # this replicates the analysis of Fredrickson et al.
  # subject is the grouping factor and this is identified in the correlation statement
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  if(fit_lme==FALSE){
    fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  }
  if(fit_lme==TRUE){
    fit1 <- lme(form, random = ~1|subject, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(maxIter=100, msMaxIter = 500, tolerance=1e-6, msVerbose = FALSE)) # tolerance=1e-6 default
  }

  return(fit1)
}

gee_with_correlated_error <- function(dt){
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  dtlong <- get_dtlong(dt,year=year,center_Y=TRUE) # returns contrast coefficients
  dtlong <- orderBy(~subject + gene, dtlong)
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  # check out interaction
  # form <- formula(paste('expression~',paste(c('gene',xcols,paste(zcols,collapse=':')),collapse='+'),sep=''))
  fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='exchangeable', std.err="san.se")
  #fit.gls <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  sum.lm <- summary(lm(form,data=dtlong))$coefficients[zcols,]
  sum.gee <- summary(fit.geeglm)$coefficients[zcols,]
  #sum.gls <- summary(fit.gls)$tTable[zcols,]
  #res <- rbind(sum.lm, sum.gee)
  return(sum.gee)
}

explore_gls_with_correlated_error <- function(){
  # note that this replicates Fredrickson 2015
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  
  # create matrix with scaled y
  #y.scale <- scale(dt[,.SD,.SDcols=ycols])
  Y <- scale(contrast_coefficients(dt[,.SD,.SDcols=ycols]))
  #Y <- scale(contrast_coefficients(dt[,.SD,.SDcols=ycols]),center=FALSE)
  dts <- cbind(dt[,.SD,.SDcols=xcols], Y)
  dts[,subject:=factor(.I)]
  
  # wide to long
  dtlong <- melt(dts,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  
  # this replicates the analysis of Fredrickson et al.
  # subject is the grouping factor and this is identified in the correlation statement
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit1)$tTable[zcols,]
  #              Value Std.Error    t-value p-value
  # zhedonia     0.08575  0.122458   0.700267  0.4838
  # zeudaimonia -0.51069  0.125713  -4.062329  0.0000
  anova(fit1)
  Anova(fit1)
  
  # gene is a ?? effect
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  fit.gls.2 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ gene | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit.gls.2)$tTable
  
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  fit.gls.3 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ gene | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit.gls.3)$tTable
  
  # unrestricted
  fit.gls.4 <- gls(form, data=dtlong, method='ML', correlation=corSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, tolerance=1e-6, msVerbose = FALSE))
  
  #lme
  form <- formula('expression ~ gene + male + age + white + bmi + alcohol + smoke + 
    illness + cd3d + cd3e + cd4 + cd8a + fcgr3a + cd19 + ncam1 + 
                  cd14 + zhedonia + zeudaimonia')
  fit5 <- lme(form,data=dtlong,random = ~1|subject)
  summary(fit5)$tTable[c('zhedonia','zeudaimonia'),]
  form <- formula('expression ~ male + age + white + bmi + alcohol + smoke + 
    illness + cd3d + cd3e + cd4 + cd8a + fcgr3a + cd19 + ncam1 + 
                  cd14 + zhedonia + zeudaimonia')
  fit5 <- lme(form,data=dtlong,random = ~1|subject)
  summary(fit5)$tTable[c('zhedonia','zeudaimonia'),]
  fit5 <- gls(form,data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject))
  summary(fit5)$tTable[c('zhedonia','zeudaimonia'),]
  fit5.gls <- gls(form,data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject),weights=varIdent(form = ~1|gene))
  summary(fit5.gls)$tTable[c('zhedonia','zeudaimonia'),]
  #********** key key key
  fit5.lme <- lme(form,data=dtlong,random = ~1|subject,cor=corCompSymm(form=~1|subject),weights=varIdent(form = ~1|gene), control=lmeControl(msMaxIter = 500, msVerbose = FALSE))
  
  summary(fit5.lme)$tTable[c('zhedonia','zeudaimonia'),]
  
  

    # gene is a random effect
  form <- formula('expression ~ male + age + white + bmi + alcohol + smoke + 
    illness + cd3d + cd3e + cd4 + cd8a + fcgr3a + cd19 + ncam1 + 
                 cd14 + zhedonia + zeudaimonia + (gene|subject)')
  form <- formula('expression ~ gene + male + age + white + bmi + alcohol + smoke + 
    illness + cd3d + cd3e + cd4 + cd8a + fcgr3a + cd19 + ncam1 + 
                  cd14 + zhedonia + zeudaimonia + (1|subject)')
  fit.lmer <- lmer(form,data=dtlong)
  form <- formula('expression ~ gene + male + age + white + bmi + alcohol + smoke + 
    illness + cd3d + cd3e + cd4 + cd8a + fcgr3a + cd19 + ncam1 + 
                  cd14 + zeudaimonia + (1|subject)')
  fit.lmer2 <- lmer(form,data=dtlong)
  anova(fit.lmer,fit.lmer2) # test zhedonia
  form <- formula('expression ~ gene + male + age + white + bmi + alcohol + smoke + 
    illness + cd3d + cd3e + cd4 + cd8a + fcgr3a + cd19 + ncam1 + 
                  cd14 + zhedonia + (1|subject)')
  fit.lmer2 <- lmer(form,data=dtlong)
  anova(fit.lmer,fit.lmer2) # test zeudaimonia
  
  # gee
  R <- cor(dt[,.SD,.SDcols=ycols])
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  fit.gee <- gee(form, data=dtlong,id=gene,corstr='unstructured')
  fit.gee <- gee(form, data=dtlong,id=gene,corstr='fixed',R=R)
  
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  fit.gee <- gee(form, data=dtlong,id=subject,corstr='unstructured')
  fit.gee <- gee(form, data=dtlong,id=subject,corstr='fixed',R=R)
  
  # unrestricted
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  fit2 <- gls(form, data=dtlong, method='ML', correlation=corSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 50, tolerance=1e-3, msVerbose = FALSE))
  summary(fit2)
  
  
  # explore blocks of Y
  R <- cor(Y)
  yvec <- R[lower.tri(R)]
  qplot(yvec,binwidth=0.1)
  y1 <- pro_inflam_genes(2015)
  y2 <- antibody_genes()
  y3 <- ifn_genes()
  R12 <- R[y1,y2]
  R13 <- R[y1,y3]
  R23 <- R[y2,y3]
  mean(R12)
  mean(R13)
  mean(R23)
}

explore_model_simplification <- function(){
  dt <- copy(dt2015)
  dtlong <- get_dtlong(dt,2015)
  # this replicates the analysis of Fredrickson et al.
  # subject is the grouping factor and this is identified in the correlation statement
  xcols2 <- setdiff(xcols,'zeudaimonia')
  form <- formula(paste('expression~',paste(c('gene',xcols2),collapse='+'),sep=''))
  fit2 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  saveRDS(fit2,'FRED15.gls.-zeudaimonia.rds')
  fitlm <- lm(form, data=dtlong) # use this further down
  
  fit1 <- readRDS('FRED15.gls.rds')
 
  fit_null <- gls(expression~1, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  
  summary(fit_null)$AIC
  summary(fit1)$AIC
  summary(fit2)$AIC
  
  # pseudo R^2 mc
  pseudo_R2 <- 1-(as.numeric(logLik(fit1)/logLik(fit_null))) #mcfadden's pesudo R^2
  
  dt <- copy(dt2013)
  dtmeans <- dt[,expression:=scale(apply(scale(dt[,.SD,.SDcols=ycols]),1,mean))]
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  fit13 <- lm(form,data=dtmeans)
  dt <- copy(dt2015)
  dtmeans <- dt[,expression:=scale(apply(scale(dt[,.SD,.SDcols=ycols]),1,mean))]
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  fit15 <- lm(form,data=dtmeans)
  dt <- copy(dtCombi)
  dtmeans <- dt[,expression:=scale(apply(scale(dt[,.SD,.SDcols=ycols]),1,mean))]
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  fitCombi <- lm(form,data=dtmeans)
  summary(fit13)$coefficients[zcols,]
  summary(fit15)$coefficients[zcols,]
  summary(fitCombi)$coefficients[zcols,]
  
  # fitlm - this is a simple lm of expression on xcols not including zeudaimonia. Now plot residuals against zeudaimonia
  dtlong[,resid.zhed:=fitlm$residuals]
  dtlong[,resid.zhed.std:=scale(resid.zhed)]
  dtlong[,fitted.zhed:=fitlm$fitted]
  qplot(x=zeudaimonia, y=resid.zhed,data=dtlong)
  qplot(x=fitted.zhed, y=resid.zhed.std,data=dtlong)
  summary(lm(resid.zhed~zeudaimonia,data=dtlong))
  means <- dtlong[,.(resid.zhed=mean(resid.zhed), zeudaimonia=mean(zeudaimonia)),by=subject]
  qplot(x=zeudaimonia, y=resid.zhed,data=means)
  summary(lm(resid.zhed~zeudaimonia,data=means))
}


bootstrap_models <- function(dt, which_file, tests=c('mv','gls','gee'), niter=200){
  # bootstrap estimates using multivariate regression, glm, and gee models
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  Y <- as.matrix(dt[, .SD, .SDcols=ycols])
  
  # remove smoke
  xcols <- setdiff(xcols, 'smoke')
  
  rows <- 1:nrow(dt)
  mv_matrix <- matrix(0,nrow=niter,ncol=2)
  colnames(mv_matrix) <- zcols
  gee_table <- data.table(NULL)
  gls_table <- data.table(NULL)
  code <- paste(sample(LETTERS,4,replace=TRUE),collapse='')
  gee.out <- paste(which_file,code,'bootstrap.gee.table','txt',sep='.')
  gls.out <- paste(which_file,code,'bootstrap.gls.table','txt',sep='.')
  samp <- 'obs'
  for(iter in 1:niter){
    if('mv' %in% tests){
      Y.samp <- scale(Y[rows,])
      form <- formula(paste('Y.samp~',paste(xcols,collapse='+'),sep=''))
      fit <- lm(form,data=dt[rows,])
      mv_matrix[iter,] <- apply(coefficients(fit)[zcols,], 1, mean)
    }
    if('gls' %in% tests){
      dtlong <- get_dtlong(dt[rows,],year=year,center_Y=TRUE) # returns contrast coefficients
      dtlong <- orderBy(~subject + gene, dtlong)
      form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
      fit.gls <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(msMaxIter = 500, msVerbose = FALSE))
      estimate <- summary(fit.gls)$tTable[zcols,'Value']
      names(estimate) <- zcols
      gls_table <- rbind(gls_table,data.table(samp=samp,t(estimate)))
      write.table(gls_table,gls.out,quote=FALSE,row.names=FALSE,sep='\t')
    }
    if('gee' %in% tests){
      dtlong <- get_dtlong(dt[rows,],year=year,center_Y=TRUE) # returns contrast coefficients
      dtlong <- orderBy(~subject + gene, dtlong)
      form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
      fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='exchangeable', std.err="san.se")
      # summary(lm(form,data=dtlong))$coefficients[zcols,]
      # summary(fit.geeglm)$coefficients[zcols,]
      estimate <- summary(fit.geeglm)$coefficients[zcols,'Estimate']
      names(estimate) <- zcols
      gee_table <- rbind(gee_table,data.table(samp=samp,t(estimate)))
      write.table(gee_table,gee.out,quote=FALSE,row.names=FALSE,sep='\t')
    }
    rows <- sample(1:nrow(dt),replace=TRUE)
    samp <- 'resample'
  }
  
  #apply(gee_table[,.SD,.SDcols=zcols],2,quantile, probs=c(0.025,0.975))
  #ci <- apply(mv_matrix,2,quantile, probs=c(0.025,0.975))
  #apply(mv_matrix,2,quantile, probs=c(0.1,0.9))
  return(NULL)
}

fake_data_gls <- function(dt,xcols,ycols,niter=1000, write_it=FALSE,fn){
  # hybrid of bootstrap and gls_with_correlated_error
  # basically use the bootstrap data to get a distribution of b_hedonia and b_zeudaimonia
  
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  
  p <- length(ycols)
  n <- nrow(dt)
  zcols <- c('zhedonia','zeudaimonia')

  Y <- scale(contrast_coefficients(dt.boot[,.SD,.SDcols=ycols]))
  dt <- cbind(dt[,.SD,.SDcols=xcols],Y)
  dt[,subject:=factor(.I)]
  
  r <- cor(dt[,zhedonia],dt[,zeudaimonia])
  b <- sqrt(r)
  res <- matrix(0,nrow=niter,ncol=2)
  colnames(res) <- c('zhedonia','zeudaimonia')
  for(iter in 1:niter){

    # replace zhedonia and zeudaimonia with fake data
    z <- rnorm(n)
    fake_hed <- scale(b*z + sqrt(1-b^2)*rnorm(n))
    fake_eud <- scale(b*z + sqrt(1-b^2)*rnorm(n))
    # cor(fake_hed, fake_eud)
    # wide to long
    dt[,zhedonia:=fake_hed]
    dt[,zeudaimonia:=fake_eud]
    
    dtlong <- melt(dt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
    dtlong[,gene:=factor(gene)]
    
    # this replicates the analysis of Fredrickson et al.
    form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
    fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(msMaxIter = 500, msVerbose = FALSE))
    res[iter,] <- coefficients(summary(fit1))[c('zhedonia','zeudaimonia')]
    write.table(res,'fake.gls.txt',quote=FALSE,row.names = FALSE)
  }
  qplot(x=res[,'zhedonia'],y=res[,'zeudaimonia'])
  return(NULL)
}

get_Y_matrix <- function(dt, scaled=TRUE, centered=TRUE){
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  Y <- scale(dt[,.SD,.SDcols=ycols], center=centered, scale=scaled)
  
  return(Y)
}

get_X_matrix <- function(dt){
  # dt <- copy(dt2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  dt[,male:=as.integer(as.character(male))]
  dt[,white:=as.integer(as.character(white))]
  dt[,alcohol:=as.integer(as.character(alcohol))]
  dt[,smoke:=as.integer(as.character(smoke))]
  X <- as.matrix(dt[,.SD,.SDcols=xcols])
  return(X)
}

get_dtlong <- function(dt,year=2015,center_Y=TRUE,factor2numeric=FALSE){
  # dt is a raw data set
  # returns scaled contrast coefficients
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  # create matrix with scaled y
  #y.scale <- scale(dt[,.SD,.SDcols=ycols])
  Y <- scale(dt[,.SD,.SDcols=ycols],center=center_Y)
  dts <- cbind(dt[,.SD,.SDcols=xcols], Y)
  dts[,subject:=factor(.I)]
  
  # wide to long
  dtlong <- melt(dts,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  if(factor2numeric==TRUE){
    dtlong[,male:=as.integer(as.character(male))]
    dtlong[,white:=as.integer(as.character(white))]
    dtlong[,alcohol:=as.integer(as.character(alcohol))]
    dtlong[,smoke:=as.integer(as.character(smoke))]
  }
  return(dtlong)
}

simulate_type1_error <- function(run_simulations=FALSE){
  # simulate either power or type I error. 
  if(run_simulations==TRUE){
    # use real data to come up with fake data
    dt <- copy(dt2015)
    if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
    ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
    xcols <- get_xcols()
    dt[,male:=as.integer(as.character(male))]
    dt[,white:=as.integer(as.character(white))]
    dt[,alcohol:=as.integer(as.character(alcohol))]
    dt[,smoke:=as.integer(as.character(smoke))]
    R <- cor(dt[,.SD,.SDcols=ycols])
    Rx <- cor(dt[,.SD,.SDcols=xcols])
    
    n_params <- length(ycols) + length(xcols)
    
    p <- 30
    multiple <- 2
    n <- multiple * (p + length(xcols))
    
    
    bigniter <- 200
    p_table <- matrix(0,nrow=bigniter,ncol=8)
    colnames(p_table) <- c('p1.boot','p2.boot','p1.mv','p2.mv','p1.gls','p2.gls','p1.lme','p2.lme')
    code <- sample(LETTERS,4,replace=TRUE)
    for(bigiter in 1:bigniter){
      inc <- sample(1:nrow(R),p)
      Y <- rmvnorm(n,sigma=R[inc,inc])
      colnames(Y) <- paste("Y",1:p,sep='')
      X <- rmvnorm(n,sigma=Rx)
      colnames(X) <- xcols
      
      
      # bootstrap t
      do_boot <- FALSE
      p1.boot <- 0
      p2.boot <- 0
      if(do_boot==TRUE){
        niter <- 2000
        b1 <- matrix(0,nrow=niter,ncol=p)
        b2 <- matrix(0,nrow=niter,ncol=p)
        rows <- 1:n
        for(iter in 1:niter){
          fit <- lm(Y[rows,]~X1[rows] + X2[rows])
          ss <- summary(fit)
          for(j in 1:p){ # save t-values instead of coefficients
            b1[iter,j] <- ss[[paste('Response ',colnames(Y)[j],sep='')]]$coefficients['X1[rows]', "Estimate"]
            b2[iter,j] <- ss[[paste('Response ',colnames(Y)[j],sep='')]]$coefficients['X2[rows]', "Estimate"]
            rows <- sample(1:n,replace=TRUE)
          }
        }
        b1means <- apply(b1,1,mean)
        b2means <- apply(b2,1,mean)
        b1bar <- mean(b1means)
        b2bar <- mean(b2means)
        b1se <- sd(b1means)
        b2se <- sd(b2means)
        t1 <- abs(b1bar)/b1se
        t2 <- abs(b2bar)/b2se
        p1.boot <- 2*pt(t1,df=(p-1),lower.tail = FALSE)
        p2.boot <- 2*pt(t2,df=(p-1),lower.tail = FALSE)
      }
      
      # gls and lme
      fd <- data.table(X,Y)
      fd[,subject:=factor(.I)]
      # wide to long
      dtlong <- melt(fd,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
      dtlong[,gene:=factor(gene)]
      form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
      fit.gls <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
      #    fit.lme <- lme(form, random=~1|subject, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(maxIter=100, msMaxIter = 500, msVerbose = FALSE))
      
      p_table[bigiter,'p1.boot'] <- p1.boot
      p_table[bigiter,'p2.boot'] <- p2.boot
      form <- formula(paste('Y~',paste(xcols,collapse='+'),sep=''))
      p_table[bigiter,'p1.mv'] <- anova(lm(form, data=fd))['zhedonia','Pr(>F)']
      p_table[bigiter,'p2.mv'] <- anova(lm(form, data=fd))['zeudaimonia','Pr(>F)']
      p_table[bigiter,'p1.gls'] <- summary(fit.gls)$tTable['zhedonia','p-value']
      p_table[bigiter,'p2.gls'] <- summary(fit.gls)$tTable['zeudaimonia','p-value']
      p_table[bigiter,'p1.lme'] <- summary(fit.lme)$tTable['zhedonia','p-value']
      p_table[bigiter,'p2.lme'] <- summary(fit.lme)$tTable['zeudaimonia','p-value']
      write.table(data.table(p=p,p_table),file=paste(paste(code,collapse=''),'power.gls.txt',sep=''),quote=FALSE,sep='\t',row.names=FALSE)
    }
    apply(p_table,2,type.1.error)
  }
  fn <- 'power_list.txt'
  p_table <- data.table(read.table(fn,header=TRUE))
  error_table <- data.table(NULL)
  for(alpha in c(0.1,0.05,0.01,0.001,0.0001)){
    error_table <- rbind(error_table,cbind(alpha=alpha,p_table[,.(p1=type.1.error(p1.gls,alpha=alpha),p2=type.1.error(p2.gls,alpha=alpha)),by=p]))
  }
  error_table[,Error:=(p1+p2)/2]
  error_table[,p:=factor(p)]
  gg <- ggplot(data=error_table,aes(x=alpha,y=Error,color=p))
  gg <- gg + geom_point()
  gg <- gg + geom_line()
  gg <- gg + scale_x_log10()
  gg
  
  error_table_2 <- data.table(NULL)
  for(alpha in c(0.1,0.05,0.01,0.001,0.0001)){
    error_table_2 <- rbind(error_table_2,cbind(alpha=alpha,p_table[,.(p1=type.1.error(p1.gls,alpha=alpha),p2=type.1.error(p2.gls,alpha=alpha))]))
  }
  error_table_2[,Error:=(p1+p2)/2]
  error_table_2[,Inflation:=Error/alpha]
  error_table_2[,alpha:=c('.1','.05','.01','.001','.0001')]
  error_table_2 <- error_table_2[,.(alpha,Error=round(Error,2),Inflation=round(Inflation,1))]
  
  return(error_table_2)
}

simulate_gls_vs_lme <- function(){
  n <- 200
  p <- 50 # number of genes
  gene <- as.factor(rep(1:p,n))
  subject <- as.factor(rep(1:n,each=p))
  niter <- 100
  p_matrix <- matrix(0,nrow=niter,ncol=2)
  colnames(p_matrix) <- c('gls','lme')
  for(iter in 1:niter){
    x <- rep(rnorm(n),each=p) # p values for subject 1, then subject 2
    y <- rnorm(n*p)
    sdt <- data.table(subject=subject,gene=gene,x=x,y=y)
    fit1.gls <- gls(y~gene + x, correlation=corCompSymm(form = ~ 1 | subject))
    p_matrix[iter,'gls'] <- summary(fit1.gls)$tTable['x','p-value']
    fit1.lme <- lme(y~gene + x, random = ~1|subject, correlation=corCompSymm(form = ~ 1 | subject))
    p_matrix[iter,'lme'] <- summary(fit1.lme)$tTable['x','p-value']
  }
 apply(p_matrix,2,type.1.error) # nothing here to suggest inflated type I error
}

explore_glh_and_gls <- function(){
  dt <- copy(dt2015)
  inc <- sample(1:length(ycols),5)
  Y <- scale(as.matrix((dt[,.SD,.SDcols=ycols[inc]])))
  xcols2 <- c('bmi','illness','cd3d','zhedonia','zeudaimonia')
  dts <- cbind(dt[,.SD,.SDcols=xcols2], Y)
  dts[,subject:=factor(.I)]
  
  # wide to long
  dtlong <- melt(dts,id.vars=c('subject',xcols2),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  
  # this replicates the analysis of Fredrickson et al.
  # subject is the grouping factor and this is identified in the correlation statement
  form <- formula(paste('expression~',paste(c('gene',xcols2),collapse='+'),sep=''))
  fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit1)
  
  form <- formula(paste('expression~',paste(c('gene',xcols2),collapse='+'),sep=''))
  fit2 <- gls(form, data=dtlong, method='ML', correlation=corSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 50, tolerance=1e-6, msVerbose = FALSE))
  summary(fit2)
  
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit3 <- lm(form, data=dt)
  linearHypothesis(fit3, "zeudaimonia = 0")
  Anova(fit3,type='2') # confirm that these are testing same
  
  # pls
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit4 <- plsr(form, data=dt)
  
}

explore_mixed_model <- function(){
  fn <- 'cole1_clean.txt'
  dt <- read_file(fn,year=2013)
  dt[,subject:=factor(.I)]
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year=2013),antibody_genes(),ifn_genes())
  dtlong <- melt(dt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  
  
  
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep='')) # doesn't make sense to model gene as a fixed effect but instead as another random intercept
  form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  # lme(CALLS ~ YEAR, random =~1 | BLOCK, data = driscoll, correlation = corSymm(form = ~1 | BLOCK))
  
 
  model1 <- lme(
    expression ~ zhedonia,
    random = ~ 1|gene,
    #random = list(~ 1|subject, ~ 1|gene),
    #random = list(subject=~1, gene=~1),
    #random = list(gene=~1, subject=~1),
    correlation = corSymm(form = ~1 | gene),
    data = dtlong, method = "ML"
  )
  summary(model1)
  
  model1 <- lme(
    form,
    #random = ~ 1|subject,
    random = list(~ 1|subject, ~1|gene),
    data = dtlong, method = "ML"
  )
  summary(model1)
  coefficients(summary(model1))
  
  model1 <- lmer(expression~zhedonia + (1|subject) + (1|gene), data=dtlong,REML=FALSE)
  summary(model1)
  model2 <- lmer(expression~1 + (1|subject) + (1|gene), data=dtlong,REML=FALSE)
  anova(model1,model2, test='Chisq')
  
  model1 <- lmer(expression~zhedonia*gene + (1|subject), data=dtlong,REML=FALSE)
  summary(model1)
  
  
}

explore_multivariate_regression <- function(){
  fn <- 'cole1_clean.txt'
  dt <- read_file(fn,year=2013)
  
  xcols <- c('bmi','cd8a','zhedonia')
  Xobs <- dt[,.SD,.SDcols=xcols]
  Yobs <- dt[,.SD,.SDcols=ycols[c(1,2,7)]]

  n <- 1000
  k <- 3 # number of X
  p <- 4 # number of Y
  X <- rmvnorm(n, sigma=cor(Xobs))
  B <- matrix(0,nrow=k,ncol=p)
  B[1,] <- c(-.5,0,0, .5)
  E <- matrix(rnorm(n*p),nrow=n)
  Y <- X%*%B + E
  fit <- lm(Y~X)
  summary(fit)
  linearHypothesis(fit, "X1 = 0",verbose = TRUE)
  
  # make x correlated to the minor axes of cor(Y1,Y2)
  n <- 1000
  r <- .8
  R <- matrix(c(1,r,r,1),nrow=2)
  Y <- rmvnorm(n,sigma=R)
  Y1 <- Y[,1]
  Y2 <- Y[,2]
  E <- eigen(cov(Y))$vectors
  pc2 <- Y%*%E[,2]
  x <- Y%*%E[,2] + sqrt(1-.5)*rnorm(n)
  cor(x,pc2)
  cor(x,Y1)
  cor(x,Y2)
  qplot(x,pc2)
  qplot(x=Y1,y=x)
  qplot(x=Y2,y=x)
  mvfit <- lm(Y~x)
  summary(mvfit)
  anova(lm(Y1~x))
  anova(lm(Y2~x))
  anova(lm(pc2~x))
  linearHypothesis(mvfit, "x = 0",verbose = TRUE)
  
  # right. The GLH is testing if x "explains" some axis through multidimensional Y space. This doesn't directly address if the mean coefficient is different from zero but does address if there is some multivariate effect. If there is no multivariate effect, then the mean coefficient is not different from zero.
  
  
  # make 2 Xs, each correlated to the minor axes of cor(Y1,Y2) so correlated with each other (like hedonia and eudaimonia). Test for difference between the two.
  n <- 100
  r <- .8 # correlation between the x
  b <- sqrt(r)
  R12 <- 0.2 # correlation in Y
  R <- matrix(c(1,R12,R12,1),nrow=2)
  Y <- rmvnorm(n,sigma=R)
  Y1 <- Y[,1]
  Y2 <- Y[,2]
  E <- eigen(cov(Y))$vectors
  pc2 <- Y%*%E[,2]
  x1 <- (b+.2)*pc2 + sqrt(1-r)*rnorm(n)
  x2 <- (b-.2)*pc2 + sqrt(1-r)*rnorm(n)
  x1 <- (.1)*pc2 + sqrt(1-r)*rnorm(n)
  x2 <- (-.1)*pc2 + sqrt(1-r)*rnorm(n)
  mvfit <- lm(Y~x1 + x2)
  coefficients(mvfit)
  linearHypothesis(mvfit, "x1 = 0",verbose = TRUE)
  linearHypothesis(mvfit, "x2 = 0",verbose = TRUE)
  linearHypothesis(mvfit, "x1 = x2",verbose = TRUE)
  
  # one can have two vectors each of which is not sig diff from zero but are sig diff from each other. 
}

explore_canonical_correlation <- function(){
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  
  Xobs <- scale(get_X_matrix(dt))
  Yobs <- scale(as.matrix(dt[,.SD,.SDcols=ycols]))
  
  # look at PCA of Yobs
  PCA <- eigen(cov(Yobs)) # this is of cor since Yobs is scaled
  eigenvectors <- PCA$vectors
  eigenvalues <- PCA$values
  # note that there isn't a vector where most/all of the coefficients have same sign so where they are co-varying together. This means that something like latent variable regression wouldn't achieve the goal of the 
  
  # canonical correlations on the residuals on xcols will capture the space through Yresid most correlated with zhedonia and zeudaimonia
  
  # get residuals on xcols2
  xcols2 <- setdiff(xcols,zcols)
  form <- formula(paste('Yobs',paste(xcols2,collapse='+'),sep='~'))
  Yresid <- lm(form,data=dt)$residuals
  fit <- cc(Xobs[,zcols], Yresid)
  fit$xcoef
  loads <- comput(Xobs[,zcols],Yresid,fit)
  cor(fit$scores$xscores,fit$scores$yscores)

}

explore_gee <- function(){
  # gee
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  dtlong <- get_dtlong(dt, year=2015, center=TRUE)
  dtlong <- orderBy(~subject + gene,dtlong)
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  #form <- formula(paste('expression~',paste(xcols,collapse='+'),sep=''))
  
  fit.lm <- lm(form,data=dtlong)
  summary(fit.lm)$coefficients[zcols,]
  
  fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject, waves=gene, corstr='independence', std.err="san.se")
  summary(fit.geeglm)$coefficients[zcols,]
  
  fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject, waves=gene, corstr='exchangeable', std.err="san.se")
  summary(fit.geeglm)$coefficients[zcols,]
  qplot(fit.geeglm$fitted.values,fit.geeglm$residuals)
  # compare this to the GLS fit
  
  
  fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject, waves=gene, corstr='unstructured', std.err="san.se")
  summary(fit.geeglm)$coefficients[zcols,]
  
  
  Y <- get_Y_matrix(dt)
  form2 <- formula(paste('Y~',paste(xcols,collapse='+'),sep=''))
  fit.mv <- lm(form2, data=dt)
  Rresid <- cor(fit.mv$residuals)
  zcor <- fixed2Zcor(Rresid, id=dtlong$subject, waves=dtlong$gene)
  fit.geeglm.fix <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='userdefined', zcor=zcor, std.err="san.se")
  summary(fit.geeglm.fix)$coefficients[zcols,]
  
  fit.gls <- gls(form, data=dtlong, method='ML', correlation=corCompSymm( form = ~ 1 | subject), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit.gls)$tTable[zcols,]

  fit.gls.hcs <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit.gls.hcs)$tTable[zcols,]
  S.hcs <- getVarCov(fit.gls.hcs)
  R.hcs <- cor(S.hcs)
  R.gls <- summary(fit.gls.hcs)$corBeta[1:52,1:52]
  # This doesn't work because the intercept. Make first row and col ~ .46
  R.gls[1,] <- .46
  R.gls[,1] <- .46
  R.gls[1,1] <- 1.0
  zcor <- fixed2Zcor(R.hcs, id=dtlong$subject, waves=dtlong$gene)
  fit.geeglm.fix <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='userdefined', zcor=zcor, std.err="san.se")
  summary(fit.geeglm.fix)$coefficients[zcols,]

  fit.lme <- lme(form, random = ~1|subject, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(maxIter=100, msMaxIter = 500, tolerance=1e-6, msVerbose = FALSE)) # tolerance=1e-6 default
  S <- getVarCov(fit.lme, type='conditional')[[1]]
  
}

explore_MRCE <- function(){
  # multivariate regression with correlated response
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  Xobs <- get_X_matrix(dt)
  Yobs <- as.matrix(dt[,.SD,.SDcols=ycols])
  fit.mrce <- mrce(X=Xobs,Y=Yobs,method='single',lam1=10^(-1.5), lam2=10^(-0.5))
  data.table(predictor=colnames(Xobs),bean_beta=apply(fit.mrce$Bhat[,],1,mean))
  
  fit.mrce$mx
  coeffs <- t(fit.mrce$Bhat[16:17,])
  apply(coeffs,2,mean)
}

explore_prediction <- function(){
  # how good is prediction of 2013 data from 2015 model? To do this we need the gls result for 2013 without IL6 in the model
  # need to get dtlong2013 manually withou IL6
  dtlong_2013 <- get_dtlong(dt2013,year=2015,factor2numeric=TRUE)
  #Yhat_2013 <- predict(fit,newdata=dtlong_2013) # why doesn't this work?
  dtlong_2015 <- get_dtlong(dt2015,year=2015,factor2numeric=TRUE)
  #Yhat_2015 <- predict(fit,dtlong_2015)
  X15 <- as.matrix(dtlong_2015[,.SD,.SDcols=get_xcols()])
  coefs <- coefficients(fit)
  inc <- (length(coefs)-16):length(coefs)
  xcoefs <- coefs[inc]
  yhat_15 <- X15%*%xcoefs + xcoefs[1]
  yhat_15_fit <- fit$fitted
  # data.table(yhat_15,yhat_15_fit)[1:10,] # check!
  # now fit FRED13 to FRED15 model
  X13 <- as.matrix(dtlong_2013[,.SD,.SDcols=get_xcols()])
  yhat_13 <- X13%*%xcoefs + coefs[1]
  # now fit FRED13 to FRED13 model
  fit13 <- readRDS(paste(which_file='FRED13.yr2015.gls.rds',sep=''))
  yhat_13_fit <- fit13$fitted
  yresid_13_fit <- fit13$residuals
  y_13 <- dtlong_2013[,expression]
  qplot(y_13,yhat_13_fit,color=dtlong_2013$gene)
  qplot(y_13,yresid_13_fit,color=dtlong_2013$subject)
  qplot(yhat_13_fit,scale(yresid_13_fit),color=dtlong_2013$subject)
  plot(fit13)
  
  sum((yresid_13_fit)^2)/sum((y_13 - mean(y_13))^2)
  cor(yhat_13_fit,dtlong_2013[,expression])
  qplot(yhat_13_fit13,yhat_13)
  cor(yhat_13_fit,yhat_13)
  cor(yhat_13_fit,dtlong_2013[,expression])
  cor(yhat_13,dtlong_2013[,expression])
  
}

explore_globalANCOVA <- function(){
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  covs <- setdiff(xcols,zcols)
  
  Y <- t(scale(dt[, .SD, .SDcols=ycols])) # GlobalAncova uses a m genes X n subjects matrix
  form.full <- formula(paste('~',paste(c(covs,zcols),collapse='+'),sep=''))
  model.dat <- dt[, .SD, .SDcols=xcols]
  GA_hed <- GlobalAncova(Y, form.full, model.dat=model.dat, test.terms='zhedonia', method='permutation')$test.result['p.perm',1]
  GA_eud <- GlobalAncova(Y, form.full, model.dat=model.dat, test.terms='zeudaimonia', method='permutation')$test.result['p.perm',1]
  
}

explore_bull <- function(){
  # gee
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  zcols <- c('zhedonia','zeudaimonia')
  dtlong <- get_dtlong(dt, year=2015, center=TRUE)
  dtlong <- orderBy(~subject + gene,dtlong)
  n <- nrow(dt)
  k <- length(ycols)
  nk <- n*k
  
  Y <- get_Y_matrix(dt)
  ols_form <- formula(paste('Y~',paste(xcols,collapse='+'),sep=''))
  fit.mv <- lm(ols_form, data=dt)
  B <- t(coefficients(fit.mv)[zcols,])
  Rresid <- cor(fit.mv$residuals) #working correlation matrix
  zcor <- fixed2Zcor(Rresid, id=dtlong$subject, waves=dtlong$gene)
  gee_form <- formula(paste('expression~',paste(c('gene',xcols,'-1'),collapse='+'),sep=''))
  fit <- geeglm(gee_form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='userdefined', zcor=zcor, std.err="san.se")
  fit1 <- geeglm(gee_form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='exchangeable', std.err="san.se")
  summary(fit)$coefficients[inc,]
  summary(fit1)$coefficients[zcols,]
  
   # the robust covariance matrix of beta is:
  inc <- 1:52
  V <- fit$geese$vbeta[inc,inc]
  V1 <- fit1$geese$vbeta[inc,inc]
  R <- cov2cor(V) # the correlation matrix of the coefficients (beta) - can this be estimated using bootstrap?
  R1 <- cov2cor(V1) # the correlation matrix of the coefficients (beta) - can this be estimated using bootstrap?
  mean(abs(R[lower.tri(R)]))
  mean(abs(R1[lower.tri(R1)]))
  
  # Bull
  rse <- sqrt(diag(V))
  C <- rse%*%solve(V)
  w <- C[1,]
  Wald <- C%*%B
  pchisq(abs(Wald), 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  
  # obrien OLS/GLS
  Ri <- cov2cor(V) # the correlation matrix of the coefficients (beta) - can this be estimated using bootstrap?
  J <- matrix(1,nrow=nrow(Ri),ncol=1)
  num <- t(J)%*%B
  denom <- c(sqrt(t(J)%*%Ri%*%J))
  Tols <- num/denom
  df <- effective_size(nrow(dt),Rresid)
  2*pt(abs(Tols),df=df,lower.tail = FALSE) #df certaintly less than nk
  
  #using bootstrap
  which_file <- 'FRED15'
  fn <- paste(which_file,'.bootstrap.txt',sep='')
  boot_coeffs <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  Ri <- cor(boot_coeffs[,.SD, .SDcols=ycols])
  mean(abs(Ri[lower.tri(Ri)]))
  J <- matrix(1,nrow=nrow(Ri),ncol=1)
  num <- t(J)%*%B
  denom <- c(sqrt(t(J)%*%Ri%*%J))
  Tols <- num/denom
  df <- effective_size(nrow(dt),Rresid)
  2*pt(abs(Tols),df=df,lower.tail = FALSE) #df certaintly less than nk
  
  V <- Ri
  rse <- sqrt(diag(V))
  C <- rse%*%solve(V)
  w <- C[1,] # notice negative weights
  Wald <- C%*%B
  pchisq(abs(Wald), 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  
  # Roast
  Yt <- t(Y)
  X <- get_X_matrix(dt)
  hedonia <- which(colnames(X)=='zhedonia')
  eudaimonia <- which(colnames(X)=='zeudaimonia')
  roast(y=Yt,design=X,contrast=hedonia)
  roast(y=Yt,design=X,contrast=eudaimonia)
  
  num <- t(J)%*%Ri%*%B
  denom <- c(sqrt(t(J)%*%Ri%*%J))
  Tgls <- num/denom # affected by negative weights

  S <- matrix(c(1,0,0,1),nrow=2)/10
  Si <- solve(S)
  rse <- sqrt(diag(S))
  C1 <- rse%*%Si
  
  S <- matrix(c(1,.5,.5,1),nrow=2)/10
  Si <- solve(S)
  rse <- sqrt(diag(S))
  C2 <- rse%*%Si
 
  S <- matrix(c(1,.9,.9,1),nrow=2)/10
  Si <- solve(S)
  rse <- sqrt(diag(S))
  C3 <- rse%*%Si
  
  S <- matrix(c(1,-.5,-.5,1),nrow=2)/10
  Si <- solve(S)
  rse <- sqrt(diag(S))
  C4 <- rse%*%Si
  
  
  
  B <- c(.6, -.4)
  W <- c(C1%*%B, C2%*%B, C3%*%B, C4%*%B)
  pchisq(W, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

  W <- B^2/diag(S)
  pchisq(W, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  
}

explore_bulls_global_hypothesis <- function(){
  n <- 200
  p <- 50 # the number of outcomes
  # set up two covariates that have mean response zero on the p outcomes
  b1 <- sample(c(rep(0.2,p/2),rep(-.2,p/2)))
  b2 <- sample(c(rep(0.4,p/2),rep(-.4,p/2)))
  B <- as.matrix(cbind(b1=b1,b2=b2))
  X <- rmvnorm(n, sigma=diag(c(1,1)))
  E <- rmvnorm(n, sigma=diag(rep((1-.2^2 -.4^2),p)))
  Y <- t(X%*%t(B) + E)
  df <- data.frame(X)
  colnames(df) <- c('x1','x2')
  
  form.full <- formula('~x1 + x2')
  fit <- GlobalAncova(Y, form.full, model.dat=df, test.terms='x1', method='permutation')
  
  # Roast
  x1 <- which(colnames(df)=='x1')
  x2 <- which(colnames(df)=='x2')
  roast(y=Y,design=df,contrast=x1)
  roast(y=Y,design=df,contrast=x2)
  
}

simulate_the_correlated_coefficients <- function(){
  # the correlation between X1 and X2 is .79. The off-diagonal element of the XtX-1 is big and negative. This means that the coefficients of Y on X1 and X2 will be negatively correlated if X1 and X2 have a positive correlation with Y.
  
  n <- 122
  r <- 0.79
  b <- sqrt(r)
  niter <- 1000
  fdres <- data.frame(matrix(0,nrow=niter,ncol=4))
  for(iter in 1:niter){
    z <- rnorm(n)
    x1 <- b*z + sqrt(1-r)*rnorm(n)
    x2 <- b*z + sqrt(1-r)*rnorm(n)
    y <- rnorm(n)
    fd <- data.table(y=y,x1=x1,x2=x2)
    fit <- lm(y~x1 + x2, data=fd)
    fdres[iter,] <- c(coefficients(fit)[c('x1','x2')],vif(fit))
  }
  colnames(fdres) <- c('x1','x2')
  qplot(x=x1,y=x2,data=fdres)
}

cite_package <- function(){
  toBibtex(citation('geepack'))
  toBibtex(citation('mvtnorm'))
  toBibtex(citation('limma'))
  
  
}
