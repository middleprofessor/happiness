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
library(MRCE)
# plot libraries
library(ggplot2)
library(showtext) # needed for eps fonts to add Arial to .eps
font.add('Arial',regular='Arial.ttf')
library(gridExtra)


# lmperm downloaded from https://github.com/kabacoff/RiA2/tree/master/lmPerm
# using R-Studio chose Tools > Install Packages and the Install From > Package Archive pop-up menu.

# originally coefficients computed raw then reversed during analysis but to compute the permutation GLH, these need to be reversed, so...
# January 13 recomputed by reversing expression levels and not coefficients

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

do_Fredrickson <- function(){
  
  # generate resampled files to compute results
  are_tests_run <- TRUE
  fn <- 'cole1_clean.txt'
  dt2013 <- read_file(fn,year=2013)

  fn <- 'cole2_clean.txt'
  dt2015 <- read_file(fn,year=2015)
  
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

  if(are_tests_run==FALSE){
    do_tests(dt2013,which_file='FRED13',niter=niter,exclude_smoke = FALSE)
    do_tests(dt2015,which_file='FRED15',niter=niter,exclude_smoke = TRUE)
    do_tests(dtCombi,which_file='FRED.Combi',niter=niter,exclude_smoke = FALSE)
  }
  if(are_gls_permutations_run==FALSE){
    
  }
  if(are_gls_parametric_tests_run==FALSE){
    fit1 <- 1
    saveRDS(gls_with_correlated_error(dt2013), "FRED13.gls.rds")
    saveRDS(gls_with_correlated_error(dt2015), "FRED15.gls.rds")
    saveRDS(gls_with_correlated_error(dtCombi), "FRED.Combi.gls.rds")
  }
  
  # print results. Have to do this semi-manually because I cannot get the glh results into the table except by looking at them.
  which_file_list <- c('FRED13', 'FRED15', 'FRED.Combi')
  naive_p_table <- data.table(NULL)
  gls_table <- data.table(NULL)
  gls_matrix <- matrix(0,nrow=1,ncol=4)
  p_matrix <- matrix(0,nrow=15,ncol=3)
  colnames(p_matrix) <- c('zhedonia','zeudaimonia','delta')
  effect_matrix <- matrix(0,nrow=3,ncol=2) # matrix of backtransformed effects
  colnames(effect_matrix) <- c('zhedonia','zeudaimonia')
  row.names(effect_matrix) <- which_file_list
  data_set <- rep(which_file_list,each=5)
  test <- rep(c('bootstrap_t','GLH','GLH_perm','t_perm','CE.LM_perm'),3)
  for(which_file in which_file_list){
    if(which_file=='FRED13'){dt <- copy(dt2013)}
    if(which_file=='FRED15'){dt <- copy(dt2015)}
    if(which_file=='FRED.Combi'){dt <- copy(dtCombi)}
    
    # get permutation data, the first row is the observed data
    fn <- paste(which_file,'.permutation.txt',sep='')
    res <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
    res <- contrast_coefficients(res) # convert to contrasts
    
    # first compute stats as in FRED13 and FRED15 for confirmation that my data/analyses are the same
    # compute back-transformed effect size
    effect_matrix[which_file,] <- effect_size(res[I==1,]) # effect size for observed data
    # do naive t test first following FRED13
    naive_t_table <- naive_t_stats(res[I==1,]) #(reported in FRED13: eudaimonic, P = 0.0045; hedonic, P = 0.0047)
    naive_p_table <- rbind(naive_p_table, data.table(data=which_file,naive_t_table[stat=='p.value',]))
    # gls coefficients and p-values
    fit <- readRDS(paste(which_file,'.gls.rds',sep=''))
    gls_matrix[1,] <- round(c(summary(fit)$tTable[c('zhedonia'),c('Value','p-value')],summary(fit)$tTable[c('zeudaimonia'),c('Value','p-value')]),3)
    colnames(gls_matrix) <- c('zhedonia','p-value','zeudaimonia','p-value')
    gls_table <- rbind(gls_table,data.table(data=which_file,gls_matrix))
    
    
    # now get new stats
    # get permutation stats
    p_matrix[data_set==which_file & test=='t_perm',] <- permutation_t(res,statistic='t')
    p_matrix[data_set==which_file & test=='GLH_perm',] <- permutation_lambda(data.table(read.table(paste('wilks.',fn,sep=''),header=TRUE)))
    
    # compute GLH
    glh <- GLH(dt)
    print(glh$GLH1,SSP=FALSE,SSPE=FALSE) # zhedonia
    print(glh$GLH2,SSP=FALSE,SSPE=FALSE) # zeudaimonia
    print(glh$GLH3,SSP=FALSE,SSPE=FALSE) # delta
    # fill in by hand. Should just use the matrices to compute my own p.value but this fails because I can't get the parameter n correct to compute df and F
    glh_FRED13 <- c(zhedonia=0.2793,zeudaimonia=0.48086,delta=0.33484)
    glh_FRED15 <- c(zhedonia=0.68856,zeudaimonia=0.44726,delta=0.70753)
    glh_Combi <- c(zhedonia=0.50954,zeudaimonia=0.1491,delta=0.32279)
    p_matrix[data_set=='FRED13' & test=='GLH',] <- glh_FRED13
    p_matrix[data_set=='FRED15' & test=='GLH',] <- glh_FRED15
    p_matrix[data_set=='FRED.Combi' & test=='GLH',] <- glh_Combi
    toBibtex(citation("ouch"))
    #compute permutation.gls stats
    # do NOT contrast_coefficient these as this was done during the gls
    fn <- paste(which_file, '.permutation.gls.list.txt',sep='')
    #fn <- 'zeudaimoniaFRED15.permutation.gls.list.txt'
    res <- read_permutation.gls_list(fn)
    p_matrix[data_set==which_file & test=='CE.LM_perm',] <- permutation.gls.p.value(res, statistic='t')
    p.adjust(permutation.gls.p.value(res, statistic='t'),method='fdr')
    
    permutation.gls.p.value(res, statistic='t',do_ci=TRUE) # get 95% CI on p-values
    
    # compute bootstrap stats
    fn <- paste(which_file,'.bootstrap.txt',sep='')
    res <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
    res <- contrast_coefficients(res) # convert to contrasts
    p_matrix[data_set==which_file & test=='bootstrap_t',] <- smart_t_stats(res)[,p.value]
    
  }
  p_table <- data.table(data=data_set,test=test,round(p_matrix,2))
  naive_p_table[,zhedonia:=round(zhedonia,3)]
  naive_p_table[,zeudaimonia:=round(zeudaimonia,3)]
  naive_p_table[,delta:=round(delta,3)]
  naive_p_table <- naive_p_table[,.SD,.SDcols=c('data','zhedonia','zeudaimonia','delta')]
  effect_matrix <- data.table(data=row.names(effect_matrix),(round(effect_matrix,1)))
  p_table
  naive_p_table
  effect_matrix
  gls_table
  
  write.table(effect_matrix,'effect_matrix.txt',sep='\t',quote=FALSE,row.names=FALSE)
  write.table(naive_p_table,'naive_p_table.txt',sep='\t',quote=FALSE,row.names=FALSE)
  write.table(p_table,'p_table.txt',sep='\t',quote=FALSE,row.names=FALSE)
  write.table(gls_table,'gls_table.txt',sep='\t',quote=FALSE,row.names=FALSE)
  
}

figure_1 <- function(){
  # A scatterplot of contrast coefficients for x=hedonia y=eudaimonia for two of the randomly permuted runs using the FRED.Combi data
  which_file <- 'FRED13'
  fn <- paste(which_file,'.permutation.txt',sep='')
  part <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  part <- contrast_coefficients(part) # convert to contrasts
  dt <- data.table(data=which_file,part)
  which_file <- 'FRED15'
  fn <- paste(which_file,'.permutation.txt',sep='')
  part <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  part <- contrast_coefficients(part) # convert to contrasts
  dt <- rbind(dt,data.table(data=which_file,part),fill=TRUE)
  which_file <- 'FRED.Combi'
  fn <- paste(which_file,'.permutation.txt',sep='')
  part <- data.table(read.table(fn,header=TRUE))  # read in raw (not multiplied by contrast coeff) regression coefficients
  part <- contrast_coefficients(part) # convert to contrasts
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
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(.85,.9))
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
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position=c(.85,.9))
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

do_tests <- function(dt,which_file,niter=2000, exclude_smoke=FALSE){
  # dt is a data.table of the responses and regressors
  # which_file is the file FRED13 or FRED15 Fredrickson et. al.
  # exclude_smoke=TRUE, see below)
  
  # get xcols and ycols
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  permutation_test_anderson(dt,xcols,ycols,niter=niter,write_it=TRUE,fn=paste(which_file,'.permutation.txt',sep=''))

  #gls permutation
  zcols <- c('zhedonia','zeudaimonia')
  permutation_gls(dt,xcols=xcols,ycols=ycols,zcols=zcols,niter=101,do_obs=TRUE,write_it=TRUE,fn=paste(which_file,'.permutation.gls',sep=''))
  
  #gls permutation on individual measures of happiness
  zcols <- c('zhedonia','zeudaimonia')
  happy <- 'zeudaimonia' # limit analysis to zcols
  red_xcols <- setdiff(xcols,setdiff(zcols,happy)) # remove the other happy from xcols
  permutation_gls(dt,xcols=red_xcols,ycols=ycols,zcols=happy,niter=80,do_obs=FALSE,write_it=TRUE,fn=paste(paste(happy,which_file,sep=''),'.permutation.gls',sep=''))
  
  # note that in FRED13, smoke has only 6/77 scored as 1 and some bootstraps will entirely miss this so
  if(exclude_smoke==TRUE){xcols <- setdiff(xcols,'smoke')}
  bootstrap_test(dt,xcols.boot,ycols,niter=niter,write_it=TRUE, fn=paste(which_file,'.bootstrap.txt',sep=''))
  
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

expression_coefficient <- function(ycols){
  # cc is the contrast coefficient (-1 or 1) set in original paper
  # 19 proinflammatory genes, which are up-regulated on average in the CTRA

  p <- length(ycols)
  if(p==53){year <- 2013}else{year <- 2015}
  cc <- rep(1.0, p)
  pro_inflam <- pro_inflam_genes(year)
  for(j in 1:length(ycols)){
    if(ycols[j] %in% pro_inflam){cc[j] <- 1}else{cc[j] <- -1}
  }
  return(cc)
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

permutation_test_anderson <- function(dt,xcols,ycols,niter=1000, write_it=FALSE,fn){
  # permutation test based on Anderson and Robinson 2001
  # dt is a data.table with the X regressors and Y responses
  # xcols are the regressors
  # ycols are the responses
  # fn is the file name to write to
  # notes: the expected association between permuted hedonic score and gene expression is zero so the expected delta is zero E(E(b.h)-E(b.e))=0-0.
  
  p <- length(ycols)
  Y <- as.matrix(dt[,.SD,.SDcols=ycols])
  
  # coefficients (save t-value instead of coefficients)
  coeffs.hed <- matrix(0,nrow=niter,ncol=p)
  coeffs.eud <- matrix(0,nrow=niter,ncol=p)
  
  # wilks GLH result
  dt.cc <- contrast_coefficients(copy(dt)) # need to do this on the reversed expression levels because I cannot reverse the coefficients later
  wilks_matrix <- matrix(0,nrow=niter,ncol=3) # matrix of Wilk's Lambda
  colnames(wilks_matrix) <- c('zhedonia=0','zeudaimonia=0','delta=0')

    
  zcols <- c('zhedonia','zeudaimonia')
  xcols2 <- setdiff(xcols, zcols)
  # Get residuals from X (so excluding zhedonia and zeudaimonia)
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit.obs <- lm(form, data=dt)
  yhat <- predict(fit.obs) # Yhat = aX
  rows <- 1:nrow(dt) # observed on first iter and permuted after
  for(iter in 1:niter){
    e <- residuals(fit.obs)[rows,] # permuted Yhat - aX
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
      coeffs.hed[iter,j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zhedonia', "t value"]
      coeffs.eud[iter,j] <- ss[[paste('Response ',ycols[j],sep='')]]$coefficients['zeudaimonia', "t value"]
    }

    # do wilks on reversed Y
    form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
    fitmv <- lm(form, data=dt.cc)
    yhat <- predict(fitmv) # Yhat = aX
    e <- residuals(fitmv)[rows,] # permuted Yhat - aX
    Ypi <- yhat + e
    form <- formula(paste('Ypi',paste(xcols,collapse='+'),sep='~'))
    fitmv.pi <- lm(form, data=dt.cc)
    glh <- linearHypothesis(fitmv.pi, "zhedonia = 0")
    wilks.zhed <- det(glh$SSPE)/det(glh$SSPE + glh$SSPH) 
    glh <- linearHypothesis(fitmv.pi, "zeudaimonia = 0")
    wilks.zeud <- det(glh$SSPE)/det(glh$SSPE + glh$SSPH) 
    glh <- linearHypothesis(fitmv.pi, "zeudaimonia = zhedonia")
    wilks.delta <- det(glh$SSPE)/det(glh$SSPE + glh$SSPH) 
    wilks_matrix[iter,] <- c(wilks.zhed,wilks.zeud,wilks.delta)
    
    # permute rows
    rows <- sample(1:nrow(dt))
  }
  
  b_matrix <- data.table(I=rep(1:niter,2), Type=rep(c('hedonic','eudaimonic'),each=niter), rbind(coeffs.hed,coeffs.eud))
  
  if(write_it==TRUE){write.table(b_matrix,fn,quote=FALSE,sep='\t',row.names=FALSE)}
  fn <- paste('wilks.',fn,sep='')
  if(write_it==TRUE){write.table(wilks_matrix,fn,quote=FALSE,sep='\t',row.names=FALSE)}
  
  return(NULL)
}

permutation_gls <- function(dt,xcols,ycols,zcols=c('zhedonia','zeudaimonia'),niter=100,do_obs=TRUE,write_it=FALSE,fn){
  # hybrid of permutation_anderson and gls_with_correlated_error
  # basically use the permuted data to get a null distribution of b_hedonia and b_zeudaimonia
  
  #fn <- 'cole2_clean.txt'
  #dt <- read_file(fn,year=2015)
  #xcols <- get_xcols()
  #ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  
  # create matrix with scaled y
  #y.scale <- scale(dt[,.SD,.SDcols=ycols])
  Y <- scale(contrast_coefficients(dt[,.SD,.SDcols=ycols]))
  dt <- cbind(dt[,.SD,.SDcols=xcols],Y)
  dt[,subject:=factor(.I)]
  
  p <- length(ycols)
  xcols2 <- setdiff(xcols, zcols)
  # Get residuals from X (so excluding zhedonia and zeudaimonia)
  form <- formula(paste('Y',paste(xcols2,collapse='+'),sep='~'))
  fit.obs <- lm(form, data=dt)
  yhat <- predict(fit.obs) # Yhat = aX
  if(do_obs==TRUE){
    rows <- 1:nrow(dt) # observed on first iter and permuted after
  }else{
    rows <- sample(1:nrow(dt))
  }
  
  # coefficients (save t-value in addition to coefficients)
  b_matrix <- matrix(0,nrow=niter,ncol=2*length(zcols))
  #colnames(b_matrix) <- c('coeff.zhedonia','coeff.zeudaimonia','t.zhedonia','t.zeudaimonia')
  colnames(b_matrix) <- c(paste('coeff',zcols,sep='.'),paste('t',zcols,sep='.'))


  code <- sample(LETTERS,4,replace=TRUE)
  fn_full <- paste(fn,'.',paste(code,collapse=''),'.txt',sep='')
  fn_temp <- paste(fn,'.',paste(code,collapse=''),'.temp.txt',sep='')
  
  
  for(iter in 1:niter){
    e <- residuals(fit.obs)[rows,] # permuted Yhat - aX
    Ypi <- scale(yhat + e)
    dt <- cbind(dt[,.SD,.SDcols=xcols],Ypi)
    dt[,subject:=factor(.I)] # need to recalc since dt is re-created
    # check apply(dt[,.SD,.SDcols=ycols],2,sd)
    
    # wide to long
    dtlong <- melt(dt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
    dtlong[,gene:=factor(gene)]
    
    # this replicates the analysis of Fredrickson et al.
    form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
    fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
    b_matrix[iter,] <- c(
      summary(fit1)$tTable[zcols,'Value'],
      summary(fit1)$tTable[zcols,'t-value'])
    # zhedonia     0.02490  0.127834  0.194792  0.8456
    # zeudaimonia -0.19396  0.131233 -1.478016  0.1395

    # permute rows for next iteration
    rows <- sample(1:nrow(dt))
    
    # partial writing so I can check results
    write.table(b_matrix,fn_temp,quote=FALSE,row.names = FALSE)
  }
  
  if(do_obs==TRUE){
    b_matrix <- data.table(permutation=c('obs',rep('perm',niter-1)),b_matrix)
  }else{
    b_matrix <- data.table(permutation=rep('perm',niter),b_matrix)
  }
  # re-write with permutation column
  write.table(b_matrix,fn_full,quote=FALSE,row.names = FALSE)
  return(NULL)
}

convert_gls_temp_file <- function(){
  # some of the gls permutation runs failed to converge and the function ended before completion so this function is converting the temp file to the final file
  fn <- 'zeudaimoniaFRED15.permutation.gls.ISQF.temp.txt'
  fn_out <- 'zeudaimoniaFRED15.permutation.gls.ISQF.txt'
  do_obs <- FALSE
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
# bootstrap 95% CI on delta - the mean difference between hed and eud coeffs
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
    form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
    Y <- as.matrix(dt[rows, .SD, .SDcols=ycols])
    fitmv <- lm(form, data=dt[rows,])
    hed_matrix[iter,] <- coefficients(fitmv)['zhedonia',]
    eud_matrix[iter,] <- coefficients(fitmv)['zeudaimonia',]
    rows <- sample(1:nrow(dt),replace=TRUE)
  }
  b_matrix <- data.table(I=rep(1:niter,2), Type=rep(c('hedonic','eudaimonic'),each=niter), rbind(hed_matrix,eud_matrix))
  
  if(write_it==TRUE){write.table(b_matrix,fn,quote=FALSE,sep='\t',row.names=FALSE)}
  return(NULL)
}

contrast_coefficients <- function(dt){
  # dt is a matrix niter * p matrix of beta coefficients or raw data
  # compute contrast coefficients AND compute mean of these
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  # cc <- expression_coefficient(ycols) # not needed
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

type.1.error <- function(x){
  # obs statistic must be in first cell
  # usage:  t <- apply(res, 2, type.1.error)
  return(length(which(x <= 0.05))/length(x))
}

effect_size <- function(res){
  # change in mean expression give 4sd change in hedonic score
  # recover numbers from Fredericksonddd
  if('IL6' %in% colnames(res)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  obs <- t(as.matrix(res[I==1,.SD, .SDcols=ycols]))
  obs2 <- (2^(4*apply(obs,2,mean)) - 1)*100
  
  names(obs2) <- c('zhedonia','zeudaimonia')
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
  dt.cc <- contrast_coefficients(copy(dt))
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv <- lm(form, data=dt.cc)
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
    dt.cc <- contrast_coefficients(copy(dt))
    fitmvr <- lm(form, data=dt.cc)
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

smart_t_stats <- function(cc){
  # cc must be the contrast coefficients and not raw regression coefficients
  if('IL6' %in% colnames(cc)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  cc[,ME:=apply(.SD,1,mean),.SDcols=ycols, by=Type]
  
  # naive bootstrap t-test
  y <- as.numeric(cc[1,.SD, .SDcols=ycols])
  p <- length(y)
  ybar <- numeric(200)
  for(iter in 1:200){
    ybar[iter] <- mean(y[sample(1:p, replace=TRUE)])
  }
  se.naive <- sd(ybar)
  
  me <- cc[I==1,ME,by=Type] # observed values
  se <- cc[,.(se=sd(ME)),by=Type] # about 10X higher than naive SE
  t.table <- merge(me,se,by='Type')
  t.table[, t:=abs(ME)/se]
  t.table[, p.value := 2*pt(t,df=(p-1),lower.tail = FALSE)]
  delta <- me[Type=='hedonic',ME] - me[Type=='eudaimonic',ME]
  delta.se <- sqrt(se[Type=='hedonic',se]^2 + se[Type=='eudaimonic',se]^2)
  delta.t <- abs(delta)/delta.se
  delta.p <- 2*pt(delta.t,df=(p-1),lower.tail = FALSE)
  t.table <- rbind(t.table, data.table(Type='delta',ME=delta,se=delta.se,t=delta.t, p.value=delta.p))
  return(t.table)
}

naive_t_stats <- function(cc){
  # compute t-statistic of regulation as in Frederickson 2013
  if('IL6' %in% colnames(cc)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  cc[,ME:=apply(.SD,1,mean),.SDcols=ycols, by=Type]
  
  # y <- cc[I==1,.SD, .SDcols=ycols] # check order of hed and eud
  
  # transpose to get contrast coefficients in columns
  y <- data.table(t(cc[I==1,.SD, .SDcols=ycols]))
  setnames(y, c('zhedonia', 'zeudaimonia'))
  p <- nrow(y)
  
  iters <- 200 # following the original
  yboot <- matrix(0,nrow=iters,ncol=2)
  rows <- 1:p
  for(iter in 1:iters){
    yboot[iter,] <- apply(y[rows],2,mean)
    rows <- sample(1:p,replace=TRUE)
  }
  colnames(yboot) <- c('zhedonia','zeudaimonia')
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
    b.hed <- cc[Type=='hedonic', ME]
    b.eud <- cc[Type=='eudaimonic', ME]
  }
  if(statistic=='t'){
    b.hed <- cc[Type=='hedonic', t]
    b.eud <- cc[Type=='eudaimonic', t]
    delta <- delta/sqrt(cc[Type=='hedonic',SD]^2/p + cc[Type=='eudaimonic',SD]^2/p)
  }
  
  niter <- length(b.hed)
  
  # is hedonic up regulated? This is two-sided despite the question
  zhedonia <- length(which(abs(b.hed)>=abs(b.hed[1])))/niter # permutation p
  # is eudamonic down regulated?
  zeudaimonia <- length(which(abs(b.eud)>=abs(b.eud[1])))/niter # permutation p
  # is there a difference in regulation?
  delta <- length(which(abs(delta)>=abs(delta[1])))/niter # permutation p
  res <- c(zhedonia,zeudaimonia,delta)
  names(res) <- c('zhedonia','zeudaimonia','delta')
  return(res)
}

permutation_lambda <- function(res){
  # res is the table of Wilk's lambda for each iteration
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

gls_with_correlated_error <- function(dt){
  # dt is a raw data set
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()

  # create matrix with scaled y
  #y.scale <- scale(dt[,.SD,.SDcols=ycols])
  Y <- scale(contrast_coefficients(dt[,.SD,.SDcols=ycols]))
  dts <- cbind(dt[,.SD,.SDcols=xcols], Y)
  dts[,subject:=factor(.I)]
  
  # wide to long
  dtlong <- melt(dts,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  
  # this replicates the analysis of Fredrickson et al.
  # subject is the grouping factor and this is identified in the correlation statement
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  # summary(fit1)
  #              Value Std.Error    t-value p-value
  # zhedonia     0.08575  0.122458   0.700267  0.4838
  # zeudaimonia -0.51069  0.125713  -4.062329  0.0000
  
  return(fit1)
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
  dts <- cbind(dt[,.SD,.SDcols=xcols], Y)
  dts[,subject:=factor(.I)]
  
  # wide to long
  dtlong <- melt(dts,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
  dtlong[,gene:=factor(gene)]
  
  # this replicates the analysis of Fredrickson et al.
  # subject is the grouping factor and this is identified in the correlation statement
  form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
  fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=glsControl(msMaxIter = 500, msVerbose = FALSE))
  summary(fit1)
  #              Value Std.Error    t-value p-value
  # zhedonia     0.08575  0.122458   0.700267  0.4838
  # zeudaimonia -0.51069  0.125713  -4.062329  0.0000
  anova(fit1)
  Anova(fit1)
  
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


bootstrap_gls <- function(dt,xcols,ycols,niter=1000, write_it=FALSE,fn){
  # hybrid of bootstrap and gls_with_correlated_error
  # basically use the bootstrap data to get a distribution of b_hedonia and b_zeudaimonia
  
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  
  p <- length(ycols)
  zcols <- c('zhedonia','zeudaimonia')
  xcols2 <- setdiff(xcols, zcols)
  rows <- 1:nrow(dt) # observed on first iter and permuted after
  res <- matrix(0,nrow=niter,ncol=2)
  colnames(res) <- c('zhedonia','zeudaimonia')
  for(iter in 1:niter){
    # create matrix with scaled y
    #y.scale <- scale(dt[,.SD,.SDcols=ycols])
    dt.boot <- dt[rows,]
    Y <- scale(contrast_coefficients(dt.boot[,.SD,.SDcols=ycols]))
    dt.boot <- cbind(dt.boot[,.SD,.SDcols=xcols],Y)
    dt.boot[,subject:=factor(.I)]
    
    # wide to long
    dtlong <- melt(dt.boot,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
    dtlong[,gene:=factor(gene)]
    
    # this replicates the analysis of Fredrickson et al.
    form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
    fit1 <- gls(form, data=dtlong, method='ML', correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene), control=lmeControl(msMaxIter = 500, msVerbose = FALSE))
    res[iter,] <- coefficients(summary(fit1))[c('zhedonia','zeudaimonia')]
    
    # resample rows
    rows <- sample(1:nrow(dt),replace=TRUE)
  }
  
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

explore_MRCE <- function(){
  # multivariate regression with correlated response
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  if('IL6' %in% colnames(dt)){year <- 2013}else{year <- 2015}
  ycols <- c(pro_inflam_genes(year),antibody_genes(),ifn_genes())
  xcols <- get_xcols()
  
  Xobs <- dt[,.SD,.SDcols=xcols]
  Yobs <- as.matrix(contrast_coefficients(dt[,.SD,.SDcols=ycols]))
  
  # convert factors to numeric
  Xobs[,male:=as.numeric(as.character(male))]
  Xobs[,white:=as.numeric(as.character(white))]
  Xobs[,smoke:=as.numeric(as.character(smoke))]
  Xobs <- as.matrix(Xobs)
  fit.mrce <- mrce(X=Xobs,Y=Yobs,method='single',lam1=10^(-1.5), lam2=10^(-0.5))
  data.table(predictor=colnames(Xobs),bean_beta=apply(fit.mrce$Bhat[,],1,mean))
  
  fit.mrce$mx
  
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
    fit <- lm(y~x1 + x2,, data=fd)
    fdres[iter,] <- c(coefficients(fit)[c('x1','x2')],vif(fit))
  }
  colnames(fdres) <- c('x1','x2')
  qplot(x=x1,y=x2,data=fdres)
}
