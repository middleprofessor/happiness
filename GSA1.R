# Script GSA1.R
# Scripts for Monte Carlo simulation of error rates of GSA methods applied to model of
# Fredrickson et al 2015 data
# Jeffrey A. Walker
# August, 22, 2016
# cleaning of code from Fredrickson_multivariate_lm.peerj.rev1
# Monte Carlo results for manuscript re-run using this and not original code
# requires functions in
  # Fredrickson_peerj.rev2.R
  # GSA_methods.R

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

start_here <- function(){
  #do_it()
  res_table <- table_error_rates()
  error_table <- clean_error_table(copy(res_table))
  print(error_table)
  error_figure(res_table[model!='gls.un',])
  write.table(error_table,'error_table.txt',row.names=FALSE,sep='\t',quote=FALSE)
}

do_it <- function(){
  # this takes about 48 hours on my macbook
  method_list <- c('R2','permF','roast','GA','obrien','gee')
  niter <- 4000
  perms <- 2000
  m_array <- c(10,30,52)
  for(m in m_array){
    res1 <- simulate_it_1(niter=niter,beta=-9999,method_list=method_list,n_array=-9999, m_array=m, perms=perms,do_power=FALSE,write_it=TRUE)
    res2 <- simulate_it_1(niter=niter,beta=-9999,method_list=method_list,n_array=-9999, m_array=m, perms=perms,do_power=TRUE,write_it=TRUE)
  }
  
  # Because of the time to model GLS, GLS was run separately and on multiple computers using
  # versions of this
  method_list <- c('gls') # office macs
  niter <- 200 # 4000, 2000, 1000 for m = 10, 30, 52
  perms <- 2000
  m_array <- c(52)
  for(m in m_array){
    res1 <- simulate_it_1(niter=niter,beta=-9999,method_list=method_list,n_array=-9999, m_array=m, perms=perms,do_power=FALSE,write_it=TRUE)
    #res2 <- simulate_it_1(niter=niter,beta=-9999,method_list=method_list,n_array=-9999, m_array=m, perms=perms,do_power=TRUE,write_it=TRUE)
  }
  
  # an attempt to run an unstructured matrix
  method_list <- c('gls.un') # office macs
  niter <- 200 # 4000, 2000, 1000 for m = 10, 30, 52
  perms <- 2000
  m_array <- c(52)
  for(m in m_array){
    #res1 <- simulate_it_1(niter=niter,beta=-9999,method_list=method_list,n_array=-9999, m_array=m, perms=perms,do_power=FALSE,write_it=TRUE)
    res2 <- simulate_it_1(niter=niter,beta=-9999,method_list=method_list,n_array=-9999, m_array=m, perms=perms,do_power=TRUE,write_it=TRUE)
  }
  

}


table_error_rates <- function(){ # and Figure!
  # script to collect the results files of the simulation in function do_it and
  # table the error types and generate the figure for the manuscript
  
  res1 <- data.table(NULL)
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.WOUI.txt',header=TRUE,sep='\t'))) 
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.JKKC.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.SMNX.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.BILN.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.BMRD.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.HQRH.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.LZAH.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.GALD.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.ROGL.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.IPXZ.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.RPYK.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.QGZB.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.NMYM.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.QCWH.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.LGIX.txt',header=TRUE,sep='\t')))
  res1 <- rbind(res1,data.table(read.table('gsa_sim_1.typeI.ETOD.txt',header=TRUE,sep='\t')))
  
  res2 <- data.table(NULL)
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.YOGT.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.RPCV.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.ELSO.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.DTEO.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.UULC.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.DMYB.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.ONKA.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.HFZQ.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.XGGD.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.IPZD.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.UWHY.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.NRNJ.txt',header=TRUE,sep='\t')))
  res2 <- rbind(res2,data.table(read.table('gsa_sim_1.typeII.AMOB.txt',header=TRUE,sep='\t')))
  
  # use only first 1000 rows of gls m=52
  inc <- which(res1[,model]=='gls' & res1[,m]==52)
  diffinc <- setdiff(1:nrow(res1),inc)
  part1 <- res1[diffinc]
  keep <- min(1000,length(inc))
  part2 <- res1[inc[1:keep]]
  res1 <- rbind(part1,part2)

  inc <- which(res2[,model]=='gls' & res2[,m]==52)
  diffinc <- setdiff(1:nrow(res2),inc)
  part1 <- res2[diffinc]
  keep <- min(1000,length(inc))
  part2 <- res2[inc[1:keep]]
  res2 <- rbind(part1,part2)
  
  alpha <- 0.05
  trueb <- 0.06716244 # the value of beta used in the simulation with an effect
  t1 <- res1[,.(N.I=.N,TypeI=length(which(prob<=alpha))/.N),by=.(model,n,m)] #MAE = mean absolute error
  # for type II
  t2 <- res2[,.(N.II=.N,Power=length(which(prob<=alpha & b>0))/.N, S=length(which(prob<=alpha & b<0))/length(which(prob<=alpha))),by=.(model,n,m)]
  t3 <- res2[prob<=alpha & b>0,.(ER=mean(b/trueb)),by=.(model,n,m)]
  
  res_table <- merge(t1,t2,by=c('model','n','m'),all=TRUE)
  res_table <- merge(res_table,t3,by=c('model','n','m'),all=TRUE)
  res_table <- orderBy(~n + m, data=res_table)
  #res_table
  
  # compute bias
  res_table[m==52,.(res_table[m==52 & model=='GA',Power]/Power),by=model]
  t4 <- res2[,.(b_hat=round((mean(b)-trueb)/trueb*100,4)),by=.(model,n,m)] # test for bias
  return(res_table)
  
}

gls_adjusted_stats <- function(res_table){
  # ust this to find the adjusted alpha to make the GLS type I error == 0.05 to see how this affects power
  alpha <- 0.00015
  t1 <- res1[,.(N.I=.N,TypeI=length(which(prob<=alpha))/.N),by=.(model,n,m)] #MAE = mean absolute error
  t2 <- res2[,.(N.II=.N,Power=length(which(prob<=alpha & b>0))/.N, S=length(which(prob<=alpha & b<0))/length(which(prob<=alpha))),by=.(model,n,m)]
  t3 <- res2[prob<=alpha & b>0,.(ER=mean(b/trueb)),by=.(model,n,m)]
  
  alt_table <- merge(t1,t2,by=c('model','n','m'),all=TRUE)
  alt_table <- merge(alt_table,t3,by=c('model','n','m'),all=TRUE)
  alt_table <- orderBy(~n + m, data=alt_table)
  alt_table[model=='gls']
  # for type II
  # how much power does GLS have when type I = 0.05?
  #m=10, alpha = 0.0155, GLS type I = 0.05000, Power = 0.156500
  #m=30, alpha = 0.0037, GLS type I = 0.05000, Power = 0.1445000
  #m=52, alpha = 0.00015, GLS type I = 0.05000, Power = 0.11500
}

clean_error_table <- function(error_table){
  # input is res_table
  error_table[,TypeI:=round(TypeI,3)]
  error_table[,Power:=round(Power,2)]
  error_table[,S:=round(S,3)]
  error_table[,ER:=round(ER,1)]
  return(error_table)
}

error_figure <- function(res_table){
  # res_table is the output from 
  gg <- ggplot(data=res_table,aes(x=m,y=TypeI,color=model))
  gg <- gg + geom_point(aes(shape=model))
  gg <- gg + scale_shape_manual(values=c(0,1,2,5,16,17,15)) 
  gg <- gg + geom_line()
  gg <- gg + ggtitle('A')
  gg <- gg + labs(x = 'genes',y = 'Type I')
  gg <- gg + scale_colour_brewer(palette = "Dark2")
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position='none',plot.margin=unit(x=c(0,0,0,0),'cm'))
  gg1 <- gg
  gg
  
  gg <- ggplot(data=res_table,aes(x=m,y=Power,color=model))
  gg <- gg + geom_point(aes(shape=model))
  gg <- gg + scale_shape_manual(values=c(0,1,2,5,16,17,15)) 
  gg <- gg + geom_line()
  gg <- gg + ggtitle('B')
  gg <- gg + labs(x = 'genes')
  gg <- gg + scale_colour_brewer(palette = "Dark2")
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position='none',plot.margin=unit(x=c(0,0,0,0),'cm'))
  gg2 <- gg
  gg
  
  gg <- ggplot(data=res_table,aes(x=m,y=S,color=model))
  gg <- gg + geom_point(aes(shape=model))
  gg <- gg + scale_shape_manual(values=c(0,1,2,5,16,17,15)) 
  gg <- gg + geom_line()
  gg <- gg + ggtitle('C')
  gg <- gg + labs(x = 'genes',y = 'Type S')
  gg <- gg + scale_colour_brewer(palette = "Dark2")
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position='none',plot.margin=unit(x=c(0,0,0,0),'cm'))
  gg3 <- gg
  gg
  
  gg <- ggplot(data=res_table,aes(x=m,y=ER,color=model))
  gg <- gg + geom_point(aes(shape=model))
  gg <- gg + scale_shape_manual(values=c(0,1,2,5,16,17,15)) 
  gg <- gg + geom_line()
  gg <- gg + ggtitle('D')
  gg <- gg + labs(x = 'genes',y = 'Exaggeration Ratio')
  gg <- gg + scale_colour_brewer(palette = "Dark2")
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.position='none',plot.margin=unit(x=c(0,0,0,0),'cm'))
  gg4 <- gg
  gg
  
  fig_name <- paste('Fig_errors.pdf',sep='')
  pdf(fig_name,paper='special',onefile=FALSE,width=5.5,height=5.5)
  #  postscript('fig_01.eps',horizontal=FALSE,onefile=FALSE,paper='special',height=3,width=6.5)
  showtext.begin()
  print(gg1)
  print(gg2)
  print(gg3)
  print(gg4)
  grid.arrange(gg1,gg2,gg3,gg4,ncol=2,nrow=2)
  showtext.end()
  dev.off()
  
  # print a legend
  gg <- ggplot(data=res_table,aes(x=m,y=ER,color=model, shape=model))
  gg <- gg + geom_point()
  gg <- gg + geom_line()
  gg <- gg + ggtitle('D')
  gg <- gg + labs(y = 'Exaggeration Ratio')
  gg <- gg + scale_colour_brewer(palette = "Dark2",labels=c('Fga','gee','obrien','Fpun','R2','roast','gls'))
  gg <- gg + scale_shape_manual(values=c(0,1,2,5,16,17,15),labels=c('Fga','gee','obrien','Fpun','R2','roast','gls')) 
  gg <- gg + theme_bw() + theme(axis.title=element_text(size=10),axis.text=element_text(size=8),plot.title=element_text(hjust=0),strip.text=element_text(size=8),legend.title=element_blank(),legend.text=element_text(size=8),legend.position='bottom',plot.margin=unit(x=c(0,0,0,0),'cm'))
  #gg <- gg + guides(color = guide_legend(label.position = "bottom"))
  gg
  
  legend <- ggplot_gtable(ggplot_build(gg))$grobs
  #dev.new()
  #pushViewport(plotViewport(rep(1, 4)))
  #grid.draw(legend[[8]])
  
  
  fig_name <- paste('Fig_errors_legend.pdf',sep='')
  pdf(fig_name,paper='special',onefile=FALSE,width=4.75,height=.5)
  #  postscript('fig_01.eps',horizontal=FALSE,onefile=FALSE,paper='special',height=3,width=6.5)
  showtext.begin()
  grid.draw(legend[[8]])
  showtext.end()
  dev.off()
  
  
}

bias_plot <- function(){
  # need res2 from table_error_rates
  trueb <- 0.06716244 # the value of beta used in the simulation with an effect
  olsdata <- res2[model=='obrien' & m==52,]
  glsdata <- res2[model=='gls' & m==52,]
  N.ols <- nrow(olsdata)
  N.gls <- nrow(glsdata)
  niter <- 10000
  iters <- as.integer(runif(niter,100,N.gls))
  b.ols <- matrix(0,nrow=niter,ncol=2)
  colnames(b.ols) <- c('iters','b')
  b.gls <- copy(b.ols)
  for(iter in 1:niter){
    n <- iters[iter] # number of rows to sample
    inc <- sample(1:N.ols,n)
    b.ols[iter,] <- c(n,mean(olsdata[inc,b]))
    inc <- sample(1:N.gls,n)
    b.gls[iter,] <- c(n,mean(glsdata[inc,b]))
  }
  dt <- rbind(data.table(method='ols',b.ols),data.table(method='gls',b.gls))
  gg <- ggplot(data=dt,aes(x=iters,y=b,color=method))
  gg <- gg + geom_point(alpha = .1)
  gg
  
  #bias
  (mean(b.gls[,'b']) - trueb)/trueb
  res2[,.(b_hat=round((mean(b)-trueb)/trueb*100,4)),by=.(model,n,m)] # test for bias
  
}

simulate_it_1 <- function(niter=2000,beta=-9999,method_list=c('permF','obrien','roast','GA','gee'),n_array=-9999, m_array=-9999, perms=2000,do_power=FALSE,write_it=FALSE){
  # simulate FRED15 (Cole et al. 2015)
  # result statistics are for eudaimonia but hedonia is kept to look at correlation in estimates when effect is zero for both
  # do_power - makes the eudaimonia effect equal to beta
  # if beta==-9999 then beta is set to the value estimated by ols using FRED15
  # if n_array=-9999 then n is set the value of FRED15
  # if m_array=-9999 then m is set to the value of FRED15


  # get the empirical data
  fn <- 'cole2_clean.txt'
  dt <- read_file(fn,year=2015)
  dt[,zhedonia:=scale(zhedonia)]
  dt[,zeudaimonia:=scale(zeudaimonia)]
  dt <- contrast_coefficients(dt)  # convert to CTRA response
  xcols <- get_xcols()
  ycols <- c(pro_inflam_genes(year=2015),antibody_genes(),ifn_genes())
  zcols <- 'zeudaimonia'
  dt[,male:=as.integer(as.character(male))]
  dt[,white:=as.integer(as.character(white))]
  dt[,alcohol:=as.integer(as.character(alcohol))]
  dt[,smoke:=as.integer(as.character(smoke))]
  X <- dt[,.SD,.SDcols=xcols]
  Y <- scale(dt[,.SD,.SDcols=ycols])
  dt <- cbind(X,Y)
  Ry <- cor(Y)
  Rx <- cor(X)
  # get mean and variance of expression levels of hedonia and eudaimonia
  form <- formula(paste('Y',paste(xcols,collapse='+'),sep='~'))
  fitmv <- lm(form, data=dt)
  b.hedm <- coefficients(fitmv)['zhedonia',]
  b.eudm <- coefficients(fitmv)['zeudaimonia',]
  if(beta==-9999){beta.mu <- abs(mean(b.eudm))}
  # beta.mu <- 0.06716244
  beta.sd <- sd(b.eudm)

  n_data <- nrow(dt)
  m_data <- ncol(Y)
  if(n_array[1]==-9999){n_array <- n_data}
  if(m_array[1]==-9999){m_array <- m_data}
  param_matrix <- expand.grid(n_array, m_array)
  colnames(param_matrix) <- c('n','m')
  
  #output file
  code <- sample(LETTERS,4,replace=TRUE)
  if(do_power==TRUE){error_type='typeII'}else{error_type='typeI'}
  fn_out <- paste('gsa_sim_1',error_type,paste(code,collapse=''),'txt',sep='.')
  
  res <- data.table(NULL)
  for(experiment in 1:nrow(param_matrix)){
    n <- param_matrix[experiment,'n']
    m <- param_matrix[experiment,'m']
    N <- m*n
    tvalue <- numeric(m)
    rycols <- paste('Y',1:m,sep='')
    
    for(iter in 1:niter){
      
      #  simulated vector of causal effects of eudaimonia on the m expression levels
      # re-center and scale so that the coefficients have precisely the mean and sd specified
      beta_vec <- rnorm(m,mean=beta.mu,sd=beta.sd)
      beta_vec <- (scale(beta_vec))[,1]*beta.sd + beta.mu
      
      X <- rmvnorm(n=n,sigma=Rx)
      colnames(X) <- xcols
      gene_inc <- sample(1:nrow(Ry),m) # subsample the genes
      E <- rmvnorm(n=n,sigma=Ry[gene_inc,gene_inc]) # matrix of expression levels for m genes
      # make Y a function of eudaimonia if do_power==TRUE otherwise Y=E
      if(do_power==TRUE){Y <- X[,'zeudaimonia']%*%t(beta_vec) + t(matrix(sqrt(1-beta_vec^2),nrow=m,ncol=n))*E}else{Y <- E}
      # scale so that true values have precise specification
      X <- scale(X)
      Y <- scale(Y)
      colnames(Y) <- rycols
      rdt <- data.table(cbind(X,Y))
      rdt[,subject:=factor(.I)]
      # wide to long for GEE/GLS
      dtlong <- melt(rdt,id.vars=c('subject',xcols),variable.name='gene',value.name='expression')
      dtlong[,gene:=factor(gene)]
      dtlong <- orderBy(~subject + gene, dtlong)
      # fit ols
      X.dm <- cbind(rep(1,n),X)
      fit.ols <- lm.fit(X.dm,Y)
      b.ols <- mean(fit.ols$coefficients[zcols,])
      b.hed <- mean(fit.ols$coefficients['zhedonia',])
      if('ols' %in% method_list){
        # not fast code using lm.fit - could use the code in obrien.fit to do this
        form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
        fit <- summary(lm(form, data=dtlong)) # used for both .ols and .roast below
        prob <- fit$coefficients['zeudaimonia','Pr(>|t|)']
        res <- rbind(res,data.table(n=n,m=m,model='ols',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if('obrien' %in% method_list){
        # get R, the correlation among Y conditional on all X but zeudaimonia
        obrien_table <- obrien.fit(rdt,xcols,rycols,zcols)
        prob <- obrien_table[,p]
        res <- rbind(res,data.table(n=n,m=m,model='obrien',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if('permt' %in% method_list){ # permutation t
        prob <- permutation_t.fit(rdt,xcols=xcols,ycols=rycols,zcols='zeudaimonia',method='t',perms=perms, write_it=FALSE,fn)
        res <- rbind(res,data.table(n=n,m=m,model='permt',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if('R2' %in% method_list){ # permutation t
        prob <- permutation_t.fit(rdt,xcols=xcols,ycols=rycols,zcols='zeudaimonia',method='R2',perms=perms, write_it=FALSE,fn)
        res <- rbind(res,data.table(n=n,m=m,model='R2',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if('permF' %in% method_list){ # permutation t
        prob <- permutation_F.fit(rdt,xcols=xcols,ycols=rycols,zcols='zeudaimonia',method='resid',perms=perms, write_it=FALSE,fn)
        res <- rbind(res,data.table(n=n,m=m,model='permF',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if('gls' %in% method_list){
        form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
        fit.gls <- gls(form, data=dtlong, method='ML',correlation=corCompSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene))
        # save gls fits
        b.gls <- summary(fit.gls)$tTable['zeudaimonia','Value']
        b.hed <- summary(fit.gls)$tTable['zhedonia','Value']
        prob <- summary(fit.gls)$tTable['zeudaimonia','p-value']
        res <- rbind(res,data.table(n=n,m=m,model='gls',b=b.gls, prob=prob,b.hed=b.hed))
      }
      if('gls.un' %in% method_list){
        form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
        fit.gls <- gls(form, data=dtlong, method='ML', correlation=corSymm(form = ~ 1 | subject), weights=varIdent(form = ~1|gene))
        # save gls fits
        b.gls <- summary(fit.gls)$tTable['zeudaimonia','Value']
        b.hed <- summary(fit.gls)$tTable['zhedonia','Value']
        prob <- summary(fit.gls)$tTable['zeudaimonia','p-value']
        res <- rbind(res,data.table(n=n,m=m,model='gls.un',b=b.gls, prob=prob,b.hed=b.hed))
      }
      if('gee' %in% method_list){
        # gee
        form <- formula(paste('expression~',paste(c('gene',xcols),collapse='+'),sep=''))
        fit.geeglm <- geeglm(form, family=gaussian, data=dtlong,id=subject,waves=gene, corstr='exchangeable', std.err="san.se")
        b.gee <- summary(fit.geeglm)$coefficients['zeudaimonia','Estimate']
        b.hed <- summary(fit.geeglm)$coefficients['zhedonia','Estimate']
        prob <- summary(fit.geeglm)$coefficients['zeudaimonia','Pr(>|W|)']
        se.gee <- summary(fit.geeglm)$coefficients['zeudaimonia','Std.err']
        res <- rbind(res,data.table(n=n,m=m,model='gee',b=b.gee, prob=prob,b.hed=b.hed))
      }
      if('roast' %in% method_list){
        # Roast format
        Yt <- t(Y) # genes as rows
        form <- formula(paste('~',paste(xcols,collapse='+'),sep=''))
        design <- model.matrix(form, data=rdt)
        eudaimonia <- which(colnames(design)=='zeudaimonia')
        roast_res <- roast(y=Yt,design=design,contrast=eudaimonia, nrot=perms)
        prob <- roast_res$p.value['UpOrDown','P.Value']
        res <- rbind(res,data.table(n=n,m=m,model='roast',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if('GA' %in% method_list){ # globalAncova
        Yt <- t(Y) # genes as rows
        form.full <- formula(paste('~',paste(xcols,collapse='+'),sep=''))
        model.dat <- rdt[, .SD, .SDcols=xcols]
        prob <- GlobalAncova(Yt, form.full, model.dat=model.dat, test.terms='zeudaimonia', method='permutation',perm=perms)$test.result['p.perm',1]
        res <- rbind(res,data.table(n=n,m=m,model='GA',b=b.ols, prob=prob,b.hed=b.hed))
      }
      if(write_it==TRUE){
        write.table(res,fn_out,sep='\t',quote=FALSE,row.names=FALSE)
      }
    }
  }
  return(res)
}
