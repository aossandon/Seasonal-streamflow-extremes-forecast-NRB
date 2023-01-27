##Calibration_Nst_GEV_1cov_sgst_best_1mlead
##last update: 16/01/23

##before running the script
##### go to >>Session >> Set Working Directory >> To source file location


#### CLEAR MEMORY
rm(list=ls())
#### Prevent Warnings
options(warn=-1)
#### Clear the R console
cat("\f")

# Inputs ------------------------------------------------------------------

LeadT=3#enter the lead time 0-3 for month lead time. For stationary model, enter 4

# Read data ---------------------------------------------------------------

if (LeadT==4){
  folder1="GEV_stationary"
  pred="St"
} else if (LeadT>=0 & LeadT<4) {
  N_sea=paste(LeadT,"mth",sep="")
  folder1=paste("GEV_",LeadT,"mth",sep="")
  if (LeadT<=1){
    pred="EPCT_SgSt"
    ind_covar=c(5)
  }else{
    pred="EPCT_PWPR_SgSt"
    ind_covar=c(5,4)
  }
}else {
  print(warning("choose a valid option please"))
}

time_elapsed=0
ptm1 = proc.time()
mainDir=getwd()

knitr::opts_knit$set(root.dir = mainDir)
setwd(mainDir)
suppressPackageStartupMessages(source("Library.R"))
source("Library.R")
# Load Packages
load_packages(c("tiff",'parallel','rstan','ismev','dplyr',"maps","ggplot2","extRemes","leaps","MPV","fields","latex2exp","copula","scales","RColorBrewer","reshape2","data.table","MASS","rlist","sf","raster","MVN","knitr","Hmisc","sm","parallel","superdiag","matrixcalc","ggspatial","foreach","doParallel","evd"),quietly=TRUE)
dir_create(paste("Results",folder1,sep = "/"))
dir_r=paste(mainDir,"Results",folder1,sep = "/")
mday=1

Year=1978:2018



setwd(mainDir)
data1=list.load(paste("data_paper_Qex_",mday,"day.rds",sep=""))
data=data1$JuAU
data=data[data$Year<=Year[length(Year)],]
data[,2:ncol(data)]=(data[,2:ncol(data)])/10000
info=data1$info
info$name[1]=colnames(data)[2]


if (LeadT>=0 & LeadT<4) {
  covar1=read.delim(paste("covariates_1978_2018_Qdaily_",N_sea,".txt",sep=""), header = TRUE, sep = "\t", dec = ".")
}
# colnames(covar1)[1:5]=c(colnames(covar1)[1:3],"PWPR","EPCT")
# write.table(covar1[,1:5], file = paste("covariates_1978_2018_Qdaily_",N_sea,".txt",sep=""),col.names = TRUE,row.names=FALSE,sep="\t")

##find years with missing values

ind_o=which(is.na(data$Garudeshwar)==F & is.na(data$Mandleshwar)==F & is.na(data$Handia)==F & is.na(data$Hoshangabad)==F & is.na(data$Sandiya)==F)
data_f=data[ind_o,]
if (LeadT>=0 & LeadT<4) {
  covar1_f=covar1[ind_o,]
}
# set the data for stan model ---------------------------------------------

# sampler iterations
iterations = 12000
warmup = 6000
rng_seed = 1111
N=nrow(data_f)
S=nrow(info)
if (LeadT<=1) {
  nam_par=c("alpha_loc0","alpha_loc1","Log_scale0","shape","cop_mat")
  gev_data = list(N=N,S=S,Nc=2, y=as.matrix(data_f[,2:ncol(data_f)]),covar1=matrix(rep(covar1_f[,ind_covar[1]],5),N,S))
  
}else if(LeadT<=3){
  nam_par=c("alpha_loc0","alpha_loc1","alpha_loc2","Log_scale0","shape","cop_mat")
  gev_data = list(N=N,S=S,Nc=2, y=as.matrix(data_f[,2:ncol(data_f)]),covar1=matrix(rep(covar1_f[,ind_covar[1]],5),N,S),covar2=matrix(rep(covar1_f[,ind_covar[2]],5),N,S))
}else{
  nam_par=c("alpha_loc0","Log_scale0","shape","cop_mat")
  #setup data
  gev_data = list(N=N,S=S,Nc=2, y=as.matrix(data_f[,2:ncol(data_f)]))
}


# running the stan model --------------------------------------------------
nchain=3
thin1=(iterations-warmup)/1000
if (thin1<=1){thin1=1}
time_elapsed=0
ptm = proc.time()
setwd(mainDir)
if (LeadT<=1) {
  model_file = "gev_multi_cop_1cov_Sgst.stan"
}else if(LeadT<=3){
  model_file = "gev_multi_cop_2cov_Sgst.stan"
}else{
  model_file = "gev_multi_cop_St.stan"
}
gev_model = stan(model_file, data = gev_data, chains = 0)
# ,max_treedepth=14
stanfit = stan(fit = gev_model, data = gev_data,
               chains =nchain, iter=iterations,
               seed=rng_seed, warmup = warmup,thin=thin1,control=list(max_treedepth=13,adapt_delta=0.9),cores=nchain)
aa1=rstan::extract(stanfit)
sprintf("elapsed time: %s minutes",round((proc.time() - ptm)[3]/60,2))

# posterior plots and performance metrics ---------------------------------



sampler_params <- get_sampler_params(stanfit, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

#the loo package to carry out Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) for purposes of model checking
loo1=loo(stanfit, save_psis = TRUE)
setwd(dir_r)
val1=rbind(loo1$estimates,c(max(summary(stanfit)$summary[,10],na.rm =T),min(summary(stanfit)$summary[,10],na.rm =T)))
rownames(val1)=c(rownames(loo1$estimates),"rhat")
colnames(val1)=colnames(loo1$estimates)
write.table(val1, file = "loo_summary_wo_cop.txt",sep="\t")
loo1
print(round(val1,4))
plot(loo1)
p1.base <- recordPlot()
invisible(dev.off())

jpeg("loo_plot_Byeasian_wo_cop.jpg", width =6, height = 4, units = 'in', res = 300)
p1.base
dev.off()
p1.base


# traceplots for each station ---------------------------------------------

#Garudeshwar
i=1
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), las = 1)
dev.off()

#Mandleshwar

i=2
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), las = 1)
dev.off()
#Handia
i=3
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), las = 1)
dev.off()
#Hoshangabad
i=4
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), las = 1)
dev.off()
#Sandiya
i=5
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = c(paste(nam_par[1:(length(nam_par)-1)],"[",i,"]",sep="")), las = 1)
dev.off()






# get the correlations parameters matrix ----------------------------------



aa=aa1
rm(stanfit,aa1)
tmp=aa$cop_mat
corr=data.frame(matrix(data=NA,nrow = dim(tmp)[1],ncol=dim(tmp)[2]*(dim(tmp)[2]-1)/2))
tmp2=c()
for (i in 1:dim(tmp)[2]) {
  for (j in 2:dim(tmp)[2]) {
    if (i < j){
      tmp2=c(tmp2,paste("rho_",i,j,sep=""))
    }
    
  }
  
}
colnames(corr)=tmp2
for (i in 1:dim(tmp)[1]) {
  tmp1=tmp[i,,]
  corr[i,]=tmp1[lower.tri(tmp1,diag = FALSE)]
}

aa$corr=corr
##save file with posterior samples of parameters
list.save(aa, paste("Par_",pred,".rds",sep=""))


# generate posterior samples of streamflow --------------------------------

## Bayesian with copula using the  inference functions for margins (IFM) approach

ptm = proc.time()
numCores <- detectCores()-1
NN=nrow(aa$shape)
NN1=nrow(data)

registerDoParallel(numCores)
res1=foreach (k=1:NN1, .packages=c("evd","copula")) %dopar% {
  tmp2=matrix(NA,NN,nrow(info))
  for (j in 1:NN) {
    n_sam=rCopula(1,normalCopula(as.numeric(aa$corr[j,]),dim=nrow(info),dispstr='un'))
    for (i in 1:nrow(info)) {
      shape=aa$shape[j,i]
      scale=exp(aa$Log_scale0[j,i])
      if (LeadT==4) {
      loc=aa$alpha_loc0[j,i]
      }else if(LeadT<=1){
        loc=aa$alpha_loc0[j,i]+covar1[k,ind_covar[1]]*aa$alpha_loc1[j,i]
      }else{
        loc=aa$alpha_loc0[j,i]+covar1[k,ind_covar[1]]*aa$alpha_loc1[j,i]+covar1[k,ind_covar[2]]*aa$alpha_loc2[j,i]
      }
      tmp2[j,i]=qgev(n_sam[i],loc=loc,scale = scale,shape = shape)
      if (tmp2[j,i]<0){tmp2[j,i]=0}
    }
  }
  
  tmp2
}
stopImplicitCluster()
post_y=array(unlist(res1),dim = c(NN,nrow(info),NN1))
dimnames(post_y)=list(1:NN,as.character(info$name),data[,1])
sprintf("elapsed time: %s minutes",round((proc.time() - ptm)[3]/60,2))


## bayesian without copula

time_elapsed=0
ptm = proc.time()
NN=nrow(aa$shape)
NN1=nrow(data)

registerDoParallel(numCores)
res1=foreach (k=1:NN1, .packages=c("evd","copula")) %dopar% {
  tmp2=matrix(NA,NN,nrow(info))
  for (j in 1:NN) {
    for (i in 1:nrow(info)) {
      n_sam=runif(1,0,1)
      shape=aa$shape[j,i]
      scale=exp(aa$Log_scale0[j,i])
      if (LeadT==4) {
        loc=aa$alpha_loc0[j,i]
      }else if(LeadT<=1){
        loc=aa$alpha_loc0[j,i]+covar1[k,ind_covar[1]]*aa$alpha_loc1[j,i]
      }else{
        loc=aa$alpha_loc0[j,i]+covar1[k,ind_covar[1]]*aa$alpha_loc1[j,i]+covar1[k,ind_covar[2]]*aa$alpha_loc2[j,i]
      }
      tmp2[j,i]=qgev(n_sam,loc=loc,scale = scale,shape = shape)
      if (tmp2[j,i]<0){tmp2[j,i]=0}
    }
  }
  
  tmp2
}
stopImplicitCluster()
post2_y=array(unlist(res1),dim = c(NN,nrow(info),NN1))
dimnames(post2_y)=list(1:NN,as.character(info$code),data[,1])

sprintf("elapsed time: %s minutes",round((proc.time() - ptm)[3]/60,2))

##convert flow data to 10^(-3)*m3/s
data[,2:6]=data[,2:6]*10
post_y=post_y*10
post2_y=post2_y*10


##compute some deterministic metrics
Val=matrix(NA,5,4)
colnames(Val)=c("R_BHM","%BIAS_BHM","R_woCOp","%BIAS_woCOp")
rownames(Val)=info$name
for (j in 1:5) {
  tmp=apply(post_y[,j,],2,mean)
  Val[j,1]=round(cor.test(tmp,data[,j+1],na.action("na.omit"),method="pearson")$estimate,3)
  Val[j,2]=round((mean(tmp,na.rm=T)-mean(data[,j+1],na.rm=T))/mean(data[,j+1],na.rm=T)*100,3)
  
  
  tmp=apply(post2_y[,j,],2,mean)
  Val[j,3]=round(cor.test(tmp,data[,j+1],na.action("na.omit"),method="pearson")$estimate,3)
  Val[j,4]=round((mean(tmp,na.rm=T)-mean(data[,j+1],na.rm=T))/mean(data[,j+1],na.rm=T)*100,3)
}
print(Val)

write.table(Val, file = paste("DeterministicMetrics_",N_sea,".txt",sep=""),col.names = TRUE,row.names=FALSE,sep="\t")



# plot of posterior streamflow at each gauge with Copula ------------------------------

jpeg(paste("TS_IFM_",info$name[1],"_",info$name[5],".jpg",sep=""), width =9, height = 12, units = 'in', res = 300)
par(mfrow=c(5,1))
par(mar = c(2, 4, 1.4, 0.3))

i=1
zz=boxplot(split(t(post_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=2
zz=boxplot(split(t(post_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=3
zz=boxplot(split(t(post_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=4
zz=boxplot(split(t(post_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=5
zz=boxplot(split(t(post_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

dev.off()



# plot of posterior streamflow at each gauge without copula ---------------
jpeg(paste("TS_Wocop_",info$name[1],"_",info$name[5],".jpg",sep=""), width =9, height = 12, units = 'in', res = 300)
par(mfrow=c(5,1))
par(mar = c(2, 4, 1.4, 0.3))

i=1
zz=boxplot(split(t(post2_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=2
zz=boxplot(split(t(post2_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=3
zz=boxplot(split(t(post2_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=4
zz=boxplot(split(t(post2_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

i=5
zz=boxplot(split(t(post2_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("Streamflow,"~Q(t) ~"("*10^3*m^3*s^-1*")"), side = 2, line = 2.1, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)

dev.off()


# Basin average flow ------------------------------------------------------

Q_m_gev_cop_IFM=as.data.frame(matrix(data=NA,nrow=NN,ncol=NN1))
Q_m_gev_bay=Q_m_gev_cop_IFM


for (i in 1:NN) {
  tmp1=matrix(data=NA,nrow=NN1,ncol=nrow(info))
  tmp_1=tmp1
  for (j in 1:nrow(info)) {
    
    tmp1[,j]=post_y[i,j,]/info$area[j]*(1000*86400)/(1000)^2*1000#mm/day
    tmp_1[,j]=post2_y[i,j,]/info$area[j]*(1000*86400)/(1000)^2*1000#mm/day
    
  }
  #with copula
  Q_m_gev_cop_IFM[i,]=apply(tmp1,1,mean)
  #without copula
  Q_m_gev_bay[i,]=apply(tmp_1,1,mean)
  
}

tmp1=matrix(data=NA,nrow=NN1,ncol=nrow(info))
for (j in 1:nrow(info)) {
  tmp1[,j]=data[,j+1]/info$area[j]*(1000*86400)/(1000)^2*1000#mm/day
}
tmp2=apply(tmp1,1,mean)

jpeg("basin_Average_flow_TwoMethods.jpg", width =6, height = 7, units = 'in', res = 300)
ran=2.3
par(mfrow=c(2,1))
par(mar = c(2, 3, 1.4, 0.4))

zz=boxplot(split(t(Q_m_gev_cop_IFM),data[,1]),plot=F,cex=1.0,range=ran)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="",xlim=c(1.5,40.5),cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext("Specific Streamflow (mm/day)", side = 2, line = 1.9, cex = 1,font = 1)
lines(z1,tmp2,lty=1,lwd=1,col="#0b5394")
points(z1,tmp2,lwd=2,col="#0b5394")
title("with copula", cex.main = 1.5,line=0.4, adj = 0.5)



zz=boxplot(split(t(Q_m_gev_bay),data[,1]),plot=F,cex=1.0,range=ran)
zz$names=rep("",length(zz$names))
z1=bxp(zz,xlab="",ylab="",xlim=c(1.5,40.5),cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext("Specific Streamflow (mm/day)", side = 2, line = 1.9, cex = 1,font = 1)
lines(z1,tmp2,lty=1,lwd=1,col="#0b5394")
points(z1,tmp2,lwd=2,col="#0b5394")
title("without copula", cex.main = 1.5,line=0.4, adj = 0.5)

dev.off()



# Save posterior flow samples ---------------------------------------------

pots_f=list(Post_wCop=post_y,Post_woCop=post2_y)
list.save(pots_f, paste("post_Q_cali_",pred,".rds",sep=""))
sprintf("elapsed time: %s minutes",round((proc.time() - ptm1)[3]/60,2))
