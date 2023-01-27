##Calibration_Nst_Poisson_1cov_sgst_best_dif_month_lead
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

LeadT=4#enter the lead time 0-3 for month lead time. For stationary model, enter 4

# Read data ---------------------------------------------------------------

if (LeadT==4){
  folder1="POI_stationary"
  pred="St"
} else if (LeadT>=0 & LeadT<4) {
  N_sea=paste(LeadT,"mth",sep="")
  folder1=paste("POI_",LeadT,"mth",sep="")
  if (LeadT==0){
    pred="EPCT_SgSt"
    ind_covar=c(5)
  }else if (LeadT==1) {
    pred="EPCT_Nino12_SgSt"
    ind_covar=c(5,2)
  }else if (LeadT==2) {
    pred="EPCT_IOD_SgSt"
    ind_covar=c(5,1)
  }else{
    pred="Nino12_PWPR_SgSt"
    ind_covar=c(2,4)
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
data=data1$OCQmax_JuAU
data=data[data$Year<=Year[length(Year)],]
info=data1$info
info$name[1]=colnames(data)[2]


if (LeadT>=0 & LeadT<4) {
  covar1=read.delim(paste("covariates_1978_2018_Qdaily_",N_sea,".txt",sep=""), header = TRUE, sep = "\t", dec = ".")
} 

##find years with missing values

ind_o=which(is.na(data$Garudeshwar)==F & is.na(data$Mandleshwar)==F & is.na(data$Handia)==F & is.na(data$Hoshangabad)==F & is.na(data$Sandiya)==F)
data_f=data[ind_o,]
if (LeadT>=0 & LeadT<4) {
  covar1_f=covar1[ind_o,]
}


# sampler iterations
iterations = 12000
warmup = 6000
rng_seed = 1111
N=nrow(data_f)
S=nrow(info)
if (LeadT==0) {
  nam_par=c("alpha_loc0","alpha_loc1")#,
  poi_data = list(N=N,S=S, y1=data_f[,2],y2=data_f[,3],y3=data_f[,4],y4=data_f[,5],y5=data_f[,5],covar1=covar1_f[,ind_covar[1]])
}else if(LeadT<=3){
  nam_par=c("alpha_loc0","alpha_loc1","alpha_loc2")
  poi_data = list(N=N,S=S, y1=data_f[,2],y2=data_f[,3],y3=data_f[,4],y4=data_f[,5],y5=data_f[,5],covar1=covar1_f[,ind_covar])
}else{
  nam_par="alpha_loc0"
  #setup data
  poi_data = list(N=N,S=S, y1=data_f[,2],y2=data_f[,3],y3=data_f[,4],y4=data_f[,5],y5=data_f[,5])
}

# running the stan model --------------------------------------------------
nchain=3
thin1=(iterations-warmup)/1000
if (thin1<=1){thin1=1}
time_elapsed=0
ptm = proc.time()
setwd(mainDir)
if (LeadT==0) {
  model_file = "Poisson_multi_1cov.stan"
}else if(LeadT<=3){
  model_file = "Poisson_multi_2cov.stan"
}else{
  model_file = "Poisson_multi_st.stan"
}
poi_model = stan(model_file, data = poi_data, chains = 0)
# ,max_treedepth=14
stanfit = stan(fit = poi_model, data = poi_data,
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
if (LeadT<=3) {
  
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""))
}else{
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""),"lp__")
}
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = pars_p, inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = pars_p, las = 1)
dev.off()

#Mandleshwar

i=2
if (LeadT<=3) {
  
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""))
}else{
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""),"lp__")
}
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = pars_p, inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = pars_p, las = 1)
dev.off()
#Handia
i=3
if (LeadT<=3) {
  
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""))
}else{
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""),"lp__")
}
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = pars_p, inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = pars_p, las = 1)
dev.off()
#Hoshangabad
i=4
if (LeadT<=3) {
  
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""))
}else{
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""),"lp__")
}
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = pars_p, inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = pars_p, las = 1)
dev.off()
#Sandiya
i=5
if (LeadT<=3) {
  
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""))
}else{
  pars_p=c(paste(nam_par[1:(length(nam_par))],"[",i,"]",sep=""),"lp__")
}
jpeg(paste("traceplot_",info$name[i],".jpg",sep=""), width =7, height = 6, units = 'in', res = 300)
traceplot(stanfit, pars = pars_p, inc_warmup = TRUE, nrow = 3)
dev.off()
jpeg(paste("pairs_",info$name[i],".jpg",sep=""), width =6, height = 6, units = 'in', res = 300)
pairs(stanfit,  pars = pars_p, las = 1)
dev.off()  

##save file with posterior samples of parameters
aa=aa1
rm(stanfit,aa1)
list.save(aa, paste("Par_",pred,".rds",sep=""))



# generate posterior samples of streamflow --------------------------------
NN=nrow(aa$alpha_loc0)
post_y=array(NA,dim = c(dim(aa$alpha_loc0),nrow(data)))
dimnames(post_y)=list(1:NN,as.character(info$code),data[,1])
for (i in 1:5) {
  post_y[,i,]=aa$y_rep[,,i]
}

##compute some deterministic metrics
Val=matrix(NA,5,2)
colnames(Val)=c("R_BHM","%BIAS_BHM")
rownames(Val)=info$name
for (j in 1:5) {
  tmp=apply(post_y[,j,],2,mean)
  Val[j,1]=round(cor.test(tmp,data[,j+1],na.action("na.omit"),method="pearson")$estimate,3)
  Val[j,2]=round((mean(tmp,na.rm=T)-mean(data[,j+1],na.rm=T))/mean(data[,j+1],na.rm=T)*100,3)
  
}
print(Val)

write.table(Val, file = paste("DeterministicMetrics_",N_sea,".txt",sep=""),col.names = TRUE,row.names=FALSE,sep="\t")

# plot of posterior streamflow at each gauge with Copula ------------------------------

jpeg(paste("TS_",info$name[1],"_",info$name[5],".jpg",sep=""), width =9, height = 12, units = 'in', res = 300)
par(mfrow=c(5,1))
par(mar = c(2, 3.2, 1.4, 0.3))

i=1
zz=boxplot(split(t(post_y[,i,]),data[,1]),plot=F,cex=1.0,range=2)
zz$names=rep("",length(zz$names))
z1=bxp(zz,ylim=c(0,range(zz$stats,data[,i+1],na.rm = T)[2]),xlim=c(1.5,40.5),xlab="",ylab="",cex=1.3,outline=FALSE,notch = FALSE,whiskcol="gray",boxcol="black",boxfill="gray90",medcol="black",show.names=FALSE,axes=FALSE)
box()
n1=as.character(data[,1])
axis(1,at=z1[seq(1,length(z1),4)],labels=n1[seq(1,length(z1),4)],cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
axis(2,cex.axis=1.2,las=1, tck = -0.02,mgp=c(0, 0.5, 0))
mtext(bquote("N days with Q>"*Q[80]), side = 2, line = 1.8, cex = 0.9,font = 1)
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
mtext(bquote("N days with Q>"*Q[80]), side = 2, line = 1.8, cex = 0.9,font = 1)
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
mtext(bquote("N days with Q>"*Q[80]), side = 2, line = 1.8, cex = 0.9,font = 1)
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
mtext(bquote("N days with Q>"*Q[80]), side = 2, line = 1.8, cex = 0.9,font = 1)
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
mtext(bquote("N days with Q>"*Q[80]), side = 2, line = 1.8, cex = 0.9,font = 1)
lines(z1,data[,i+1],lty=1,lwd=1,col="#0b5394")
points(z1,data[,i+1],lwd=2,col="#0b5394")
title(info$name[i], cex.main = 1.7,line=0.3, adj = 0.5)


dev.off()

# Save posterior flow samples ---------------------------------------------

list.save(post_y, paste("post_Q_cali_",pred,".rds",sep=""))
sprintf("elapsed time: %s minutes",round((proc.time() - ptm1)[3]/60,2))