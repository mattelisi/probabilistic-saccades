# This script fir the asymmetric-cost model, as well as alternative models
# store the parameters, and make some plots.
# Matteo Lisi, 2019

rm(list=ls())
setwd("~/git_local/OSF_repo_probsaccades/")
library(mlisi)
d <- read.table("./data/saccades_xp123_allfit.txt",sep="\t",header=T)

### --------------------------------------------------------------------------------------- ###
# this code estimate the asymmetry by splitting the data according to session
source("asymCostFitFunctions_multi.R")
d$fit_level <- factor(paste(d$sess_n, d$posu,sep="_"))
d$id_i <- factor(paste(d$id, d$exp,sep="_"))
d$cond <- as.numeric(d$fit_level)
d$err_g <- d$gain - 1 # saccade error in gain units

# run fit on all subjects
startpar <- c(0.12, 0.12, 0.14, 0.14, 0.18, 0.18, 0.7) # initial guesses
startparR <- c(0.12, 0.12, 0.14, 0.14, 0.18, 0.18, 0.634483, 0.8) # robust model
par_LB_R <- c(1e-3, 1e-3, 1e-3,1e-3, 1e-3, 1e-3, 1e-4, 0)
par_UB_R <-  c(0.5, 0.5, 0.5,0.5, 0.5, 0.5, 10, 1)
dfit <- {}
dfit_2 <- {}
dtestfit <- {}

bootStdSE <- function (v, nsim = 1000, ...){
  bootFOO <- function(v, i) sd(v[i], na.rm = T, ...)
  bootRes <- boot::boot(v, bootFOO, nsim)
  return(sd(bootRes$t, na.rm = T))
}

cat("Begin...\n")
for(i in unique(d$id_i)){
  
  # fit assuming quadratic exponent
  fit0 <- optim(par=startpar, fn=negloglik_multi, d=d[d$id_i==i,], hessian=T, method="L-BFGS-B", lower=c(1e-3, 1e-3, 1e-3,1e-3, 1e-3, 1e-3, 0), upper=c(0.5, 0.5, 0.5,0.5, 0.5, 0.5, 1))
  fit0 <- optim(par=fit0$par, fn=negloglik_multi, d=d[d$id_i==i,], hessian=T)
  par <- fit0$par
  par_se <- sqrt(diag(MASS::ginv(fit0$hessian)))
  
  ideal_gain <- 1 + idealAimpoint_v(fit0$par[1:6],fit0$par[7])
  
  # observed value of gain and standard deviations
  observed_gain <- tapply(d$gain[d$id_i==i],d$fit_level[d$id_i==i], mean)
  observed_gain_se <- tapply(d$gain[d$id_i==i],d$fit_level[d$id_i==i], bootMeanSE)	
  observed_sd <- tapply(d$gain[d$id_i==i],d$fit_level[d$id_i==i], sd)		
  observed_sd_se <- tapply(d$gain[d$id_i==i],d$fit_level[d$id_i==i], bootStdSE)		
  
  # fit robust model
  fitR <- optim(par=startparR, fn=negloglik_multi_R, d=d[d$id_i==i,], hessian=T, method="L-BFGS-B", lower=par_LB_R, upper=par_UB_R)
  par_R <- fitR$par
  par_se_R <- sqrt(diag(MASS::ginv(fitR$hessian)))
  ideal_gain_R <- 1 + idealAimpoint_v_R(fitR$par[1:6],fitR$par[7:8])
  
  # collect info
  id <- rep(unique(d$id[d$id_i==i]),6) 
  if(unique(d$exp[d$id_i==i])=="exp1"){
    posu <- c(1,2,3,1,2,3)
    sigma <- c(0.3, 0.9, 1.5, 0.3, 0.9, 1.5)
    PL <- c(146, 57, 50, 146, 57, 50)
    session <- c(rep("large",3),rep("small",3))
    exp <- rep("exp1",6) 
    
  }else if(unique(d$exp[d$id_i==i])=="exp2"){
    posu <- c(1,2,3,1,2,3)
    sigma <- c(0.9, 0.9, 0.9, 0.3, 0.9, 1.5)
    PL <- c(146, 57, 50, 146, 146, 146)
    session <- c(rep("fixed-size",3),rep("fixed-luminance",3))
    exp <- rep("exp2",6) 
    
  }else if(unique(d$exp[d$id_i==i])=="exp3"){
    posu <- c(1,2,3,1,2,3) 
    sigma <- c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9)
    PL <- c(146, 57, 50, 146, 57, 50)
    session <- c(rep("large",3),rep("small",3))
    exp <- rep("exp3",6) 
  }
  
  # store results
  dfit <- rbind(dfit, data.frame(id, exp,session,PL,sigma,posu, STD=par[1:6], STD_se=par_se[1:6], alpha=par[7], alpha_se=par_se[7], id=i, ideal_gain, ideal_gain_R, observed_gain, observed_gain_se, observed_sd, STD_R=par_R[1:6], STD_se=par_se_R[1:6], Lshp=par_R[7],  Lshp_se=par_se_R[7], alpha_R=par_R[8], alpha_se_R=par_se_R[8],observed_sd_se))
  
  # store likelihood values
  n_i <- nrow(d[d$id_i==i,])
  L_a <- -fit0$value # store likelihood
  L_R <- -fitR$value # minus sign
  
  # compute information theoretic criteria
  # these are an approximation of the predictive quality of the model
  # so for the model comparison we will use the cross-validation instead
  AIC_a   <- 2*7 - 2*L_a
  AIC_R   <- 2*8 - 2*L_R
  
  msea <- MSE_a(fit0$par, d[d$id==i,])
  mser <- MSE_R(fitR$par, d[d$id==i,])
  
  AICc_a   <- AIC_a + (2*7^2 + 2*7)/(n_i - 7 - 1)
  AICc_R   <- AIC_R + (2*8^2 + 2*8)/(n_i - 8 - 1)
  
  BIC_a   <- log(n_i)*7 - 2*L_a
  BIC_R   <- log(n_i)*8 - 2*L_R
  
  dtestfit <- rbind(dtestfit, data.frame(id=unique(d$id[d$id_i==i]), exp=unique(d$exp[d$id_i==i]),L_a, msea, mser, L_R, L_R, AIC_a, AIC_R, AICc_a, AICc_R, BIC_a, BIC_R))
  
  cat("dataset ",i, "processed. moving on...\n")
  
}
# dfit$exp <- substr(dfit$id.1,4,7)

### --------------------------------------------------------------------------------------- ###
# save results
#write.table(dfit, file="./data/fit_all_multi_v2.txt")
#write.table(dtestfit, file="./data/test_fit_all_multi_v2.txt")
cat("...done!\n")

dfit <- read.table("./data/fit_all_multi_v2.txt")
dtestfit <- read.table("./data/test_fit_all_multi_v2.txt")

# 
dag_alpha <- aggregate(alpha~id+exp, dfit, mean)
median(dag_alpha$alpha[dag_alpha$exp=="exp1"])
median(dag_alpha$alpha[dag_alpha$exp=="exp2"])
median(dag_alpha$alpha[dag_alpha$exp=="exp3"])

library(mlisi)
bootFooCI(dag_alpha$alpha[dag_alpha$exp=="exp1"], nsim=10^4, foo="median")
bootFooCI(dag_alpha$alpha[dag_alpha$exp=="exp2"], nsim=10^4, foo="median")
bootFooCI(dag_alpha$alpha[dag_alpha$exp=="exp3"], nsim=10^4, foo="median")


### --------------------------------------------------------------------------------------- ###
# some summary plots

# graphical ggplot parameters
library(ggplot2) # nicer theme
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

### --------------------------------------------------------------------------------------- ###
# plot individual fits with predictions
nd <- expand.grid(STD=seq(0.02,0.4,length.out=100),id.1=unique(dfit$id.1),gain=NA, KEEP.OUT.ATTRS = F)
R_i <- vector("numeric",length(unique(dfit$id.1)))
err_i <- vector("numeric",length(unique(dfit$id.1)))
alpha_i <- vector("numeric",length(unique(dfit$id.1)))

i_c <- 0
for(i in unique(dfit$id.1)){
  i_c <- i_c + 1
  a_i <- mean(dfit$alpha[dfit$id.1==i])
  nd$gain[nd$id.1==i] <- 1+idealAimpoint_v(nd$STD[nd$id.1==i], a_i)
  R_i[i_c] <- summary(lm(ideal_gain~observed_gain,dfit[dfit$id.1==i,]))$r.squared
  err_i[i_c] <- with(dfit[dfit$id.1==i,], mean(abs(ideal_gain-observed_gain)))
  alpha_i[i_c] <- mean(dfit$alpha[dfit$id.1==i])
}
names(R_i)<- unique(dfit$id.1)
names(err_i)<- unique(dfit$id.1)

nd$observed_gain <- nd$gain
nd$observed_gain_se <- NA

dfit2 <- aggregate(alpha~id,dfit,mean)
dfit2$label <- paste("cost ratio =",round(dfit2$alpha/(1-dfit2$alpha),digits=0))

dfit$id_n <- paste(as.numeric(dfit$id.1)," (",substr(dfit$id.1,4,7),")",sep="")
nd$id_n <- paste(as.numeric(nd$id.1)," (",substr(nd$id.1,4,7),")",sep="")

### --------------------------------------------------------------------------------------- ###
# plot for paper

# shape code factor
dfit$shape <- ifelse(dfit$exp=="exp1","fixed-energy",ifelse(dfit$exp=="exp2",as.character(dfit$session),"fixed-size"))
dfit$shape <- ifelse(dfit$exp=="exp2",dfit$shape, paste(dfit$shape," (",dfit$session,")",sep=""))
dfit_avg <- aggregate(cbind(ideal_gain,observed_gain)~exp+posu,dfit,mean)

## 1
pdf("predicted_gain_split_R1.pdf",width=5,height=3)
ggplot(dfit, aes(x=ideal_gain,y=observed_gain,ymin=observed_gain-observed_gain_se,ymax=observed_gain+observed_gain_se,color=factor(posu),shape=shape))+geom_abline(intercept=0,slope=1,lty=2,size=0.4)+geom_point(size=1,alpha=0.8)+coord_equal(xlim=c(0.7,1.12),ylim=c(0.7,1.12))+scale_color_manual(values=c("black","dark grey","blue"),guide=F)+nice_theme+labs(x="predicted gain",y="observed gain")+facet_grid(exp~posu)+scale_shape_manual(values=c(21,19,17,15,22,15),name=NULL)+geom_hline(data=dfit_avg,aes(yintercept=observed_gain,color=factor(posu)),size=0.4,lty=3)+geom_vline(data=dfit_avg,aes(xintercept=ideal_gain,color=factor(posu)),size=0.4,lty=3)### --------------------------------------------------------------------------------------- ###

#+geom_errorbar(width=0,color="grey")
dev.off()

## 1.2
dfit_SD_avg <- aggregate(cbind(STD,observed_sd)~exp+posu,dfit,mean)
pdf("predicted_SD_split_R1.pdf",width=5,height=3)
ggplot(dfit, aes(x=STD,y=observed_sd,color=factor(posu),shape=shape))+geom_abline(intercept=0,slope=1,lty=2,size=0.4)+geom_point(size=1,alpha=0.8)+scale_color_manual(values=c("black","dark grey","blue"),guide=F)+nice_theme+labs(x="predicted gain SD",y="observed gain SD")+facet_grid(exp~posu)+scale_shape_manual(values=c(21,19,17,15,22,15),name=NULL)+geom_hline(data=dfit_SD_avg,aes(yintercept=observed_sd,color=factor(posu)),size=0.4,lty=3)+geom_vline(data=dfit_SD_avg,aes(xintercept=STD,color=factor(posu)),size=0.4,lty=3)+coord_equal(xlim=c(0.01,0.34),ylim=c(0.01,0.34))
#+geom_errorbar(width=0,color="grey")
dev.off()

## selected some example subjects
## and plot individual bias-variance relationship
selsjs <- c("10 (exp1)","1 (exp1)","25 (exp2)","15 (exp2)","39 (exp3)","21 (exp3)")
nd$id <- substr(nd$id.1,1,2)
nd$exp <- substr(nd$id.1,4,7)

dfit_i <- dfit[is.element(dfit$id_n,selsjs),]
nd_i <- nd[is.element(nd$id_n,selsjs),]

dfit_i$id_old <- dfit_i$id_n 
dfit_i$id_n <- factor(ifelse(is.element(dfit_i$id,c("aa","db","gh")),1,2))
nd_i$id_n <- factor(ifelse(is.element(nd_i$id,c("aa","db","gh")),1,2))
nd_i$shape <- NA

pdf("all_individual_fits_R1.pdf",width=3,height=2.4)
ggplot(dfit_i, aes(x=STD,y=observed_gain,color=factor(posu),shape=shape,group=id))+geom_line(data=nd_i,aes(x=STD,y=gain),col="black")+geom_errorbar(aes(ymin=observed_gain-1.96*observed_gain_se,ymax=observed_gain+1.96*observed_gain_se),width=0)+geom_point(size=1.5)+scale_color_manual(values=c("black","dark grey","blue"),guide=F)+nice_theme+labs(x="gain SD (model-based estimate)",y="observed gain")+scale_shape_manual(values=c(21,19,17,15,22,15),guide=F)+facet_grid(id_n~exp)+coord_cartesian(xlim=c(0.05,0.23),ylim=c(0.7,0.95))+scale_x_continuous(breaks=seq(0,0.3,0.1))
dev.off()

