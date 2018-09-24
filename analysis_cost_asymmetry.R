### --------------------------------------------------------------------------------------- ###
# analysis of asymmetrical costs in saccade targeting

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
  dfit <- rbind(dfit, data.frame(id, exp,session,PL,sigma,posu, STD=par[1:6], STD_se=par_se[1:6], alpha=par[7], alpha_se=par_se[7], id=i, ideal_gain, ideal_gain_R, observed_gain, observed_gain_se, observed_sd, STD_R=par_R[1:6], STD_se=par_se_R[1:6], Lshp=par_R[7],  Lshp_se=par_se_R[7], alpha_R=par_R[8], alpha_se_R=par_se_R[8]))
  
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
write.table(dfit, file=".data/fit_all_multi.txt")
write.table(dtestfit, file=".data/test_fit_all_multi.txt")
cat("...done!\n")

### --------------------------------------------------------------------------------------- ###
### --------------------------------------------------------------------------------------- ###
# some summary plots

library(ggplot2) # nicer theme
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

### --------------------------------------------------------------------------------------- ###
# predicted and observed saccadic gain
pdf("predicted_gain_split.pdf",width=5,height=1.8)
ggplot(dfit, aes(x=ideal_gain,y=observed_gain,ymin=observed_gain-observed_gain_se,ymax=observed_gain+observed_gain_se,color=factor(posu),shape=session))+geom_abline(intercept=0,slope=1,lty=2,size=0.4)+geom_errorbar(width=0)+geom_point(size=1)+coord_equal(xlim=c(0.7,1.12),ylim=c(0.7,1.12))+scale_color_manual(values=c("black","dark grey","blue"),guide=F)+nice_theme+labs(x="predicted gain",y="observed gain")+facet_grid(.~exp)+scale_shape_manual(values=c(21,19,21,19))
dev.off()

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

# 
dfit$id_n <- factor(dfit$id_n, unique(dfit$id_n)[order(alpha_i)],ordered=T)
nd$id_n <- factor(nd$id_n, unique(nd$id_n)[order(alpha_i)],ordered=T)

pdf("all_individual_fits.pdf",width=6,height=5)
ggplot(dfit, aes(x=STD,y=observed_gain,color=factor(posu),shape=session))+geom_line(data=nd,aes(x=STD,y=gain,shape="large"),col="black")+geom_errorbar(aes(ymin=observed_gain-1.96*observed_gain_se,ymax=observed_gain+1.96*observed_gain_se),width=0)+geom_point(size=1)+scale_color_manual(values=c("black","dark grey","blue"),guide=F)+nice_theme+labs(x="gain variability [SD]",y="gain")+coord_cartesian(xlim=c(0.05,0.33),ylim=c(0.7,1.13))+scale_shape_manual(values=c(21,19,21,19),guide=F)+facet_wrap(~id_n, ncol=10)#+scale_x_log10()
dev.off()

### --------------------------------------------------------------------------------------- ###
# secondary saccades (old cost analysis)

d2 <- read.table("./data/secondary_saccades_xp123_allfit.txt",sep="\t",header=T)

# explo plot
ggplot(d2, aes(x=rt,y=errorChange, color=factor(sacN_1)))+geom_hline(yintercept=0,size=0.2,lty=2)+nice_theme+labs(x="latency", y="change in position error")+geom_point(pch=21,size=1)+facet_grid(.~posu)

# cut offs
d2$corrective <- ifelse(d2$errorChange< 2.5 & d2$errorChange>-2.5, 1, 0)
d2 <- d2[d2$sacN_1==2,]

ggplot(d2, aes(x=rt,y=errorChange, color=corrective))+geom_hline(yintercept=0,size=0.2,lty=2)+nice_theme+labs(x="latency", y="change in position error")+geom_point(pch=21,size=1)+facet_grid(.~posu)

d2 <- d2[d2$corrective==1,]
d2 <- d2[d2$rt>=30,]

# plot overall latencies distributions
dag <- aggregate(rt ~ posu + dir, d2, mean)
ggplot(d2, aes(x=rt, fill=factor(posu)))+geom_histogram(binwidth=25,alpha=0.85)+geom_vline(data=dag, aes(xintercept=rt, color=factor(posu)))+nice_theme+labs(x="secondary saccade latency [ms]")+scale_fill_manual(values=c("black","dark grey", "blue"),name=expression(paste(sigma[blob]," [deg]")))+scale_color_manual(values=c("black","dark grey", "blue"),name=expression(paste(sigma[blob]," [deg]")))+facet_grid(dir~posu)


# given that there are relatively few observations, and that latency differences are generally noisy, I use the median which is more robust and less influenced by single extreme points (call them "outliers" if you want)
medianDiff <- function(dd) median(dd$res_RT[dd$dir=="opposite to primary"])-median(dd$res_RT[dd$dir=="same as primary"])
bootmedianDiffSE <- function(dd,nsim=1000){
  bootF <- function(dd,i) medianDiff(dd[i,])
  bootRes <- boot::boot(dd,bootF,nsim)
  return(sd(bootRes$t,na.rm=T))
}

# fit model at group level, including quadratic component
library(lme4)
library(MASS)

m2all <- lmer(rt ~ ampli + I(ampli^2) + (ampli+I(ampli^2)|id) + exp, d2)

summary(m2all)
d2$res_RT <- residuals(m2all)

# pair values (alpha & sacc. pars.)
d2$id <- paste(d2$id, d2$exp,sep="_")
d_a <- aggregate(alpha ~ id + exp, dfit, mean)
d_a$alpha_se <- aggregate(alpha_se ~ id + exp, dfit, mean)$alpha_se
d_a$id <- paste(d_a$id, d_a$exp,sep="_")

d_a$RTdiff <- NA
d_a$RTdiff_se <- NA
d_a$errdiff <- NA
d_a$errChangediff <- NA

for(i in unique(d_a$id)){
  
  dd <- d2[d2$id==i,]
  
  d_a$RTdiff[d_a$id==i] <- medianDiff(dd)
  d_a$RTdiff_se[d_a$id==i] <- bootmedianDiffSE(dd,nsim=1000)
}

d_a$cost_ratio <- d_a$alpha / (1-d_a$alpha)
d_a$log_cost_ratio <- log(d_a$cost_ratio)

# measure correlation
mean(d_a$RTdiff)
bootMeanCI(d_a$RTdiff)
sum(d_a$RTdiff_se>=60)
with(d_a[d_a$RTdiff_se<60,], cor.test(log_cost_ratio, RTdiff))

pdf("cost_corr_1.pdf",width=5,height=1.8)
ggplot(d_a[d_a$RTdiff_se<60,], aes(x=log_cost_ratio,y=RTdiff,ymax=RTdiff+RTdiff_se,ymin=RTdiff-RTdiff_se))+geom_hline(yintercept=0,lty=3,col="grey",size=0.4)+stat_ellipse(type="norm",level=0.75,col="grey",size=0.3)+stat_ellipse(type="norm",level=0.95,col="grey",size=0.3)+geom_errorbar(width=0,color="dark grey")+geom_point(size=1.5,pch=19)+nice_theme+labs(x=expression(paste("log"," ", bgroup("(",frac(italic("cost overshoot"),italic("cost undershoot")),")"))), y="secondary saccades\nlatency difference [ms]")+facet_grid(.~exp)+coord_cartesian(ylim=c(-150,200))
dev.off()

pdf("cost_corr_2.pdf",width=2,height=2)
ggplot(d_a[d_a$RTdiff_se<60,], aes(x=log_cost_ratio,y=RTdiff,ymax=RTdiff+RTdiff_se,ymin=RTdiff-RTdiff_se,color=exp))+geom_hline(yintercept=0,lty=2,col="grey",size=0.5)+stat_ellipse(type="norm",level=0.75,col="grey",size=0.3)+stat_ellipse(type="norm",level=0.95,col="grey",size=0.3)+geom_errorbar(width=0,color="dark grey")+geom_point(size=2,pch=19)+geom_point(size=2.05,pch=21,color="black")+nice_theme+labs(x=expression(paste("log"," ", bgroup("(",frac(italic("cost overshoot"),italic("cost undershoot")),")"))), y="secondary saccades\nlatency difference [ms]")+scale_color_viridis_d(begin=0,end=0.7)+scale_y_continuous(breaks=seq(-200,200,50))+scale_x_continuous(breaks=seq(-3,6,1))
dev.off()
