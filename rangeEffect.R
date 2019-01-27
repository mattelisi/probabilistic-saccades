# This script document the analysis of the central bias in Exp. 1 and 3
# This require that Stan and RStan are installed
# The lines that run the MCMC sampling are commented (see below)
# will need to be run at least once to reproduce the plots in the article
# Matteo Lisi, last update: 27/1/2019
#

# clear workspace, set working directory, load data table
rm(list=ls())
setwd("~/git_local/OSF_repo_probsaccades/")
library(mlisi)
d <- read.table("./data/saccades_xp123_allfit.txt",sep="\t",header=T)

library(mlisi)
library(ggplot2) # nicer theme
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

# select only experiments where the range was manipulated
d <- d[d$exp!="exp2",]

# apply corrections for differences across testing day in Exp. 1 data (see Supporting Information for details)
d_exp1 <- d[d$exp=="exp1",]
d_exp1$sacXresp_c <- d_exp1$sacXresp
d_exp1$sigmaf <- factor(round(d_exp1$sigma,digits=1))
barG <- aggregate(gain ~ id, d_exp1[d_exp1$sigmaf=="0.3",], mean)
Gs1 <- aggregate(gain ~ id, d_exp1[d_exp1$sigmaf=="0.3" & d_exp1$sess_n==1,], mean)
Gs2 <- aggregate(gain ~ id, d_exp1[d_exp1$sigmaf=="0.3" & d_exp1$sess_n==2,], mean)
for(i in 1:nrow(d_exp1)){
  if(d_exp1$sess_n[i]==1){
    d_exp1$sacXresp_c[i] <- d_exp1$sacXresp[i] * (1+barG$gain[barG$id==d_exp1$id[i]]-Gs1$gain[Gs1$id==d_exp1$id[i]])
  }else{
    d_exp1$sacXresp_c[i] <- d_exp1$sacXresp[i] * (1+barG$gain[barG$id==d_exp1$id[i]]-Gs2$gain[Gs2$id==d_exp1$id[i]])
  }
}

# substitue corrected amplitude values for experiment 1
d$sacXresp[d$exp=="exp1"] <- d_exp1$sacXresp_c

# dataset in stan-friendly format
id_<- data.frame(labels=unique(d_exp1$id),id_n=1:length(unique(d_exp1$id)))
d_exp1$id_n <- NA
for(i in 1:nrow(d_exp1)){
  d_exp1$id_n[i] <- id_$id_n[id_$labels==d_exp1$id[i]]
}
dxp1 <- list(N=nrow(d_exp1),
             J=length(unique(d_exp1$id)),
             posu=d_exp1$posu,
             E=d_exp1$tarX,
             S=d_exp1$sacXresp_c,
             id=d_exp1$id_n,
             meanE=ifelse(d_exp1$session=="small",9,11))

# run sampling - exp1
library(rstan)
options(mc.cores = parallel::detectCores()) # tell stan to use multiple cores
# m_stan <- stan(file = "range_model.stan", data = dxp1, iter = 2000, chains = 4) # This run MCMC sampling
# saveRDS(m_stan,file="rangemodel_exp1.rds")
m_stan <- readRDS("rangemodel_exp1.rds") # load results if already run

# some plots to evaluate convergence
traceplot(m_stan, pars = c("beta"), inc_warmup = F)
traceplot(m_stan, pars = c("alpha"), inc_warmup = F)
print(m_stan, pars = c("alpha"), probs = c(0.025, 0.975),digits=3)

# data from experiment 3
d_exp3 <- d[d$exp=="exp3",]
id_<- data.frame(labels=unique(d_exp3$id),id_n=1:length(unique(d_exp3$id)))
d_exp3$id_n <- NA
for(i in 1:nrow(d_exp3)){
  d_exp3$id_n[i] <- id_$id_n[id_$labels==d_exp3$id[i]]
}
dxp3 <- list(N=nrow(d_exp3),
             J=length(unique(d_exp3$id)),
             posu=d_exp3$posu,
             E=d_exp3$tarX,
             S=d_exp3$sacXresp,
             id=d_exp3$id_n,
             meanE=ifelse(d_exp3$session=="small",8,11))

# run sampling - exp3
#m_stan3 <- stan(file = "range_model.stan", data = dxp3, iter = 2000, chains = 4) # This run MCMC sampling
#saveRDS(m_stan3,file="rangemodel_exp3.rds")
m_stan3 <- readRDS("rangemodel_exp3.rds") 

# convergence
traceplot(m_stan3, pars = c("beta"), inc_warmup = F)
traceplot(m_stan3, pars = c("alpha"), inc_warmup = F)

# confidence intervals for the group-level alpha parameters
print(m_stan3, pars = c("alpha"), probs = c(0.025, 0.975),digits=3)

# create an new dataframe to compute the models predictions
nd1 <- expand.grid(E = seq(8,12,length.out=100), posu=1:3, meanE=c(9,11),exp="exp1")
nd1 <- nd1[(nd1$E<=10 & nd1$meanE==9) | (nd1$E>=10 & nd1$meanE==11),]
nd3 <- expand.grid(E = seq(6,13,length.out=100), posu=1:3, meanE=c(8,11),exp="exp3")
nd3 <- nd3[(nd3$E<=10 & nd3$meanE==8) | (nd3$E>=9 & nd3$meanE==11),]

# reproduce here the fixed-effects structure of the model
# to compute the mean predictions
predict_range_model <- function(m,nd){
  beta <- extract(m, pars = "beta")$beta 
  alpha <- extract(m, pars = "alpha")$alpha 
  nd$mu<-rep(NA,nrow(nd))
  nd$lb <- nd$mu; nd$ub <- nd$mu
  for(i in 1:nrow(nd)){
    y_ <- beta[,nd$posu[i]] + beta[,3+nd$posu[i]]*(alpha[,nd$posu[i]] * nd$meanE[i] + (1-alpha[,nd$posu[i]])*nd$E[i])
    nd$mu[i] <- mean(y_)
    nd$lb[i] <- quantile(y_, probs = 0.05/2)
    nd$ub[i] <- quantile(y_, probs = 1-0.05/2)
  }
  nd$mu_gain <- nd$mu / nd$E
  nd$lb_gain <- nd$lb / nd$E
  nd$ub_gain <- nd$ub / nd$E
  return(nd)
}

# predict
nd3 <- predict_range_model(m_stan3, nd3)
nd1 <- predict_range_model(m_stan, nd1)
nd <- rbind(nd3, nd1)
nd$session <- ifelse(nd$meanE<10,"small","large")

# means
d$group_f <- paste(d$tarX, d$exp, d$posu, d$session, sep="_")
dag <- aggregate(gain~group_f + tarX + exp + posu + session, d, mean)

# standard errors (pooled, taking imbalances in number of trials into account)
bootPoolMean <- function (M, N, nsim = 1000) {
  d <- data.frame(M, N)
  poolFoo <- function(d){
    M_pool <- sum(d$M * (d$N)) / sum(d$N)
    return(M_pool)
  }
  bootFoo <- function(d, i) poolFoo(d[i,])
  bootRes <- boot::boot(d, bootFoo, nsim)
  return(sd(bootRes$t, na.rm = T))
}
dag_V <- aggregate(gain~group_f + tarX + exp + posu + session +id, d, var)
dag_V$n <- aggregate(gain~group_f + tarX + exp + posu + session +id, d, length)$gain
dag_M <- aggregate(gain~group_f + tarX + exp + posu + session +id, d, mean)
dag_M$n <- aggregate(gain~group_f + tarX + exp + posu + session +id, d, length)$gain
dag$gain_2 <- NA
dag$se <- NA
for(i in 1:nrow(dag)){
  M <- dag_M$gain[dag_M$group_f==dag$group_f[i] & dag_M$exp==dag$exp[i] & dag_M$posu==dag$posu[i]]
  N <- dag_M$n[dag_M$group_f==dag$group_f[i] & dag_M$exp==dag$exp[i] & dag_M$posu==dag$posu[i]]
  dag$gain_2[i] <- sum(M * (N)) / sum(N)
  dag$se[i] <- bootPoolMean(M,N)
}

# plot with model predictions
dag$mu_gain <- dag$gain_2
dag$lb_gain <- dag$gain_2 - dag$se
dag$ub_gain <- dag$gain_2 + dag$se
dag$E <- dag$tarX

pdf("range_withPredictions.pdf",width=4,height=2.5)
ggplot(dag, aes(x=E, y=mu_gain, ymin=lb_gain, ymax=ub_gain, group=session, color=session,fill=session))+geom_ribbon(data=nd,alpha=0.4,color=NA)+geom_line(data=nd,size=0.3)+geom_errorbar(width=0,size=0.4)+geom_point(pch=21,size=0.8,fill=NA)+facet_grid(exp~posu)+nice_theme+scale_color_manual(values=c("black","dark grey"))+scale_fill_manual(values=c("black","dark grey"))+labs(x="target distance [deg]", y="saccadic gain")
dev.off()

# now plot alpha parameters (proportion of compression)
print(m_stan, pars = c("alpha"), probs = c(0.025, 0.975),digits=2)
print(m_stan3, pars = c("alpha"), probs = c(0.025, 0.975),digits=2)
alpha_1 <- extract(m_stan, pars = "alpha")$alpha 
alpha_3 <- extract(m_stan3, pars = "alpha")$alpha 
x1 <- sort(unique(d_exp1$sigma))
x3 <- sort(unique(d_exp3$PL),decreasing=T)
da1 <- data.frame(x=x1, alpha=apply(alpha_1,2,mean), se=apply(alpha_1,2,sd),x_f=factor(round(x1,digits=1)))
da3 <- data.frame(x=x3, alpha=apply(alpha_3,2,mean), se=apply(alpha_3,2,sd),x_f=factor(round(x3)))

pdf(file="compression_exp1.pdf",height=1.3,width=1.2)
ggplot(da1, aes(x=x_f,y=alpha,ymin=alpha-se,ymax=alpha+se,group=1,color=x_f))+geom_hline(yintercept=0,lty=2,size=0.4)+geom_line(color="black")+geom_errorbar(width=0,alpha=1,color="black")+geom_point(size=2.3)+nice_theme+labs(x=expression(paste(sigma," [deg]")), y=expression(paste("central bias [",alpha,"]")))+scale_color_manual(values=c("black","dark grey","blue"),guide=F)
dev.off()

da3$x_f <- factor(da3$x_f, levels=c("146","57","50"))
pdf(file="compression_exp3.pdf",height=1.3,width=1.2)
ggplot(da3, aes(x=x_f,y=alpha,ymin=alpha-se,ymax=alpha+se,group=1,color=x_f))+geom_hline(yintercept=0,lty=2,size=0.4)+geom_line(color="black")+geom_errorbar(width=0,alpha=1,color="black")+geom_point(size=2.3)+nice_theme+labs(x=expression(paste("lum. [cd/",m^2,"]")), y=expression(paste("central bias [",alpha,"]")))+scale_color_manual(values=c("black","dark grey","blue"),guide=F)
dev.off()
