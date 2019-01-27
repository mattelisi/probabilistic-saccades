rm(list=ls())
setwd("~/git_local/OSF_repo_probsaccades/")

library(mlisi)
library(lme4)
library(ggplot2) # nicer theme

nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))


# new secondary saccade analysis
d <- read.table("./data/secondary_saccades_xp123_allfit.txt", header=T, sep="\t")
d$corrective <- ifelse(d$errorChange< 2.5 & d$errorChange>-5, 1, 0)
d <- d[d$sacN_1==2,] # only the first secondary saccade after the primary
d <- d[d$corrective==1,]
d <- d[d$rt>=30,] # latency > 30 ms
tapply(d$ampli, d$exp, mean)



# plot as cumulative distributions
d$dir2 <- ifelse(d$dir=="opposite to primary","backward","forward")
pdf("secondary_RT_ECDF.pdf",width=4.8,height=1.8)
ggplot(d, aes(x=rt,group=dir2,color=factor(posu),linetype=dir2))+geom_hline(yintercept=0.5,size=0.2,color="light grey",lty=1)+geom_vline(xintercept=seq(100,1000,100),size=0.2,color="light grey",lty=1)+ stat_ecdf(geom = "step",size=0.6)+facet_grid(.~posu)+nice_theme+scale_color_manual(values=c("black","dark grey","blue"),guide=F)+labs(y="cumulative probability",x="secondary saccade latency [ms]")+scale_linetype_manual(values=c(2,1),name="direction")+scale_x_continuous(breaks=seq(100,1000,100))+coord_cartesian(xlim=c(50,450))
dev.off()

# fit parameters
dfit <- read.table("./data/fit_all_multi.txt",header=T)

# bn data for plotting
d$X_bin <- cut(d$sacXresp, c(-10,-1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75 , 1, 1.25, 1.5, 10))
dag <- aggregate(cbind(sacXresp,ampli,rt) ~ X_bin, d, mean)
dag_i <- aggregate(rt ~ X_bin + id, d, mean)
dag_i$n <- aggregate(rt ~ X_bin + id, d, length)$rt
bootPoolMN <- function (M, N, nsim = 1000) {
  d <- data.frame(M, N)
  poolFoo <- function(d){
    M_pool <- sum(d$M * (d$N)) / sum(d$N)
    return(M_pool)
  }
  bootFoo <- function(d, i) poolFoo(d[i,])
  bootRes <- boot::boot(d, bootFoo, nsim)
  return(sd(bootRes$t, na.rm = T)) 
}
dag$rt_2 <- NA
dag$se <- NA
for(i in 1:nrow(dag)){
  M <- dag_i$rt[dag_i$X_bin==dag$X_bin[i]]
  N <- dag_i$n[dag_i$X_bin==dag$X_bin[i]]
  dag$rt_2[i] <- sum(M * N) / sum(N)
  dag$se[i] <- bootPoolMN(M,N)
}
dag$rt_median_lb <- dag$rt_2 - dag$se
dag$rt_median_ub <- dag$rt_2 + dag$se

# multilevel quadratic model
d$rt_sec <- d$rt/1000
d$id2 <- paste(d$id, d$exp,sep="_")
m2all <- lmer(rt_sec ~ ampli + I(ampli^2) + I(ampli^3) + (ampli + I(ampli^2)| id2), d)
summary(m2all)

#
d_plot <- expand.grid(ampli=seq(0,3,length.out=500),id2=NA) 
d_plot$rt <- predict(m2all, newdata=d_plot, re.form=NA)

#
mfun <- function(.) predict(., newdata=d_plot, re.form=NA)
#boot_m2all <- bootMer(m2all, mfun, nsim=100,.progress="txt")
#saveRDS(boot_m2all, file="boot_s2_model.rds")
boot_m2all <- readRDS("boot_s2_model.rds")
d_plot$rt_median_lb <- d_plot$rt - apply(boot_m2all$t,2,sd)
d_plot$rt_median_ub <- d_plot$rt + apply(boot_m2all$t,2,sd)
d_plot$side <- NA

dag$side <- ifelse(dag$sacXresp>0,"forward","backward")

dag$rt_median_lb <- dag$rt_median_lb/1000
dag$rt_median_ub <- dag$rt_median_ub/1000
dag$rt_2 <- dag$rt_2/1000
d_plot <- d_plot[d_plot$ampli<=max(dag$ampli),]

pdf("secondary_latency_binned_v2.pdf",width=1.6,height=2.4)
ggplot(dag, aes(x=ampli, y=rt_2, ymin=rt_median_lb,ymax=rt_median_ub,group=side, color=side))+geom_ribbon(data=d_plot, aes(x=ampli,y=rt),fill="black",size=0,alpha=0.1)+geom_line(data=d_plot, aes(x=ampli,y=rt),size=0.3,color="black",lty=2)+geom_line(size=0.3)+geom_errorbar(width=0,size=0.4)+geom_point(pch=19,size=1)+nice_theme+labs(x="secondary saccade\namplitude [deg]",y="mean latency [sec]")+theme(legend.justification = c(0, 0), legend.position = c(0, 0.6), legend.box.margin=margin(rep(0, times=4)),legend.background = element_blank(), legend.spacing = unit(0.1, 'cm'))+scale_fill_manual(guide=F)+scale_color_manual(values=c("red","dark green"),name="direction")+coord_cartesian(ylim=c(0.14,0.4),xlim=c(0,3))+scale_y_continuous(breaks=seq(0,1,0.05))
dev.off()


## check and save residuals
d$res_RT <- residuals(m2all)

## compute correlations
# pair values (alpha & sacc. pars.)
d$id2 <- paste(d$id, d$exp,sep="_")
d_a <- aggregate(alpha ~ id + exp, dfit, mean)
d_a$alpha_se <- aggregate(alpha_se ~ id + exp, dfit, mean)$alpha_se
d_a$id2 <- paste(d_a$id, d_a$exp,sep="_")

d_a$RTdiff <- NA
d_a$RTdiff_se <- NA
d_a$errdiff <- NA
d_a$errChangediff <- NA

meanDiff <- function(dd) mean(dd$res_RT[dd$dir=="opposite to primary"])-mean(dd$res_RT[dd$dir=="same as primary"])
bootmeanDiffSE <- function(dd,nsim=1000){
  bootF <- function(dd,i) meanDiff(dd[i,])
  bootRes <- boot::boot(dd,bootF,nsim)
  return(sd(bootRes$t,na.rm=T))
}

for(i in unique(d_a$id2)){
  dd <- d[d$id2==i,]
  d_a$RTdiff[d_a$id2==i] <- meanDiff(dd)
  d_a$RTdiff_se[d_a$id2==i] <- bootmeanDiffSE(dd,nsim=1000)
}
d_a$cost_ratio <- d_a$alpha / (1-d_a$alpha)
d_a$log_cost_ratio <- log(d_a$cost_ratio)

# measure correlation
round(mean(d_a$RTdiff)*1000)
round(bootMeanCI(d_a$RTdiff)*1000)

# cut-off on the precision of the estimated latency cost
rt_se_cut_off <- 0.03
sum(d_a$RTdiff_se>=rt_se_cut_off)
mean(d_a$RTdiff_se[d_a$RTdiff_se>=rt_se_cut_off])
mean(d_a$RTdiff_se[d_a$RTdiff_se<rt_se_cut_off])

with(d_a, cor.test(log_cost_ratio, RTdiff)) # correlation on all participants

with(d_a[d_a$RTdiff_se<rt_se_cut_off,], cor.test(log_cost_ratio, RTdiff)) # correlation after cut-off

# separate for experiments
with(d_a[d_a$exp=="exp1" & d_a$RTdiff_se<rt_se_cut_off,], cor.test(log_cost_ratio, RTdiff))
with(d_a[d_a$exp=="exp2" & d_a$RTdiff_se<rt_se_cut_off,], cor.test(log_cost_ratio, RTdiff))
with(d_a[d_a$exp=="exp3" & d_a$RTdiff_se<rt_se_cut_off,], cor.test(log_cost_ratio, RTdiff))

d_a <- d_a[seq(nrow(d_a),1,-1),] # reorder for plotting layers
pdf("cost_corr_2.2.pdf",width=3.1,height=2.5)
ggplot(d_a[d_a$RTdiff_se<rt_se_cut_off,], aes(x=log_cost_ratio,y=RTdiff,ymax=RTdiff+RTdiff_se,ymin=RTdiff-RTdiff_se,color=exp))+geom_hline(yintercept=0,lty=3,col="grey",size=0.5)+geom_vline(xintercept=0,lty=3,col="grey",size=0.5)+stat_ellipse(type="norm",level=0.75,col="black",size=0.3)+stat_ellipse(type="norm",level=0.95,col="black",size=0.3)+geom_errorbar(width=0)+geom_point(size=2,pch=19)+nice_theme+labs(x=expression(paste("log"," ", bgroup("(",frac(italic("cost overshoot"),italic("cost undershoot")),")"))), y="secondary saccades\nlatency difference [sec]")+scale_color_manual(values=c("black","red","dark grey"),name="Exp.")+scale_y_continuous(breaks=seq(-200,200,50)/1000)+scale_x_continuous(breaks=seq(-3,6,1))+coord_cartesian(ylim=c(-160,220)/1000,xlim=c(-2.8,5.8))
dev.off()

