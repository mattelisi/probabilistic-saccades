### quadratic loss function fitting
# this version use 3 * 2 conditions for fitting

### quadratic loss // closed-form for expected loss
loss <- function(err, a){
  (a + (1-2*a)*(err<0)) * abs(err)^2
}

expected_loss <- function(aim_p, sig, a){
  (aim_p^2+sig^2)/2 + ((1-2*a)*(aim_p^2+sig^2)*erf(-aim_p/(sqrt(2)*sig)))/2 + ((2*a-1)*aim_p*sig *exp(-(aim_p^2)/(2*sig^2)))/ sqrt(2*pi)
}

# derivative of expected loss with respect to the aim point
ELderiv <- function(aim_p, sig, a){
  u <- aim_p; s<-sig;
  e_1 <- (2*a-1)*u*erf(u/(sqrt(2)*s)) +u +2*s*(2*a-1)*exp(-(u^2)/(2*s^2))/(sqrt(2*pi))
  return(e_1)
}

# compute ideal aimpoint using gradient descent # unnecessary (Brent seems actually faster)
# idealAimpoint <- function(sig, a){
# 	stepSize<-1
# 	convTh<-1e-5
# 	maxIt<-1e8
# 	converged <- F
# 	it <- 0
# 	g <- -0.1
# 	while(converged == F & it<maxIt){
# 		it <- it+1 
# 		g_new <- g - stepSize*ELderiv(g, sig, a )
# 		if(abs(g_new - g)<=convTh){
# 			converged<-T
# 		}
# 		g <- g_new
# 	}
# 	if(converged == F){
# 		warning(paste(it,">",maxIt,"; max number of iterations exceeded!", sep=""))
# 	}
# 	return(g_new)
# }

idealAimpoint <- function(sig, a){
  exp_cost <- function(aim_p) expected_loss(aim_p, sig, a)
  return(optim(-0.05, exp_cost, method = "Brent", lower = -2, upper = 2)$par)
}

MSE_a <- function(par, d){
  n <- length(par)-1
  MSE <- 0
  for(i in 1:n){
    MSE <- MSE + sum(d$err_g[d$cond==i]-idealAimpoint(par[i], par[n+1]))^2
  }
  return(MSE/nrow(d))
}

idealAimpoint_v <- function(sig, a){
  aimp <- vector("numeric",length(sig))
  for(i in 1:length(sig)){
    aimp[i] <- idealAimpoint(sig[i],a)
  }
  return(aimp)
}

# negative log likelihood function for a single subject
# par = c(sig1, sig2, ... , alpha)
negloglik_multi<- function(par, d){
  n <- length(par)-1
  aimpoint <- vector("numeric",n)
  L <- 0
  for(i in 1:n){
    aimpoint[i] <- idealAimpoint(par[i], par[n+1])
    L <- L - sum(dnorm(d$err_g[d$cond==i], mean=aimpoint[i], sd= par[i], log=TRUE))
  }
  return(L)
}

### ------------------------------------------------------------------------------------- ###
# log likelihood for parametric (null) models

# two parameters (mean+sd) for each condition
# par = c(sd1, sd2, ...,  intercept, slope)
negloglik_m00 <-function(d, par){
  n <- length(par)-2
  L <- 0
  for(i in 1:n){
    L <- L - sum(dnorm(d$err_g[d$cond==i], mean=par[n+1]+par[i]*par[n+2], sd= par[i], log=TRUE))
    
  }
  return(L)
}

MSE_00 <- function(par, d){
  n <- length(par)-2
  MSE <- 0
  for(i in 1:n){
    MSE <- MSE + sum(d$err_g[d$cond==i]-(par[n+1]+par[i]*par[n+2]))^2
  }
  return(MSE/nrow(d))
}


### ------------------------------------------------------------------------------------- ###
# alternative version with robust loss function

# Robust loss function (downweight extreme errors)
# typical values for k, which provide high efficiency also in the normal case, are then
# k = sigma_hat * 4.685 # for Bisquare function
# (these values are taken from John Fox notes on robust regression,
# avalable here: http://users.stat.umn.edu/~sandy/courses/8053/handouts/robust.pdf)
bisquare_loss <- function(err, k){
  # this start quadratic then flattens to a constant value
  ifelse(abs(err)<=k, ((k^2)/6) * (1 - (1-(err/k)^2)^3),  (k^2)/6)
}

loss_R <- function(err, par){
  (par[2] + (1-2*par[2])*(err<0)) * bisquare_loss(err, k=par[1])
}

idealAimpoint_R <- function(sig, par){
  integrand <-function(x, aim_p) loss_R(x, par)*dnorm(x, mean=aim_p, sd= sig)
  #exp_cost <- function(aim_p) integrate(integrand, lower=-Inf, upper=Inf, aim_p=aim_p)$value
  exp_cost <- function(aim_p) integrate(integrand, lower=-Inf, upper=Inf, aim_p=aim_p, rel.tol = 1e-5)$value
  return(optim(0, exp_cost, method = "Brent", lower = -2, upper = 2)$par)
}

idealAimpoint_v_R <- function(sig, par){
  aimp <- vector("numeric",length(sig))
  for(i in 1:length(sig)){
    aimp[i] <- idealAimpoint_R(sig[i],par)
  }
  return(aimp)
}

# negative log likelihood function for a single subject
# par = c(sig1, sig2, ... , k, alpha)
negloglik_multi_R <- function(par, d){
  n <- length(par)-2
  aimpoint <- vector("numeric",n)
  L <- 0
  for(i in 1:n){
    integrand <-function(x, aim_p) loss_R(x, par[c(n+1,n+2)])*dnorm(x, mean=aim_p, sd= par[i])
    #exp_cost <- function(aim_p) integrate(integrand, lower=-Inf, upper=Inf, aim_p=aim_p)$value
    exp_cost <- function(aim_p) integrate(integrand, lower=-Inf, upper=Inf, aim_p=aim_p, rel.tol = 1e-5)$value
    aimpoint[i] <- optim(0, exp_cost, method = "Brent", lower = -2, upper = 2)$par
    L <- L - sum(dnorm(d$err_g[d$cond==i], mean=aimpoint[i], sd= par[i], log=TRUE))
  }
  return(L)
}

MSE_R <- function(par, d){
  n <- length(par)-2
  MSE <- 0
  for(i in 1:n){
    MSE <- MSE + sum(d$err_g[d$cond==i]-idealAimpoint_R(par[i], par[c(n+1,n+2)]))^2
  }
  return(MSE/nrow(d))
}

### ------------------------------------------------------------------------------------- ###
# perform cross-validation of the quadratic-asymmetric model
# and of the linear null-model above, byt holding one condition out at a time
hcoCrossValidation <- function(d){
  
  # split based on condition
  nCs <- length(unique(d$cond))
  Cs <- unique(d$cond)
  dout <- {}
  
  startpar <- c(0.13, 0.13, 0.13, 0.13, 0.13, 0.8)
  par_LB <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0)
  par_UB <-  c(0.5, 0.5, 0.5, 0.5, 0.5, 1)
  # startpar <- c(0.12, 0.12, 0.14, 0.14, 0.18, 0.18, 0.7) # initial guesses
  # startparR <- c(0.12, 0.12, 0.14, 0.14, 0.18, 0.18, 0.634483, 0.8) # robust model
  # par_LB_R <- c(1e-3, 1e-3, 1e-3,1e-3, 1e-3, 1e-3, 1e-4, 0)
  # par_UB_R <-  c(0.5, 0.5, 0.5,0.5, 0.5, 0.5, 10, 1)
  
  # robust model
  startparR <- c(0.13, 0.13, 0.13, 0.13, 0.13, 0.634483, 0.8)
  par_LB_R <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-4, 0)
  par_UB_R <-  c(0.5, 0.5, 0.5, 0.5, 0.5, 10, 1)
  
  # start par for null models
  startpar00 <- c(0.13, 0.13, 0.13, 0.13, 0.13, 0, -0.5)
  LB00 <- c(1e-2, 1e-2, 1e-2, 1e-2, 1e-2, -100, -100)
  UB00 <- c(0.5, 0.5, 0.5, 0.5, 0.5, 100, 100)
  
  cat(paste("condition-split cross-validation begin"))
  
  # begin
  t<-0
  for(i in 1:nCs){
    t <- t+1
    c_test <- Cs[i]
    
    dtest <- d[which(d$cond==c_test), ]
    dtrain <- d[which(d$cond!=c_test), ]
    
    # fit quadratic model
    fit1 <- optim(par=startpar, fn=negloglik_multi, d=dtrain, method="L-BFGS-B", lower=par_LB, upper=par_UB)
    #fit1 <- optim(par=fit1$par, fn=negloglik_multi, d=dtrain)
    L1 <- -negloglik_multi(fit1$par, dtest) # note the minus sign
    MSE1 <- MSE_a(fit1$par, dtest)
    
    # fit robust model
    fitR <- optim(par=startparR, fn=negloglik_multi_R, d=dtrain, method="L-BFGS-B", lower=par_LB_R, upper=par_UB_R)
    #fitR <- optim(par=fitR$par, fn=negloglik_multi_R, d=dtrain)
    L_R <- -negloglik_multi_R(fitR$par, dtest) # note the minus sign
    MSE_R <- MSE_R(fitR$par, dtest)
    
    # null models
    fit00 <- optim(par=startpar00, fn=negloglik_m00, d=dtrain, method="L-BFGS-B", lower=LB00, upper=UB00)
    #fit00 <- optim(par=fit00$par, fn=negloglik_m00, d=dtrain)
    L00 <- -negloglik_m00(dtest, fit00$par)
    MSE00 <- MSE_00(fit00$par, dtest)
    
    # store data
    dout <- rbind(dout, data.frame(t, c_test, sig1=fit1$par[1], sig2=fit1$par[2], sig3=fit1$par[3], alpha=fit1$par[4], L1, L00, L_R, MSE1, MSE00, MSE_R, sig1R=fitR$par[1], sig2R=fitR$par[2], sig3R=fitR$par[3], alphaR=fitR$par[5], lhspR=fitR$par[4]))
    
  }
  cat(" ... completed!\n")
  return(dout)
}