
#包含线性交互模拟数据
library(brms)
library(rstan)
library(BhGLM)
library(survival)
library(survminer)
library(Hmisc)

set.seed(23456)
beta<-rnorm(11,0.7,0.05)
beta
deta<-c(-1,1,1,-1,-1,
        1,-1,1,
        -1,1,1)
beta<-beta*deta  
beta

nz=c(1,2,3,4,5,201,202,205,206,207,208)





testsimsurv<- function(n=1, size=500, nsnp=200, nclin=5, corr=0.5, 
                       sigma,nz=nz){
  simdatasum<-data.frame()
  for(j in 1:n){
    X=sim.x(n=size,m=nsnp,corr=corr)
    clin=sim.x(n=size,m=nclin,corr = corr)
    names(clin)<-c("x10001","x10002","x10003","x10004","x10005")
    var=cbind(X,clin)
    
    for(i in 1:nsnp){
      x10005x<-var[,i]*clin[,5]
      var<-cbind(var,x10005x)
    }
    aa=nsnp+nclin+1
    bb=nsnp+nclin+nsnp
    names(var)[aa:bb]<-as.list(paste("x10005x", c(1:nsnp), sep = ""))
    
    nn<-rep(j, times=size)
    eta = sim.eta(x=var[,nz],mu=0,coefs=beta) 
    etaa<-eta$eta
    y.normal <- rnorm(size, etaa, sigma)
    y <- rexp(size, exp(etaa))
    c <- rexp(size, exp(y.normal - etaa))
    time <- ifelse(y > c, c, y) # min(y, c)
    status <- ifelse(y > c, 0, 1) # 1: uncensored; 0: censored
    clidata <- data.frame(nn,time,status)
    
    total<-cbind(clidata,var)
    simdatasum<-rbind(simdatasum,total)
  }
  return(simdatasum)
}


set.seed(12345)
data<-testsimsurv(n=100, size=500, nsnp=200, nclin=5, 
                  corr=0.5, sigma=1,
                  nz=nz)

min=min(data$time)
max=max(data$time)
table(data$status)



set.seed(1234)
train<-testsimsurv(n=100, size=100, nsnp=200, nclin=5, 
                   corr=0.5, sigma=1,
                   nz=nz)









setwd("pathway")
#horseshoe
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  horseshoe=set_prior("horseshoe(df=1,df_global=1)")
  set.seed(20220902)
  f <- brm(time|cens(1-status) ~ .,
           data = datanew, 
           prior=horseshoe,
           chain=3,iter=4000, 
           #control = list(adapt_delta = 0.99,max_treedepth=12),
           family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))))
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  a<-result[["fixed"]]
  beta.hat<-a[2:406,1]
  beta.hatt<-as.vector(beta.hat)
  
  yrep_mean<-as.matrix(trainnew[,3:407]) %*% beta.hatt+a[1,1]
  C_index=rcorr.cens(yrep_mean,Surv(trainnew$time,trainnew$status))
  
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}












setwd("pathway")
#horseshoe
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  ###horseshoe
  nsnp=200
  nvar=205
  nz=205
  stanvars=stanvar(x=nsnp,name="nsnp",scode="int nsnp;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
    stanvar(x=nz,name="nz",scode="int nz;",block="data")
  horseshoe=set_prior("horseshoe(df=1,df_global=1)")
  
  sdata<-make_standata(time|cens(1-status) ~ .,
                       data = datanew,
                       #family = cox(),
                       family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
                       prior=horseshoe,stanvars=stanvars)
  set.seed(20220902)
  kf=rstan::stan(file ="pathway\\sshorseshoe_cox.stan",
                 data=sdata,
                 #control = list(adapt_delta = 0.99,max_treedepth=12),
                 chain=3,iter=4000)
  
  # feed the Stan model back into brms
  f <- brm(time|cens(1-status) ~ ., 
           data = datanew, 
           prior=horseshoe,stanvars=stanvars,
           family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
           empty = TRUE)
  f$fit <- kf
  f <- rename_pars(f)
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  a<-result[["fixed"]]
  beta.hat<-a[2:406,1]
  beta.hatt<-as.vector(beta.hat)
  
  yrep_mean<-as.matrix(trainnew[,3:407]) %*% beta.hatt+a[1,1]
  C_index=rcorr.cens(yrep_mean,Surv(trainnew$time,trainnew$status))
  
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}







setwd("pathway")
#ssde
#Multiple theta
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  nvar=405
  s0=0.05
  ssde=set_prior("for (j in 1:nvar)
               target += log_sum_exp(log(1-gamma[j])+normal_lpdf(b[j]|0,s0*tau),
                                     log(gamma[j])+normal_lpdf(b[j]|0,tau))",check=F)+
    set_prior("target += inv_gamma_lpdf(tau | 1, 10)",check=F)+
    set_prior("for (j in 1:nvar)
               target += beta_lpdf(gamma[j]|1,1)",check=F)
  
  stanvars=stanvar(scode="real<lower=0> tau;",block="parameters")+
    stanvar(x=s0,name="s0",scode="real s0;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
    stanvar(scode="vector<lower=0,upper=1>[nvar] gamma;",block="parameters")
  
  #set.seed(20220902)
  f=brm(time|cens(1-status) ~ .,
        data = datanew, 
        family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
        prior=ssde,stanvars=stanvars,
        chain=3,iter=5000)
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  a<-result[["fixed"]]
  beta.hat<-a[2:406,1]
  beta.hatt<-as.vector(beta.hat)
  
  yrep_mean<-as.matrix(trainnew[,3:407]) %*% beta.hatt+a[1,1]
  C_index=rcorr.cens(yrep_mean,Surv(trainnew$time,trainnew$status))
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}






setwd("pathway")
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  ###effect heredity2
  nsnp=200 
  nvar=205
  s0=0.05
  
  eefssde=set_prior("for (j in 1:nsnp)
               target += log_sum_exp(log(1-gamma[j])+normal_lpdf(b[j]|0,s0*tau),
                                     log(gamma[j])+normal_lpdf(b[j]|0,tau))",check=F)+
    set_prior("for (j in (nsnp+1):nvar)
               target += normal_lpdf(b[j]|0,tau)",check=F)+
    set_prior("for (j in (nvar+1):(nvar+nsnp))
               target += log_sum_exp(log(1-gamma[j-nvar])+normal_lpdf(b[j]|0,s0*tau),
                                     log(gamma[j-nvar])+normal_lpdf(b[j]|0,tau))",check=F)+
    set_prior("target += inv_gamma_lpdf(tau | 1, 10)",check=F)+
    set_prior("for (j in 1:nsnp)
               target += beta_lpdf(gamma[j]|1,1)",check=F)
  
  eefstanvars=stanvar(scode="vector<lower=0,upper=1>[nsnp] gamma;",block="parameters")+
    stanvar(scode="real<lower=0> tau;",block="parameters")+
    stanvar(x=s0,name="s0",scode="real s0;",block="data")+
    stanvar(x=nsnp,name="nsnp",scode="int nsnp;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")
  
  set.seed(20220902)
  f=brm(time|cens(1-status) ~ .,
        data = datanew, 
        family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
        prior=eefssde,stanvars=eefstanvars,
        chain=3,iter=5000)
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  a<-result[["fixed"]]
  beta.hat<-a[2:406,1]
  beta.hatt<-as.vector(beta.hat)
  
  yrep_mean<-as.matrix(trainnew[,3:407]) %*% beta.hatt+a[1,1]
  C_index=rcorr.cens(yrep_mean,Surv(trainnew$time,trainnew$status))
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}




