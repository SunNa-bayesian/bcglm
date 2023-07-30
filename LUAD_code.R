
library(brms)
library(rstan)
library(BhGLM)
library(survival)
library(survminer)
library(Hmisc)

#set.seed(12345)
horseshoe=set_prior("horseshoe(df=1,df_global=1)")
hf <- brm(Time|cens(1-OS) ~ .,
          data = data,
          prior=horseshoe,
          chain=3,iter=5000, 
          #control = list(adapt_delta = 0.99,max_treedepth=12),
          family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))))
result<-summary(hf)
resultt<-result[["fixed"]]
hsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

hkfold <- kfold(hf, chains = 3, K=5,iter=5000,
                save_fits = TRUE)
hreskfold<-hkfold[["estimates"]]





#set.seed(12345)
nsnp=56 #Adjustments based on data
nvar=63 #Adjustments based on data
nz=62   #Adjustments based on data
stanvars=stanvar(x=nsnp,name="nsnp",scode="int nsnp;",block="data")+
  stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
  stanvar(x=nz,name="nz",scode="int nz;",block="data")
horseshoe=set_prior("horseshoe(df=1,df_global=1)")

sdata<-make_standata(Time|cens(1-OS) ~ .,
                     data = data, 
                     family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
                     prior=horseshoe,stanvars=stanvars)

kf=rstan::stan(file ="pathway\\sshorseshoe_cox.stan",
               data=sdata,
               #control = list(adapt_delta = 0.999,max_treedepth=12),
               #control = list(adapt_delta = 0.99),
               chain=3,iter=4000)

# feed the Stan model back into brms
fit <- brm(Time|cens(1-OS) ~ .,
           data = data, 
           prior=horseshoe,stanvars=stanvars,
           family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
           empty = TRUE)
fit$fit <- kf
fit <- rename_pars(fit)

result<-summary(fit)
resultt<-result[["fixed"]]
ehsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

ehkfold <- kfold(fit, chains = 3, K=5,iter=4000,
                 save_fits = TRUE)
ehreskfold<-ehkfold[["estimates"]]







nvar=119
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
hpf=brm(Time|cens(1-OS) ~ .,
        data = data, 
        family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
        prior=ssde,stanvars=stanvars,
        chain=3,iter=5000)

result<-summary(hpf)
resultt<-result[["fixed"]]
hpsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

hpkfold <- kfold(hpf, chains = 3, K=5,iter=5000,
                 save_fits = TRUE)
hpreskfold<-hpkfold[["estimates"]]







nsnp=56
nvar=63
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

#set.seed(20220902)
hphf=brm(Time|cens(1-OS) ~ .,
         data = data, 
         family = brmsfamily("cox", bhaz = list(Boundary.knots = c(min, max))),
         prior=eefssde,stanvars=eefstanvars,
         chain=3,iter=5000)

result<-summary(hphf)
resultt<-result[["fixed"]]
hphsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

hphkfold <- kfold(hphf, chains = 3, K=5,iter=5000,
                  save_fits = TRUE)
hphreskfold<-hphkfold[["estimates"]]





