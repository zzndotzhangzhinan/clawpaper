
## ------------------------------------------------------------------------------------------
source("coordinate_func2.R")


## ----include=FALSE-------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(lemon)
library(doParallel)
library(foreach)
library(adaptMT)
library(splines)
library(IHW)


## ----warning=FALSE-------------------------------------------------------------------------
m = 3000
al = 0.05
f0 = c(0,1)
h = 150 #weights bandwidth

ppvec = rep(0.02, m)
ppvec[201:350]<-0.6
ppvec[1501:1650]<-0.6
ppvec[801:1000]<-0.3
ppvec[2101:2300]<-0.3
#muvec = seq(2.4,3,length.out=10)
rvec = seq(-0.9,0.9,length.out=10)

fdr.new = NA
tdp.new = NA
fdr.bh = NA
tdp.bh = NA
fdr.ad = NA
tdp.ad = NA
fdr.law = NA
tdp.law = NA
fdr.sab = NA
tdp.sab = NA
fdr.adapt = NA
tdp.adapt = NA

for (k in 1:length(rvec) ){
  mu = 3#muvec[k]
  r = rvec[k]
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.bh = NA
  subtdp.bh = NA
  subfdr.ad = NA
  subtdp.ad = NA
  subfdr.law = NA
  subtdp.law = NA
  subfdr.sab = NA
  subtdp.sab = NA
  subfdr.adapt = NA
  subtdp.adapt = NA
  
  for (l in 1:200){
    theta = rbinom(m, size=1, prob=ppvec)
    ts = arima.sim(list(order=c(1,0,0), ar=r), n=m, sd=sqrt(1-r^2))
    noi = rnorm(m,0,0.1)
    
    x = (ts+noi)/sqrt(1.01) + theta*mu
    #x = rep(0,m)
    #x[1:m/2] = (1-theta[1:m/2])*rnorm(m/2,f0[1],f0[2]) + theta[1:m/2]*rnorm(m/2,mu,1)
    #x[(1+m/2):m] = (1-theta[(1+m/2):m])*rnorm(m/2,f0[1],f0[2]) +
    #  theta[(1+m/2):m]*rnorm(m/2,0.8*mu,1.5)
    y = (ts+rnorm(m,0,0.1))/sqrt(1.01)
    
    #new mwthod
    pwscore = cor_1D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #BH
    pv = 2*(1-pnorm(abs(x),0,1))
    
    de.bh = bh(pv,al)
    subfdr.bh[l] = sum(de.bh*(1-theta)) / max(1, sum(de.bh) )
    subtdp.bh[l] = sum(de.bh*theta) / max(1, sum(theta) )
    
    #Adadetect
    score.pool = exscore.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    subfdr.ad[l] = sum(de.ad*(1-theta)) / max(1, sum(de.ad) )
    subtdp.ad[l] = sum(de.ad*theta) / max(1, sum(theta) )
    
    #LAWS
    pis = pis_1D.func(pv,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    subfdr.law[l] = sum(de.law*(1-theta)) / max(1, sum(de.law) )
    subtdp.law[l] = sum(de.law*theta) / max(1, sum(theta) )
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    subfdr.sab[l] = sum(de.sab*(1-theta)) / max(1, sum(de.sab) )
    subtdp.sab[l] = sum(de.sab*theta) / max(1, sum(theta) )
    
    ## AdaPT
    pvals = pv
    dist <- beta_family()
    formulas = paste("ns(x, df = ",c(2,4,6,8,10)," )")
    models <- lapply(formulas, function(formula){
    piargs <- muargs <- list(formula = formula)
    gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)})
    adapt.res <- adapt(x = data.frame(x = 1:m), pvals = pvals, models = models, 
                       alphas =c(0.05),dist = dist, nfits = 10)
    adapt.de = (adapt.res$qvals<=al)
    subfdr.adapt[l] = sum(adapt.de*(1-theta))/ max(sum(adapt.de),1)
    subtdp.adapt[l] = sum(adapt.de*theta)/ max(sum(theta),1)
    
    
    #cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.bh[k] = mean(subfdr.bh)
  tdp.bh[k] = mean(subtdp.bh)
  fdr.ad[k] = mean(subfdr.ad)
  tdp.ad[k] = mean(subtdp.ad)
  fdr.law[k] = mean(subfdr.law)
  tdp.law[k] = mean(subtdp.law)
  fdr.sab[k] = mean(subfdr.sab)
  tdp.sab[k] = mean(subtdp.sab)
  fdr.adapt[k] = mean(subfdr.adapt)
  tdp.adapt[k] = mean(subtdp.adapt)
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.1fdr = data.frame(rvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.1fdr = melt(ex1.1fdr,id="rvec")
colnames(ex1.1fdr) = c("rho","METHOD","value")


ex1.1tdp = data.frame(rvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.1tdp = melt(ex1.1tdp,id="rvec")
colnames(ex1.1tdp) = c("rho","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.1fdr = ggplot(ex1.1fdr)+aes(x=rho, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.1tdp = ggplot(ex1.1tdp)+aes(x=rho, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.1fdr)
print(line1.1tdp)


## ----warning=FALSE-------------------------------------------------------------------------
m = 3000
al = 0.05
f0 = c(0,1)
h = 150 #weights bandwidth

pivec = rep(0.3,10)#seq(0.15,0.4,length.out=10)

rvec = seq(-0.9,0.9,length.out=10)
#mu = 3

fdr.new = NA
tdp.new = NA
fdr.bh = NA
tdp.bh = NA
fdr.ad = NA
tdp.ad = NA
fdr.law = NA
tdp.law = NA
fdr.sab = NA
tdp.sab = NA
fdr.adapt = NA
tdp.adapt = NA

for (k in 1:length(rvec) ){
  ppvec = rep(0.02,m)
  ppvec[201:350]<-2*pivec[k]
  ppvec[1501:1650]<-2*pivec[k]
  ppvec[801:1000]<-pivec[k]
  ppvec[2101:2300]<-pivec[k]
  
  r = rvec[k]
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.bh = NA
  subtdp.bh = NA
  subfdr.ad = NA
  subtdp.ad = NA
  subfdr.law = NA
  subtdp.law = NA
  subfdr.sab = NA
  subtdp.sab = NA
  subfdr.adapt = NA
  subtdp.adapt = NA
  
  for (l in 1:200){
    theta = rbinom(m, size=1, prob=ppvec)
    
    ts = arima.sim(list(order=c(1,0,0), ar=r), n=m, sd=sqrt(1-r^2))
    noi = rnorm(m,0,0.1)
    
    x = rep(0,m)
    x[1:m/2] = (ts[1:m/2]+noi[1:m/2])/sqrt(1.01) + theta[1:m/2]*(-2.5)
    x[(1+m/2):m] = (1-theta[(1+m/2):m])*(ts[(1+m/2):m]+noi[(1+m/2):m])/sqrt(1.01) +
      theta[(1+m/2):m]*rnorm(m/2,3.6,1.5)
    
    y = (ts+rnorm(m,0,0.1))/sqrt(1.01)
    
    #new mwthod
    pwscore = cor_1D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #BH
    pv = 2*(1-pnorm(abs(x),0,1))
    
    de.bh = bh(pv,al)
    subfdr.bh[l] = sum(de.bh*(1-theta)) / max(1, sum(de.bh) )
    subtdp.bh[l] = sum(de.bh*theta) / max(1, sum(theta) )
    
    #Adadetect
    score.pool = exscore.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    subfdr.ad[l] = sum(de.ad*(1-theta)) / max(1, sum(de.ad) )
    subtdp.ad[l] = sum(de.ad*theta) / max(1, sum(theta) )
    
    #LAWS
    pis = pis_1D.func(pv,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    subfdr.law[l] = sum(de.law*(1-theta)) / max(1, sum(de.law) )
    subtdp.law[l] = sum(de.law*theta) / max(1, sum(theta) )
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    subfdr.sab[l] = sum(de.sab*(1-theta)) / max(1, sum(de.sab) )
    subtdp.sab[l] = sum(de.sab*theta) / max(1, sum(theta) )
    
    ## AdaPT
    pvals = pv
    dist <- beta_family()
    formulas = paste("ns(x, df = ",c(2,4,6,8,10)," )")
    models <- lapply(formulas, function(formula){
    piargs <- muargs <- list(formula = formula)
    gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)})
    adapt.res <- adapt(x = data.frame(x = 1:m), pvals = pvals, models = models, 
                       alphas =c(0.05),dist = dist, nfits = 10)
    adapt.de = (adapt.res$qvals<=al)
    subfdr.adapt[l] = sum(adapt.de*(1-theta))/ max(sum(adapt.de),1)
    subtdp.adapt[l] = sum(adapt.de*theta)/ max(sum(theta),1)
    
    #cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.bh[k] = mean(subfdr.bh)
  tdp.bh[k] = mean(subtdp.bh)
  fdr.ad[k] = mean(subfdr.ad)
  tdp.ad[k] = mean(subtdp.ad)
  fdr.law[k] = mean(subfdr.law)
  tdp.law[k] = mean(subtdp.law)
  fdr.sab[k] = mean(subfdr.sab)
  tdp.sab[k] = mean(subtdp.sab)
  fdr.adapt[k] = mean(subfdr.adapt)
  tdp.adapt[k] = mean(subtdp.adapt)
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.2fdr = data.frame(rvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.2fdr = melt(ex1.2fdr,id="rvec")
colnames(ex1.2fdr) = c("rho","METHOD","value")


ex1.2tdp = data.frame(rvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.2tdp = melt(ex1.2tdp,id="rvec")
colnames(ex1.2tdp) = c("rho","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.2fdr = ggplot(ex1.2fdr)+aes(x=rho, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.2tdp = ggplot(ex1.2tdp)+aes(x=rho, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.2fdr)
print(line1.2tdp)


## ----warning=FALSE-------------------------------------------------------------------------
m = 3000
al = 0.05
f0 = c(0,1)
h = 150 #weights bandwidth

indx = 1:m
ppvec = rep(0.02,m)
pdex = 0.4*(1+sin(0.02*indx))
ppvec[201:500]<-pdex[201:500]
ppvec[1501:1800]<-pdex[1501:1800]
ppvec[801:1100]<-pdex[801:1100]
ppvec[2101:2400]<-pdex[2101:2400]

#muvec = seq(2.4,3,length.out=10)
rvec = seq(-0.9,0.9,length.out=10)

fdr.new = NA
tdp.new = NA
fdr.bh = NA
tdp.bh = NA
fdr.ad = NA
tdp.ad = NA
fdr.law = NA
tdp.law = NA
fdr.sab = NA
tdp.sab = NA
fdr.adapt = NA
tdp.adapt = NA

for (k in 1:length(rvec) ){
  mu = 3#muvec[k]
  mumu = mu+0.15*sin(0.6*indx) 
  
  r = rvec[k]
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.bh = NA
  subtdp.bh = NA
  subfdr.ad = NA
  subtdp.ad = NA
  subfdr.law = NA
  subtdp.law = NA
  subfdr.sab = NA
  subtdp.sab = NA
  subfdr.adapt = NA
  subtdp.adapt = NA
  
  for (l in 1:200){
    theta = rbinom(m, size=1, prob=ppvec)
    ts = arima.sim(list(order=c(1,0,0), ar=r), n=m, sd=sqrt(1-r^2))
    noi = rnorm(m,0,0.1)
    
    x = (ts+noi)/sqrt(1.01) + theta*mumu
    #x = rep(0,m)
    #x[1:m/2] = (1-theta[1:m/2])*rnorm(m/2,f0[1],f0[2]) + theta[1:m/2]*rnorm(m/2,mu,1)
    #x[(1+m/2):m] = (1-theta[(1+m/2):m])*rnorm(m/2,f0[1],f0[2]) +
    #  theta[(1+m/2):m]*rnorm(m/2,0.8*mu,1.5)
    y = (ts+rnorm(m,0,0.1))/sqrt(1.01)

    
    #new mwthod
    pwscore = cor_1D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #BH
    pv = 2*(1-pnorm(abs(x),0,1))
    
    de.bh = bh(pv,al)
    subfdr.bh[l] = sum(de.bh*(1-theta)) / max(1, sum(de.bh) )
    subtdp.bh[l] = sum(de.bh*theta) / max(1, sum(theta) )
    
    #Adadetect
    score.pool = exscore.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    subfdr.ad[l] = sum(de.ad*(1-theta)) / max(1, sum(de.ad) )
    subtdp.ad[l] = sum(de.ad*theta) / max(1, sum(theta) )
    
    #LAWS
    pis = pis_1D.func(pv,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    subfdr.law[l] = sum(de.law*(1-theta)) / max(1, sum(de.law) )
    subtdp.law[l] = sum(de.law*theta) / max(1, sum(theta) )
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    subfdr.sab[l] = sum(de.sab*(1-theta)) / max(1, sum(de.sab) )
    subtdp.sab[l] = sum(de.sab*theta) / max(1, sum(theta) )
    
    ## AdaPT
    pvals = pv
    dist <- beta_family()
    formulas = paste("ns(x, df = ",c(2,4,6,8,10)," )")
    models <- lapply(formulas, function(formula){
    piargs <- muargs <- list(formula = formula)
    gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)})
    adapt.res <- adapt(x = data.frame(x = 1:m), pvals = pvals, models = models, 
                       alphas =c(0.05),dist = dist, nfits = 10)
    adapt.de = (adapt.res$qvals<=al)
    subfdr.adapt[l] = sum(adapt.de*(1-theta))/ max(sum(adapt.de),1)
    subtdp.adapt[l] = sum(adapt.de*theta)/ max(sum(theta),1)
    
    #cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.bh[k] = mean(subfdr.bh)
  tdp.bh[k] = mean(subtdp.bh)
  fdr.ad[k] = mean(subfdr.ad)
  tdp.ad[k] = mean(subtdp.ad)
  fdr.law[k] = mean(subfdr.law)
  tdp.law[k] = mean(subtdp.law)
  fdr.sab[k] = mean(subfdr.sab)
  tdp.sab[k] = mean(subtdp.sab)
  fdr.adapt[k] = mean(subfdr.adapt)
  tdp.adapt[k] = mean(subtdp.adapt)
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.3fdr = data.frame(rvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.3fdr = melt(ex1.3fdr,id="rvec")
colnames(ex1.3fdr) = c("rho","METHOD","value")


ex1.3tdp = data.frame(rvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.3tdp = melt(ex1.3tdp,id="rvec")
colnames(ex1.3tdp) = c("rho","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.3fdr = ggplot(ex1.3fdr)+aes(x=rho, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+geom_line(size=0.7)+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.3tdp = ggplot(ex1.3tdp)+aes(x=rho, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.3fdr)
print(line1.3tdp)

