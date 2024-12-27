
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
muvec = seq(2.4,3,length.out=10)

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

for (k in 1:length(muvec) ){
  mu = muvec[k]
  
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
    x = (1-theta)*rnorm(m,f0[1],f0[2]) + theta*rnorm(m,mu,1)
    #x = rep(0,m)
    #x[1:m/2] = (1-theta[1:m/2])*rnorm(m/2,f0[1],f0[2]) + theta[1:m/2]*rnorm(m/2,mu,1)
    #x[(1+m/2):m] = (1-theta[(1+m/2):m])*rnorm(m/2,f0[1],f0[2]) +
    #  theta[(1+m/2):m]*rnorm(m/2,0.8*mu,1.5)
    y = rnorm(m,f0[1],f0[2])
    
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
ex1.1fdr = data.frame(muvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.1fdr = melt(ex1.1fdr,id="muvec")
colnames(ex1.1fdr) = c("mu","METHOD","value")


ex1.1tdp = data.frame(muvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.1tdp = melt(ex1.1tdp,id="muvec")
colnames(ex1.1tdp) = c("mu","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.1fdr = ggplot(ex1.1fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.1tdp = ggplot(ex1.1tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
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

pivec = seq(0.15,0.4,length.out=10)

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

for (k in 1:length(pivec) ){
  ppvec = rep(0.02,m)
  ppvec[201:350]<-2*pivec[k]
  ppvec[1501:1650]<-2*pivec[k]
  ppvec[801:1000]<-pivec[k]
  ppvec[2101:2300]<-pivec[k]
  
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
    x = rep(0,m)
    x[1:m/2] = (1-theta[1:m/2])*rnorm(m/2,f0[1],f0[2]) + theta[1:m/2]*rnorm(m/2,-2.5,1)
    x[(1+m/2):m] = (1-theta[(1+m/2):m])*rnorm(m/2,f0[1],f0[2]) +
      theta[(1+m/2):m]*rnorm(m/2,3.6,1.5)
    y = rnorm(m,f0[1],f0[2])
    
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
ex1.2fdr = data.frame(pivec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.2fdr = melt(ex1.2fdr,id="pivec")
colnames(ex1.2fdr) = c("pi","METHOD","value")


ex1.2tdp = data.frame(pivec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.2tdp = melt(ex1.2tdp,id="pivec")
colnames(ex1.2tdp) = c("pi","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.2fdr = ggplot(ex1.2fdr)+aes(x=pi, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.2tdp = ggplot(ex1.2tdp)+aes(x=pi, y=value, group=METHOD, color=METHOD, shape=METHOD)+
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

muvec = seq(2.4,3,length.out=10)

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

for (k in 1:length(muvec) ){
  mu = muvec[k]
  mumu = mu+0.15*sin(0.6*indx) 
  
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
    #x = rep(0,m)
    #x[1:m/2] = (1-theta[1:m/2])*rnorm(m/2,f0[1],f0[2]) + 
    #  theta[1:m/2]*rnorm(m/2,mumu[1:m/2],1)
    #x[(1+m/2):m] = (1-theta[(1+m/2):m])*rnorm(m/2,f0[1],f0[2]) +
    #  theta[(1+m/2):m]*rnorm(m/2,-mumu[1:m/2],1.5)
    x = (1-theta)*rnorm(m,f0[1],f0[2]) + theta*rnorm(m,mumu,1)
    y = rnorm(m,f0[1],f0[2])
    
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
ex1.3fdr = data.frame(muvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.3fdr = melt(ex1.3fdr,id="muvec")
colnames(ex1.3fdr) = c("mu","METHOD","value")


ex1.3tdp = data.frame(muvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.3tdp = melt(ex1.3tdp,id="muvec")
colnames(ex1.3tdp) = c("mu","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.3fdr = ggplot(ex1.3fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.3tdp = ggplot(ex1.3tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.3fdr)
print(line1.3tdp)


## ------------------------------------------------------------------------------------------
grid_arrange_shared_legend(line1.1fdr,line1.2fdr,line1.3fdr,
                           line1.1tdp,line1.2tdp,line1.3tdp,nrow=2,ncol=3,position = "top")


## ------------------------------------------------------------------------------------------
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
pi0 = 1

# A and B: 2D region with signals
# A is a circle
# B is a square 

isA.func = function(i, j)
{
	# decide whether (i, j) is in circle A
	if((sqrt((i-30)^2+(j-70)^2)<=20) && (sqrt((i-30)^2+(j-70)^2)>=10))
	y=1
	else
	y=0
	return(y)
}

isB.func = function(i, j)
{
	# decide whether (i, j) is in the square B
	if(i>=62 && i<=90 && j>=10 && j<=38)
	y=1
	else
	y=0
	return(y)
}

for(i in 1:m1)
{
	for (j in 1:m1)
	{
		if (isA.func(i,j) || isB.func(i, j))
		thetaS[i,j] = rbinom(1, 1, pi0)
		else
		thetaS[i,j] = rbinom(1, 1, 0)
	}
}
#image(thetaS, main="True States, mu=2.0")
thetaS.vec = c(thetaS)


## ----warning=FALSE-------------------------------------------------------------------------
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

muvec = seq(2.4,3,length.out=10)

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

for (k in 1:length(muvec) ){
  mu = muvec[k]
  
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
  for (l in 1:200){
    thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 0.75)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0.02)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,mu,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-thetaS.vec)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    subfdr.bh[l] = sum(de.bh*(1-thetaS.vec)) / max(1, sum(de.bh) )
    subtdp.bh[l] = sum(de.bh*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    subfdr.ad[l] = sum(de.ad*(1-thetaS.vec)) / max(1, sum(de.ad) )
    subtdp.ad[l] = sum(de.ad*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    subfdr.law[l] = sum(de.law*(1-thetaS.vec)) / max(1, sum(de.law) )
    subtdp.law[l] = sum(de.law*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    subfdr.sab[l] = sum(de.sab*(1-thetaS.vec)) / max(1, sum(de.sab) )
    subtdp.sab[l] = sum(de.sab*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    
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
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.4fdr = data.frame(muvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab)
ex1.4fdr = melt(ex1.4fdr,id="muvec")
colnames(ex1.4fdr) = c("mu","METHOD","value")


ex1.4tdp = data.frame(muvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab)
ex1.4tdp = melt(ex1.4tdp,id="muvec")
colnames(ex1.4tdp) = c("mu","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.4fdr = ggplot(ex1.4fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.4tdp = ggplot(ex1.4tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.4fdr)
print(line1.4tdp)


## ----warning=FALSE-------------------------------------------------------------------------
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

pivec = seq(0.2,0.9,length.out=10)

mu = 2.8

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

for (k in 1:length(pivec) ){
  
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
  for (l in 1:200){
    thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, pivec[k])
		    else
		    thetaS[i,j] = rbinom(1, 1, 0.02)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,mu,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-thetaS.vec)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    subfdr.bh[l] = sum(de.bh*(1-thetaS.vec)) / max(1, sum(de.bh) )
    subtdp.bh[l] = sum(de.bh*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    subfdr.ad[l] = sum(de.ad*(1-thetaS.vec)) / max(1, sum(de.ad) )
    subtdp.ad[l] = sum(de.ad*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    subfdr.law[l] = sum(de.law*(1-thetaS.vec)) / max(1, sum(de.law) )
    subtdp.law[l] = sum(de.law*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    subfdr.sab[l] = sum(de.sab*(1-thetaS.vec)) / max(1, sum(de.sab) )
    subtdp.sab[l] = sum(de.sab*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    
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
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.5fdr = data.frame(pivec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab)
ex1.5fdr = melt(ex1.5fdr,id="pivec")
colnames(ex1.5fdr) = c("pi","METHOD","value")


ex1.5tdp = data.frame(pivec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab)
ex1.5tdp = melt(ex1.5tdp,id="pivec")
colnames(ex1.5tdp) = c("pi","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.5fdr = ggplot(ex1.5fdr)+aes(x=pi, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.5tdp = ggplot(ex1.5tdp)+aes(x=pi, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.5fdr)
print(line1.5tdp)


## ----warning=FALSE-------------------------------------------------------------------------
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth
mu = 2.5

Rvec = seq(13,25,length.out=10)

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

for (k in 1:length(Rvec) ){
  isA.func2 = function(i, j)
    {
	# decide whether (i, j) is in circle A
	 if((sqrt((i-30)^2+(j-70)^2)<=Rvec[k]) && (sqrt((i-30)^2+(j-70)^2)>=0.5*Rvec[k]))
	 y=1
	 else
	 y=0
	 return(y)
  }
  
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
  for (l in 1:200){
    thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func2(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 0.75)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0.02)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,mu,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-thetaS.vec)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    subfdr.bh[l] = sum(de.bh*(1-thetaS.vec)) / max(1, sum(de.bh) )
    subtdp.bh[l] = sum(de.bh*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    subfdr.ad[l] = sum(de.ad*(1-thetaS.vec)) / max(1, sum(de.ad) )
    subtdp.ad[l] = sum(de.ad*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    subfdr.law[l] = sum(de.law*(1-thetaS.vec)) / max(1, sum(de.law) )
    subtdp.law[l] = sum(de.law*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    subfdr.sab[l] = sum(de.sab*(1-thetaS.vec)) / max(1, sum(de.sab) )
    subtdp.sab[l] = sum(de.sab*thetaS.vec) / max(1, sum(thetaS.vec) )
    
    
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
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.6fdr = data.frame(Rvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab)
ex1.6fdr = melt(ex1.6fdr,id="Rvec")
colnames(ex1.6fdr) = c("R","METHOD","value")


ex1.6tdp = data.frame(Rvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab)
ex1.6tdp = melt(ex1.6tdp,id="Rvec")
colnames(ex1.6tdp) = c("R","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.6fdr = ggplot(ex1.6fdr)+aes(x=R, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.6tdp = ggplot(ex1.6tdp)+aes(x=R, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.6fdr)
print(line1.6tdp)


## ------------------------------------------------------------------------------------------
grid_arrange_shared_legend(line1.4fdr,line1.5fdr,line1.6fdr,
                           line1.4tdp,line1.5tdp,line1.6tdp,nrow=2,ncol=3,position = "top")


## ------------------------------------------------------------------------------------------
set.seed(123456)
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 1)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,1.5,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    thetaS.hat.new1 = matrix(de.new, dims[1], dims[2])
    
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    thetaS.hat.bh1 = matrix(de.bh, dims[1], dims[2])
    
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    thetaS.hat.ad1 = matrix(de.ad, dims[1], dims[2])
    
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    thetaS.hat.law1 = matrix(de.law, dims[1], dims[2])
    
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    thetaS.hat.sab1 = matrix(de.sab, dims[1], dims[2])
    


## ------------------------------------------------------------------------------------------
par(mfrow=c(2, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

image(thetaS, main="True States, mu=1.5")
image(thetaS.hat.new1, main="CLAW")
image(thetaS.hat.bh1, main="BH")
image(thetaS.hat.ad1, main="AdaDetect")
image(thetaS.hat.law1, main="LAWS")
image(thetaS.hat.sab1, main="SABHA")


## ------------------------------------------------------------------------------------------
set.seed(123456)
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 1)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,1.8,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    thetaS.hat.new2 = matrix(de.new, dims[1], dims[2])
    
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    thetaS.hat.bh2 = matrix(de.bh, dims[1], dims[2])
    
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    thetaS.hat.ad2 = matrix(de.ad, dims[1], dims[2])
    
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    thetaS.hat.law2 = matrix(de.law, dims[1], dims[2])
    
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    thetaS.hat.sab2 = matrix(de.sab, dims[1], dims[2])
    


## ------------------------------------------------------------------------------------------
par(mfrow=c(2, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

image(thetaS, main="True States, mu=1.8")
image(thetaS.hat.new2, main="CLAW")
image(thetaS.hat.bh2, main="BH")
image(thetaS.hat.ad2, main="AdaDetect")
image(thetaS.hat.law2, main="LAWS")
image(thetaS.hat.sab2, main="SABHA")


## ------------------------------------------------------------------------------------------
set.seed(123456)
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 1)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,2,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    thetaS.hat.new3 = matrix(de.new, dims[1], dims[2])
    
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    thetaS.hat.bh3 = matrix(de.bh, dims[1], dims[2])
    
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    thetaS.hat.ad3 = matrix(de.ad, dims[1], dims[2])
    
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    thetaS.hat.law3 = matrix(de.law, dims[1], dims[2])
    
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    thetaS.hat.sab3 = matrix(de.sab, dims[1], dims[2])
    


## ------------------------------------------------------------------------------------------
par(mfrow=c(2, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

image(thetaS, main="True States, mu=2")
image(thetaS.hat.new3, main="CLAW")
image(thetaS.hat.bh3, main="BH")
image(thetaS.hat.ad3, main="AdaDetect")
image(thetaS.hat.law3, main="LAWS")
image(thetaS.hat.sab3, main="SABHA")


## ------------------------------------------------------------------------------------------
set.seed(123456)
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 1)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,2.2,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    thetaS.hat.new4 = matrix(de.new, dims[1], dims[2])
    
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    thetaS.hat.bh4 = matrix(de.bh, dims[1], dims[2])
    
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    thetaS.hat.ad4 = matrix(de.ad, dims[1], dims[2])
    
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    thetaS.hat.law4 = matrix(de.law, dims[1], dims[2])
    
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    thetaS.hat.sab4 = matrix(de.sab, dims[1], dims[2])
    


## ------------------------------------------------------------------------------------------
par(mfrow=c(2, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

image(thetaS, main="True States, mu=2.2")
image(thetaS.hat.new4, main="CLAW")
image(thetaS.hat.bh4, main="BH")
image(thetaS.hat.ad4, main="AdaDetect")
image(thetaS.hat.law4, main="LAWS")
image(thetaS.hat.sab4, main="SABHA")


## ------------------------------------------------------------------------------------------
set.seed(123456)
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 1)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,2.5,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    thetaS.hat.new5 = matrix(de.new, dims[1], dims[2])
    
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    thetaS.hat.bh5 = matrix(de.bh, dims[1], dims[2])
    
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    thetaS.hat.ad5 = matrix(de.ad, dims[1], dims[2])
    
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    thetaS.hat.law5 = matrix(de.law, dims[1], dims[2])
    
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    thetaS.hat.sab5 = matrix(de.sab, dims[1], dims[2])
    


## ------------------------------------------------------------------------------------------
par(mfrow=c(2, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

image(thetaS, main="True States, mu=2.5")
image(thetaS.hat.new5, main="CLAW")
image(thetaS.hat.bh5, main="BH")
image(thetaS.hat.ad5, main="AdaDetect")
image(thetaS.hat.law5, main="LAWS")
image(thetaS.hat.sab5, main="SABHA")


## ------------------------------------------------------------------------------------------
set.seed(123456)
dims = c(100, 100)
m1 = 100
m = dims[1]*dims[2]
al = 0.05
f0 = c(0,1)
h = 15 #weights bandwidth

thetaS = matrix(rep(0, m), dims[1], dims[2]) # true states of nature
    for(i in 1:m1)
    {
	    for (j in 1:m1)
	    {
		    if (isA.func(i,j) || isB.func(i, j))
		    thetaS[i,j] = rbinom(1, 1, 1)
		    else
		    thetaS[i,j] = rbinom(1, 1, 0)
	    }
    }
    thetaS.vec = c(thetaS)
    
    x.vec = (1-thetaS.vec)*rnorm(m,f0[1],f0[2]) + thetaS.vec*rnorm(m,3,1)
    x = matrix(x.vec, dims[1], dims[2])
    
    y.vec = rnorm(m,f0[1],f0[2])
    y = matrix(y.vec, dims[1], dims[2])
    
    #new mwthod
    pwscore = cor_2D.func(x, y, tau=0.5, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    thetaS.hat.new6 = matrix(de.new, dims[1], dims[2])
    
    
    #BH
    pv = 2*(1-pnorm(abs(x.vec),0,1))
    
    de.bh = bh(pv,al)
    thetaS.hat.bh6 = matrix(de.bh, dims[1], dims[2])
    
    
    #Adadetect
    score.pool = exscore.func(x.vec,y.vec,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.ad = Adadetect(pools,poolst,al)
    thetaS.hat.ad6 = matrix(de.ad, dims[1], dims[2])
    
    
    #LAWS
    pis = pis_2D.func(x,tau=0.5, h)
    
    de.law = law.func(pv, pis, al)
    thetaS.hat.law6 = matrix(de.law, dims[1], dims[2])
    
    
    #SABHA
    de.sab = sab.func(pv, pis, al)
    thetaS.hat.sab6 = matrix(de.sab, dims[1], dims[2])
    


## ------------------------------------------------------------------------------------------
par(mfrow=c(2, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

image(thetaS, main="True States, mu=3")
image(thetaS.hat.new6, main="CLAW")
image(thetaS.hat.bh6, main="BH")
image(thetaS.hat.ad6, main="AdaDetect")
image(thetaS.hat.law6, main="LAWS")
image(thetaS.hat.sab6, main="SABHA")


## ------------------------------------------------------------------------------------------
par(mfrow=c(3, 6), mgp=c(2, 0.5, 0), mar=c(1.5, 1.5, 1, 1))
image(thetaS, main="True States, mu=1.5")
image(thetaS.hat.new1, main="CLAW")
image(thetaS.hat.bh1, main="BH")
image(thetaS.hat.ad1, main="AdaDetect")
image(thetaS.hat.law1, main="LAWS")
image(thetaS.hat.sab1, main="SABHA")
image(thetaS, main="True States, mu=2")
image(thetaS.hat.new3, main="CLAW")
image(thetaS.hat.bh3, main="BH")
image(thetaS.hat.ad3, main="AdaDetect")
image(thetaS.hat.law3, main="LAWS")
image(thetaS.hat.sab3, main="SABHA")
image(thetaS, main="True States, mu=2.5")
image(thetaS.hat.new5, main="CLAW")
image(thetaS.hat.bh5, main="BH")
image(thetaS.hat.ad5, main="AdaDetect")
image(thetaS.hat.law5, main="LAWS")
image(thetaS.hat.sab5, main="SABHA")


## ------------------------------------------------------------------------------------------
################ CLAW vs PLIS : true HMM
source("PLIS_func.R")
library(doParallel)
library(foreach)

n = 3000
m = n
al = 0.05
mu = seq(2.5, 3.5, length.out=6)
pii = c(1,0) #initial distribution
A = matrix(c(0.95, 0.05, 0.5, 0.5), 2, 2, byrow=T) #trans. matrix
#B = matrix(c(0.95, 0.05, 0.1, 0.9),2,2,byrow = T) #a fake trans. matrix
B=A
pp = A[1,2]/(A[1,2]+A[2,1]) #stationay prob. of non-null
f0 = c(0,1)
h = 150
fdr.new = NA
tdp.new = NA
fdr.2 = NA
tdp.2 = NA

for (k in 1:length(mu) ){
  f1 = c(mu[k],1)
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.2 = NA
  subtdp.2 = NA
  for (l in 1:200){
    f1 = c(mu[k],1)
    
    rdata = rdata1.hmm(n, pii, A, f0, f1)
    
    x = rdata$o
    a = rdata$s
    
    mir = rnorm(m, mean=f0[1], sd=f0[2])
    xx = NA
    for (q in 1:m) {
      if (abs(x[q]) >= abs(mir[q])){
        xx[q] = x[q]
      }else{
        xx[q] = mir[q]
      }
    }
    
    ####PLIS
    hmmest = em1.hmm(xx,maxiter=200)
    estf1 = hmmest$f1
    estA = hmmest$A
    
    cores = 10
    cl = makeCluster(cores)
    registerDoParallel(cl, cores = cores)
    two.plis = foreach(s=1:m, .combine = cbind) %dopar% Func2(s, xx, x, mir)
    stopCluster(cl)
    estlisx = two.plis[1,]
    estlism = two.plis[2,]
    
    estlis.de = conf.q(estlisx,estlism,al)
    subfdr.2[l] = sum(estlis.de*(1-a))/ max(sum(estlis.de),1)
    subtdp.2[l] = sum(estlis.de*a)/ max(sum(a),1)
    
    ###### CLAW
    pwscore = cor_1D.func(x, mir, tau=0.2, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-a)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*a) / max(1, sum(a) )
    
    cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.2[k] = mean(subfdr.2)
  tdp.2[k] = mean(subtdp.2)
  
}


ex2.2fdr = data.frame(mu, CLAW=fdr.new, PLIS=fdr.2)
ex2.2fdr = melt(ex2.2fdr,id="mu")
colnames(ex2.2fdr) = c("mu","METHOD","value")

ex2.2tdp = data.frame(mu, CLAW=tdp.new, PLIS=tdp.2)
ex2.2tdp = melt(ex2.2tdp,id="mu")
colnames(ex2.2tdp) = c("mu","METHOD","value")


line2.2fdr = ggplot(ex2.2fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0,0.2)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line2.2tdp = ggplot(ex2.2tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


print(line2.2fdr)
print(line2.2tdp)


################ CLAW vs PLIS : sin, not HMM
#pi = 0.4*[1+sin(0.2t)] for [201,500], [1501,1800], [801,1100], [2101,2400]
#pi = 0.02 otherwise
#f0=N(0,1), f1=N(nu,1) where  nu=mu+0.2sin(0.6t) 

m = 3000
al = 0.05
f0 = c(0,1)
h = 150 #weights bandwidth

pii = c(1,0)

indx = 1:m
ppvec = rep(0.02,m)
pdex = 0.4*(1+sin(0.2*indx))
ppvec[201:500]<-pdex[201:500]
ppvec[1501:1800]<-pdex[1501:1800]
ppvec[801:1100]<-pdex[801:1100]
ppvec[2101:2400]<-pdex[2101:2400]

muvec = seq(3,4,length.out=6)

fdr.new = NA
tdp.new = NA
fdr.2 = NA
tdp.2 = NA

for (k in 1:length(muvec) ){
  mu = muvec[k]
  mumu = mu+0.2*sin(0.6*indx) 
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.2 = NA
  subtdp.2 = NA
  for (l in 1:200){
    theta = rbinom(m, size=1, prob=ppvec)
    x = (1-theta)*rnorm(m,f0[1],f0[2]) + theta*rnorm(m,mumu,1)
    y = rnorm(m,f0[1],f0[2])
    
    mir = y
    xx = NA
    for (q in 1:m) {
      if (abs(x[q]) >= abs(mir[q])){
        xx[q] = x[q]
      }else{
        xx[q] = mir[q]
      }
    }
    
    ####PLIS
    hmmest = em1.hmm(xx,maxiter=200)
    estf1 = hmmest$f1
    estA = hmmest$A
    
    cores = 10
    cl = makeCluster(cores)
    registerDoParallel(cl, cores = cores)
    two.plis = foreach(s=1:m, .combine = cbind) %dopar% Func2(s, xx, x, mir)
    stopCluster(cl)
    estlisx = two.plis[1,]
    estlism = two.plis[2,]
    
    estlis.de = conf.q(estlisx,estlism,al)
    subfdr.2[l] = sum(estlis.de*(1-theta))/ max(sum(estlis.de),1)
    subtdp.2[l] = sum(estlis.de*theta)/ max(sum(theta),1)
    
    ###### CLAW
    pwscore = cor_1D.func(x, mir, tau=0.2, h)
    s = pwscore$drx
    st = pwscore$dry
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.2[k] = mean(subfdr.2)
  tdp.2[k] = mean(subtdp.2)
  
}


ex2.3fdr = data.frame(muvec, CLAW=fdr.new, PLIS=fdr.2)
ex2.3fdr = melt(ex2.3fdr,id="muvec")
colnames(ex2.3fdr) = c("mu","METHOD","value")

ex2.3tdp = data.frame(muvec, CLAW=tdp.new, PLIS=tdp.2)
ex2.3tdp = melt(ex2.3tdp,id="muvec")
colnames(ex2.3tdp) = c("mu","METHOD","value")


line2.3fdr = ggplot(ex2.3fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0,0.2)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line2.3tdp = ggplot(ex2.3tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


print(line2.3fdr)
print(line2.3tdp)


#############. plot
grid_arrange_shared_legend(line2.1fdr,line2.1tdp,nrow=1,ncol=2,position = "top")

grid_arrange_shared_legend(line2.2fdr,line2.3fdr,
                           line2.2tdp,line2.3tdp,nrow=2,ncol=2,position = "top")








################ CLAW vs PLIS : two-group model
source("PLIS_func.R")
library(doParallel)
library(foreach)

n = 3000
m = n
al = 0.05
mu = seq(3, 3.5, length.out=6)
pii = c(1,0) #initial distribution
A = matrix(c(0.9, 0.1, 0.9, 0.1), 2, 2, byrow=T) #trans. matrix
#B = matrix(c(0.95, 0.05, 0.1, 0.9),2,2,byrow = T) #a fake trans. matrix
B=A
pp = A[1,2]/(A[1,2]+A[2,1]) #stationay prob. of non-null
f0 = c(0,1)
h = 150
fdr.new = NA
tdp.new = NA
fdr.2 = NA
tdp.2 = NA

for (k in 1:length(mu) ){
  f1 = c(mu[k],1)
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.2 = NA
  subtdp.2 = NA
  for (l in 1:200){
    f1 = c(mu[k],1)
    
    rdata = rdata1.hmm(n, pii, A, f0, f1)
    
    x = rdata$o
    a = rdata$s
    
    mir = rnorm(m, mean=f0[1], sd=f0[2])
    xx = NA
    for (q in 1:m) {
      if (abs(x[q]) >= abs(mir[q])){
        xx[q] = x[q]
      }else{
        xx[q] = mir[q]
      }
    }
    
    ####PLIS
    hmmest = em1.hmm(xx,maxiter=200)
    estf1 = hmmest$f1
    estA = hmmest$A
    
    cores = 10
    cl = makeCluster(cores)
    registerDoParallel(cl, cores = cores)
    two.plis = foreach(s=1:m, .combine = cbind) %dopar% Func2(s, xx, x, mir)
    stopCluster(cl)
    estlisx = two.plis[1,]
    estlism = two.plis[2,]
    
    estlis.de = conf.q(estlisx,estlism,al)
    subfdr.2[l] = sum(estlis.de*(1-a))/ max(sum(estlis.de),1)
    subtdp.2[l] = sum(estlis.de*a)/ max(sum(a),1)
    
    ###### CLAW
    pwscore = exscore.func(x,mir,f0)#cor_1D.func(x, mir, tau=0.2, h)
    s = pwscore$s
    st = pwscore$st
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-a)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*a) / max(1, sum(a) )
    
    cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.2[k] = mean(subfdr.2)
  tdp.2[k] = mean(subtdp.2)
  
}


ex2.4fdr = data.frame(mu, CLAW=fdr.new, PLIS=fdr.2)
ex2.4fdr = melt(ex2.4fdr,id="mu")
colnames(ex2.4fdr) = c("mu","METHOD","value")

ex2.4tdp = data.frame(mu, CLAW=tdp.new, PLIS=tdp.2)
ex2.4tdp = melt(ex2.4tdp,id="mu")
colnames(ex2.4tdp) = c("mu","METHOD","value")


line2.4fdr = ggplot(ex2.4fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0,0.2)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line2.4tdp = ggplot(ex2.4tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


print(line2.4fdr)
print(line2.4tdp)



grid_arrange_shared_legend(line2.4fdr,line2.4tdp,nrow=1,ncol=2,position = "top")

grid_arrange_shared_legend(line2.2fdr,line2.3fdr,line2.4fdr,
                           line2.2tdp,line2.3tdp,line2.4tdp,nrow=2,ncol=3,position = "top")


################ CLAW vs PLIS : two-group model. (PLIS_TG)

n = 3000
m = n
al = 0.05
mu = seq(3.6, 3.8, length.out=6)
pii = c(1,0) #initial distribution
A = matrix(c(0.8, 0.2, 0.8, 0.2), 2, 2, byrow=T) #trans. matrix
#B = matrix(c(0.95, 0.05, 0.1, 0.9),2,2,byrow = T) #a fake trans. matrix
B=A
pp = A[1,2]/(A[1,2]+A[2,1]) #stationay prob. of non-null
f0 = c(0,1)
h = 150
fdr.new = NA
tdp.new = NA
fdr.2 = NA
tdp.2 = NA

for (k in 1:length(mu) ){
  f1 = c(mu[k],1)
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.2 = NA
  subtdp.2 = NA
  for (l in 1:200){
    f1 = c(mu[k],1)
    
    rdata = rdata1.hmm(n, pii, A, f0, f1)
    
    x = rdata$o
    a = rdata$s
    
    mir = rnorm(m, mean=f0[1], sd=f0[2])
    xx = NA
    for (q in 1:m) {
      if (abs(x[q]) >= abs(mir[q])){
        xx[q] = x[q]
      }else{
        xx[q] = mir[q]
      }
    }
    
    ####PLIS
    zv.ds = density(xx, from=min(x)-10, to=max(x)+10, n=1000)
    x.ds = lin.itp(x, zv.ds$x, zv.ds$y)
    estdrx = dnorm(x,f0[1],f0[2])/x.ds
    estdrx = estdrx/(1+estdrx)
    m.ds = lin.itp(mir, zv.ds$x, zv.ds$y)
    estdrm = dnorm(mir,0,1)/ m.ds
    estdrm = estdrm/(1+estdrm)
    
    estlis.de = conf.q(estdrx,estdrm,al)
    subfdr.2[l] = sum(estlis.de*(1-a))/ max(sum(estlis.de),1)
    subtdp.2[l] = sum(estlis.de*a)/ max(sum(a),1)
    
    ###### CLAW
    pwscore = exscore.func(x,mir,f0)#cor_1D.func(x, mir, tau=0.2, h)
    s = pwscore$s
    st = pwscore$st
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-a)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*a) / max(1, sum(a) )
    

  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.2[k] = mean(subfdr.2)
  tdp.2[k] = mean(subtdp.2)
  
}


ex2.5fdr = data.frame(mu, CLAW=fdr.new, PLIS=fdr.2)
ex2.5fdr = melt(ex2.5fdr,id="mu")
colnames(ex2.5fdr) = c("mu","METHOD","value")

ex2.5tdp = data.frame(mu, CLAW=tdp.new, PLIS=tdp.2)
ex2.5tdp = melt(ex2.5tdp,id="mu")
colnames(ex2.5tdp) = c("mu","METHOD","value")


line2.5fdr = ggplot(ex2.5fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0,0.2)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line2.5tdp = ggplot(ex2.5tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


print(line2.5fdr)
print(line2.5tdp)



grid_arrange_shared_legend(line2.5fdr,line2.5tdp,nrow=1,ncol=2,position = "top")



