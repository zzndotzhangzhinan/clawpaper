
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
library(Compositional)


## ----warning=FALSE-------------------------------------------------------------------------
m = 3000 
al = 0.05
f0 = c(0,1)
#h = 150 #weights bandwidth
#ss = rbeta(m,2,5)

ppvec = 0.1#rep(0.02, m)
#ppvec[201:350]<-0.6
#ppvec[1501:1650]<-0.6
#ppvec[801:1000]<-0.3
#ppvec[2101:2300]<-0.3

muvec = seq(1.8,2.4,length.out=6)

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
    mux = (1-theta)*(0)+theta*(1)
    muy = (1-theta)*(0)+theta*(-mu)
    x1 = rnorm(m,mux,1)
    x2 = rnorm(m,muy,0.7)
    sigp = sqrt(0.5*1^2+0.5*0.7^2)
    kapp = 1^2/(0.7^2)
    x = sqrt(0.5)*(x1-x2)/sigp
    ss = sqrt(0.5)*(x1+kapp*x2)/(sqrt(kapp)*sigp)
    y = rnorm(m,f0[1],f0[2])
    
    
    pv = 2*(1-pnorm(abs(x),0,1))
    pis = pis_1D.func(pv,tau=0.5, ss)
    
    
    #new mwthod
    pvy = 2*(1-pnorm(abs(y),0,1))
    pisy = pis_1D.func(pvy,tau=0.5, ss)
    piclaw = (pis+pisy)/2
    pdf1 = dnorm(c(x,y),0,1)
    sdensity = density(ss,from=min(ss)-10, to=max(ss)+10, n=1000)
    pdf2 = lin.itp(ss, sdensity$x, sdensity$y)
    mat = matrix(c(x,y,ss,ss), ncol=2, byrow = F)
    bidens = mkde(mat)
    s = (1-piclaw)*pdf1[1:m]*pdf2/bidens[1:m]
    st = (1-piclaw)*pdf1[-(1:m)]*pdf2/bidens[-(1:m)]

    
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #BH
    #pv = 2*(1-pnorm(abs(x),0,1))
    
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
    #pis = pis_1D.func(pv,tau=0.5, ss)
    
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
    adapt.res <- adapt(x = data.frame(x = ss), pvals = pvals, models = models, 
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
#h = 150 #weights bandwidth
#ss = rbeta(m,2,5)

pivec = seq(0.1,0.2,length.out=6)#rep(0.02, m)
#ppvec[201:350]<-0.6
#ppvec[1501:1650]<-0.6
#ppvec[801:1000]<-0.3
#ppvec[2101:2300]<-0.3

#muvec = seq(1.6,3,length.out=6)

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
  #mu = muvec[k]
  
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
    theta = rbinom(m, size=1, prob=pivec[k])
    mux = (1-theta)*(0)+theta*(1)
    muy = (1-theta)*(0)+theta*(-2.2)
    x1 = rnorm(m,mux,1)
    x2 = rnorm(m,muy,0.7)
    sigp = sqrt(0.5*1^2+0.5*0.7^2)
    kapp = 1^2/(0.7^2)
    x = sqrt(0.5)*(x1-x2)/sigp
    ss = sqrt(0.5)*(x1+kapp*x2)/(sqrt(kapp)*sigp)
    y = rnorm(m,f0[1],f0[2])
    
    
    pv = 2*(1-pnorm(abs(x),0,1))
    pis = pis_1D.func(pv,tau=0.5, ss)
    
    
    #new mwthod
    pvy = 2*(1-pnorm(abs(y),0,1))
    pisy = pis_1D.func(pvy,tau=0.5, ss)
    piclaw = (pis+pisy)/2
    pdf1 = dnorm(c(x,y),0,1)
    sdensity = density(ss,from=min(ss)-10, to=max(ss)+10, n=1000)
    pdf2 = lin.itp(ss, sdensity$x, sdensity$y)
    mat = matrix(c(x,y,ss,ss), ncol=2, byrow = F)
    bidens = mkde(mat)
    zero.ind = which(bidens==0)
    bidens[zero.ind] = min(bidens[-zero.ind])
    s = (1-piclaw)*pdf1[1:m]*pdf2/bidens[1:m]
    st = (1-piclaw)*pdf1[-(1:m)]*pdf2/bidens[-(1:m)]

    
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #BH
    #pv = 2*(1-pnorm(abs(x),0,1))
    
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
    #pis = pis_1D.func(pv,tau=0.5, ss)
    
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
    adapt.res <- adapt(x = data.frame(x = ss), pvals = pvals, models = models, 
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
#h = 150 #weights bandwidth
#ss = rbeta(m,2,5)

#pivec = seq(0.1,0.4,length.out=6)#rep(0.02, m)
#ppvec[201:350]<-0.6
#ppvec[1501:1650]<-0.6
#ppvec[801:1000]<-0.3
#ppvec[2101:2300]<-0.3

#muvec = seq(1.6,3,length.out=6)
NNvec = c(500,600,700,800,900,1000)

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

for (k in 1:length(NNvec) ){
  NN = NNvec[k]
  
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
    #theta = rbinom(m, size=1, prob=0.2)
    mux = c(rep(0,800),rep(3,200),rep(3,NN),rep(0,m-1000-NN))
    muy = c(rep(0,1000),rep(-0.5,1000),rep(0,1000))
    theta = rep(0,m)
    theta[mux!=muy] = 1

    x1 = rnorm(m,mux,1)
    x2 = rnorm(m,muy,0.7)
    sigp = sqrt(0.5*1^2+0.5*0.7^2)
    kapp = 1^2/(0.7^2)
    x = sqrt(0.5)*(x1-x2)/sigp
    ss = sqrt(0.5)*(x1+kapp*x2)/(sqrt(kapp)*sigp)
    y = rnorm(m,f0[1],f0[2])
    
    
    pv = 2*(1-pnorm(abs(x),0,1))
    pis = pis_1D.func(pv,tau=0.5, ss)
    
    
    #new mwthod
    pvy = 2*(1-pnorm(abs(y),0,1))
    pisy = pis_1D.func(pvy,tau=0.5, ss)
    piclaw = (pis+pisy)/2
    pdf1 = dnorm(c(x,y),0,1)
    sdensity = density(ss,from=min(ss)-10, to=max(ss)+10, n=1000)
    pdf2 = lin.itp(ss, sdensity$x, sdensity$y)
    mat = matrix(c(x,y,ss,ss), ncol=2, byrow = F)
    bidens = mkde(mat)
    s = (1-piclaw)*pdf1[1:m]*pdf2/bidens[1:m]
    st = (1-piclaw)*pdf1[-(1:m)]*pdf2/bidens[-(1:m)]

    
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #BH
    #pv = 2*(1-pnorm(abs(x),0,1))
    
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
    #pis = pis_1D.func(pv,tau=0.5, ss)
    
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
    adapt.res <- adapt(x = data.frame(x = ss), pvals = pvals, models = models, 
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
ex1.3fdr = data.frame(NNvec, CLAW = fdr.new, BH=fdr.bh, AdaDetect=fdr.ad,
                      LAWS = fdr.law, SABHA = fdr.sab, AdaPT=fdr.adapt)
ex1.3fdr = melt(ex1.3fdr,id="NNvec")
colnames(ex1.3fdr) = c("N","METHOD","value")


ex1.3tdp = data.frame(NNvec, CLAW = tdp.new, BH=tdp.bh, AdaDetect=tdp.ad,
                      LAWS = tdp.law, SABHA = tdp.sab, AdaPT=tdp.adapt)
ex1.3tdp = melt(ex1.3tdp,id="NNvec")
colnames(ex1.3tdp) = c("N","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.3fdr = ggplot(ex1.3fdr)+aes(x=N, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.3tdp = ggplot(ex1.3tdp)+aes(x=N, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.3fdr)
print(line1.3tdp)


## ------------------------------------------------------------------------------------------
grid_arrange_shared_legend(line1.1fdr,line1.2fdr,line1.3fdr,
                           line1.1tdp,line1.2tdp,line1.3tdp,nrow=2,ncol=3,position = "top")

