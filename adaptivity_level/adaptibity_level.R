## ----setup, include=FALSE------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----include=FALSE-------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(lemon)
library(doParallel)
library(foreach)
library(adaptMT)
library(splines)
library(IHW)


## ------------------------------------------------------------------------------------------
source("group_func_storey.R")


## ----warning=FALSE-------------------------------------------------------------------------
m = 4500

m1 = 3000
m2 = 1500

al = 0.05

f0 = c(0,1)

fdr.new = NA
tdp.new = NA
fdr.pbh = NA
tdp.pbh = NA
fdr.sbh = NA
tdp.sbh = NA
fdr.pad = NA
tdp.pad = NA
fdr.sad = NA
tdp.sad = NA
fdr.ihw = NA
tdp.ihw = NA

muvec = seq(2,4,length.out=9)

for (k in 1:length(muvec) ){
  mu = muvec[k]
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.pbh = NA
  subtdp.pbh = NA
  subfdr.sbh = NA
  subtdp.sbh = NA
  subfdr.pad = NA
  subtdp.pad = NA
  subfdr.sad = NA
  subtdp.sad = NA
  subfdr.ihw = NA
  subtdp.ihw = NA
  
  for (l in 1:200){
    theta1 = rbinom(m1,1,0.5)
    theta2 = rbinom(m2,1,0.1)
    theta = c(theta1,theta2)
    
    x1 = (1-theta1)*rnorm(m1,f0[1],f0[2]) + theta1*rnorm(m1,mu,1)
    x2 = (1-theta2)*rnorm(m2,f0[1],f0[2]) + theta2*rnorm(m2,-2,0.5)
    x = c(x1,x2)
    
    y1 = rnorm(m1,f0[1],f0[2])
    y2 = rnorm(m2,f0[1],f0[2])
    y = c(y1,y2)
    
    score.g1 = group.func(x1,y1,f0)
    score.g2 = group.func(x2,y2,f0)
    
    s1 = score.g1$s
    st1 = score.g1$st
    s2 = score.g2$s
    st2 = score.g2$st
    
    s = c(s1,s2)
    st = c(st1,st2)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #pooled BH
    de.pbh = bh(x,al)
    subfdr.pbh[l] = sum(de.pbh*(1-theta)) / max(1, sum(de.pbh) )
    subtdp.pbh[l] = sum(de.pbh*theta) / max(1, sum(theta) )
    
    #BH-storey
    pv.s = 2*(1-pnorm(abs(x),0,1))
    pi01 = pi0.storey(pv.s,0.1)
    de.sbh = bh(x,al/pi01)
    subfdr.sbh[l] = sum(de.sbh*(1-theta)) / max(1, sum(de.sbh) )
    subtdp.sbh[l] = sum(de.sbh*theta) / max(1, sum(theta) )
    
    #pooled Adadetect
    score.pool = group.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.pad = Adadetect(pools,poolst,al)
    subfdr.pad[l] = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
    subtdp.pad[l] = sum(de.pad*theta) / max(1, sum(theta) )
    
    #Adadetect-storey
    pv.a = build.conf.p(pools,poolst)
    pi02 = pi0.storey(pv.a,0.1)
    de.sad = Adadetect(pools,poolst,al/pi02)
    subfdr.sad[l] = sum(de.sad*(1-theta)) / max(1, sum(de.sad) )
    subtdp.sad[l] = sum(de.sad*theta) / max(1, sum(theta) )
    
    
    
    #cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.pbh[k] = mean(subfdr.pbh)
  tdp.pbh[k] = mean(subtdp.pbh)
  fdr.sbh[k] = mean(subfdr.sbh)
  tdp.sbh[k] = mean(subtdp.sbh)
  fdr.pad[k] = mean(subfdr.pad)
  tdp.pad[k] = mean(subtdp.pad)
  fdr.sad[k] = mean(subfdr.sad)
  tdp.sad[k] = mean(subtdp.sad)
  fdr.ihw[k] = mean(subfdr.ihw)
  tdp.ihw[k] = mean(subtdp.ihw)
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.4fdr = data.frame(muvec, CLAW = fdr.new, BH=fdr.pbh, BH_Storey=fdr.sbh,
                      AdaDetect = fdr.pad, AdaDetect_Storey = fdr.sad)
ex1.4fdr = melt(ex1.4fdr,id="muvec")
colnames(ex1.4fdr) = c("mu","METHOD","value")


ex1.4tdp = data.frame(muvec, CLAW = tdp.new, BH=tdp.pbh, BH_Storey=tdp.sbh,
                      AdaDetect = tdp.pad, AdaDetect_Storey = tdp.sad)
ex1.4tdp = melt(ex1.4tdp,id="muvec")
colnames(ex1.4tdp) = c("mu","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.4fdr = ggplot(ex1.4fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0.02, 0.08)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.4tdp = ggplot(ex1.4tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.4fdr)
print(line1.4tdp)


## ------------------------------------------------------------------------------------------
#grid_arrange_shared_legend(line1.1fdr,line1.1tdp,nrow=1,ncol=2,position = "top")


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
ppvec[201:350]<-0.9
ppvec[1501:1650]<-0.9
ppvec[801:1000]<-0.6
ppvec[2101:2300]<-0.6
muvec = seq(2.1,3.6,length.out=9)

fdr.new = NA
tdp.new = NA
fdr.bh = NA
tdp.bh = NA
fdr.sbh = NA
tdp.sbh = NA
fdr.ad = NA
tdp.ad = NA
fdr.sad = NA
tdp.sad = NA

for (k in 1:length(muvec) ){
  mu = muvec[k]
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.bh = NA
  subtdp.bh = NA
  subfdr.sbh = NA
  subtdp.sbh = NA
  subfdr.ad = NA
  subtdp.ad = NA
  subfdr.sad = NA
  subtdp.sad = NA

  
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
    
    
    #BH-storey
    pi01 = pi0.storey(pv,0.1)
    de.sbh = bh(pv,al/pi01)
    subfdr.sbh[l] = sum(de.sbh*(1-theta)) / max(1, sum(de.sbh) )
    subtdp.sbh[l] = sum(de.sbh*theta) / max(1, sum(theta) )
    
    
    #Adadetect-storey
    pv.a = build.conf.p(pools,poolst)
    pi02 = pi0.storey(pv.a,0.1)
    de.sad = Adadetect(pools,poolst,al/pi02)
    subfdr.sad[l] = sum(de.sad*(1-theta)) / max(1, sum(de.sad) )
    subtdp.sad[l] = sum(de.sad*theta) / max(1, sum(theta) )
    
    
    #cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.bh[k] = mean(subfdr.bh)
  tdp.bh[k] = mean(subtdp.bh)
  fdr.ad[k] = mean(subfdr.ad)
  tdp.ad[k] = mean(subtdp.ad)
  fdr.sbh[k] = mean(subfdr.sbh)
  tdp.sbh[k] = mean(subtdp.sbh)
  fdr.sad[k] = mean(subfdr.sad)
  tdp.sad[k] = mean(subtdp.sad)
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.1fdr = data.frame(muvec, CLAW = fdr.new, BH=fdr.bh, BH_Storey=fdr.sbh,
                      AdaDetect = fdr.ad, AdaDetect_Storey = fdr.sad)
ex1.1fdr = melt(ex1.1fdr,id="muvec")
colnames(ex1.1fdr) = c("mu","METHOD","value")


ex1.1tdp = data.frame(muvec, CLAW = tdp.new, BH=tdp.bh, BH_Storey=tdp.sbh,
                      AdaDetect = tdp.ad, AdaDetect_Storey = tdp.sad)
ex1.1tdp = melt(ex1.1tdp,id="muvec")
colnames(ex1.1tdp) = c("mu","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.1fdr = ggplot(ex1.1fdr)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0.02, 0.08)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.1tdp = ggplot(ex1.1tdp)+aes(x=mu, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.1fdr)
print(line1.1tdp)


## ----message=FALSE, warning=FALSE----------------------------------------------------------
grid_arrange_shared_legend(line1.4fdr,line1.1fdr,
                           line1.4tdp,line1.1tdp,nrow=2,ncol=2,position = "top")


## ------------------------------------------------------------------------------------------
source("group_func_storey.R")


## ----warning=FALSE-------------------------------------------------------------------------
m = 4500

m1 = 3000
m2 = 1500

al = 0.05

f0 = c(0,1)

fdr.new = NA
tdp.new = NA
fdr.pbh = NA
tdp.pbh = NA
fdr.sbh = NA
tdp.sbh = NA
fdr.pad = NA
tdp.pad = NA
fdr.sad = NA
tdp.sad = NA
fdr.ihw = NA
tdp.ihw = NA

#muvec = seq(2,4,length.out=9)
pivec = seq(0.1,0.6,length.out=9)

for (k in 1:length(pivec) ){
  mu = 3.6
  pp = pivec[k]
  
  subfdr.new = NA
  subtdp.new = NA
  subfdr.pbh = NA
  subtdp.pbh = NA
  subfdr.sbh = NA
  subtdp.sbh = NA
  subfdr.pad = NA
  subtdp.pad = NA
  subfdr.sad = NA
  subtdp.sad = NA
  subfdr.ihw = NA
  subtdp.ihw = NA
  
  for (l in 1:200){
    theta1 = rbinom(m1,1,pp)
    theta2 = rbinom(m2,1,0.1)
    theta = c(theta1,theta2)
    
    x1 = (1-theta1)*rnorm(m1,f0[1],f0[2]) + theta1*rnorm(m1,mu,1)
    x2 = (1-theta2)*rnorm(m2,f0[1],f0[2]) + theta2*rnorm(m2,-2,0.5)
    x = c(x1,x2)
    
    y1 = rnorm(m1,f0[1],f0[2])
    y2 = rnorm(m2,f0[1],f0[2])
    y = c(y1,y2)
    
    score.g1 = group.func(x1,y1,f0)
    score.g2 = group.func(x2,y2,f0)
    
    s1 = score.g1$s
    st1 = score.g1$st
    s2 = score.g2$s
    st2 = score.g2$st
    
    s = c(s1,s2)
    st = c(st1,st2)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
    
    #pooled BH
    de.pbh = bh(x,al)
    subfdr.pbh[l] = sum(de.pbh*(1-theta)) / max(1, sum(de.pbh) )
    subtdp.pbh[l] = sum(de.pbh*theta) / max(1, sum(theta) )
    
    #BH-storey
    pv.s = 2*(1-pnorm(abs(x),0,1))
    pi01 = pi0.storey(pv.s,0.1)
    de.sbh = bh(x,al/pi01)
    subfdr.sbh[l] = sum(de.sbh*(1-theta)) / max(1, sum(de.sbh) )
    subtdp.sbh[l] = sum(de.sbh*theta) / max(1, sum(theta) )
    
    #pooled Adadetect
    score.pool = group.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.pad = Adadetect(pools,poolst,al)
    subfdr.pad[l] = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
    subtdp.pad[l] = sum(de.pad*theta) / max(1, sum(theta) )
    
    #Adadetect-storey
    pv.a = build.conf.p(pools,poolst)
    pi02 = pi0.storey(pv.a,0.1)
    de.sad = Adadetect(pools,poolst,al/pi02)
    subfdr.sad[l] = sum(de.sad*(1-theta)) / max(1, sum(de.sad) )
    subtdp.sad[l] = sum(de.sad*theta) / max(1, sum(theta) )
    
    
    
    #cat(k, l, "\n")
  }
  fdr.new[k] = mean(subfdr.new)
  tdp.new[k] = mean(subtdp.new)
  fdr.pbh[k] = mean(subfdr.pbh)
  tdp.pbh[k] = mean(subtdp.pbh)
  fdr.sbh[k] = mean(subfdr.sbh)
  tdp.sbh[k] = mean(subtdp.sbh)
  fdr.pad[k] = mean(subfdr.pad)
  tdp.pad[k] = mean(subtdp.pad)
  fdr.sad[k] = mean(subfdr.sad)
  tdp.sad[k] = mean(subtdp.sad)
  fdr.ihw[k] = mean(subfdr.ihw)
  tdp.ihw[k] = mean(subtdp.ihw)
}


## ----warning=FALSE-------------------------------------------------------------------------
ex1.5fdr = data.frame(pivec, CLAW = fdr.new, BH=fdr.pbh, BH_Storey=fdr.sbh,
                      AdaDetect = fdr.pad, AdaDetect_Storey = fdr.sad)
ex1.5fdr = melt(ex1.5fdr,id="pivec")
colnames(ex1.5fdr) = c("pi","METHOD","value")


ex1.5tdp = data.frame(pivec, CLAW = tdp.new, BH=tdp.pbh, BH_Storey=tdp.sbh,
                      AdaDetect = tdp.pad, AdaDetect_Storey = tdp.sad)
ex1.5tdp = melt(ex1.5tdp,id="pivec")
colnames(ex1.5tdp) = c("pi","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.5fdr = ggplot(ex1.5fdr)+aes(x=pi, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0.02, 0.08)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.5tdp = ggplot(ex1.5tdp)+aes(x=pi, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.5fdr)
print(line1.5tdp)


## ----message=FALSE, warning=FALSE----------------------------------------------------------
grid_arrange_shared_legend(line1.4fdr,line1.1fdr,line1.5fdr,
                           line1.4tdp,line1.1tdp,line1.5tdp,nrow=2,ncol=3,position = "top")

