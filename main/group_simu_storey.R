
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

muvec = seq(3.8,5,length.out=11)

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
    theta1 = rbinom(m1,1,0.2)
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
    
    #separate BH
    de.sbh1 = bh(x1,al)
    de.sbh2 = bh(x2,al)
    de.sbh = c(de.sbh1,de.sbh2)
    subfdr.sbh[l] = sum(de.sbh*(1-theta)) / max(1, sum(de.sbh) )
    subtdp.sbh[l] = sum(de.sbh*theta) / max(1, sum(theta) )
    
    #pooled Adadetect
    score.pool = group.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.pad = Adadetect(pools,poolst,al)
    subfdr.pad[l] = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
    subtdp.pad[l] = sum(de.pad*theta) / max(1, sum(theta) )
    
    #separate Adadetect
    de.sad1 = Adadetect(s1,st1,al)
    de.sad2 = Adadetect(s2,st2,al)
    de.sad = c(de.sad1,de.sad2)
    subfdr.sad[l] = sum(de.sad*(1-theta)) / max(1, sum(de.sad) )
    subtdp.sad[l] = sum(de.sad*theta) / max(1, sum(theta) )
    
    #IHW
    pv = 2*(1-pnorm(abs(x),0,1))
    ihw_fdr = ihw(pv, covariates=c(rep(1,m1),rep(2,m2)), alpha=al, nbins=2)
    ihw.de = (adj_pvalues(ihw_fdr)<=al)
    subfdr.ihw[l] = sum(ihw.de*(1-theta))/ max(sum(ihw.de),1)
    subtdp.ihw[l] = sum(ihw.de*theta)/ max(sum(theta),1)
    
    
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
ex1.4fdr = data.frame(muvec, CLAW = fdr.new, PooledBH=fdr.pbh, SeparateBH=fdr.sbh,
                      PooledAD = fdr.pad, SeparateAD = fdr.sad, IHW=fdr.ihw)
ex1.4fdr = melt(ex1.4fdr,id="muvec")
colnames(ex1.4fdr) = c("mu","METHOD","value")


ex1.4tdp = data.frame(muvec, CLAW = tdp.new, PooledBH=tdp.pbh, SeparateBH=tdp.sbh,
                      PooledAD = tdp.pad, SeparateAD = tdp.sad, IHW=tdp.ihw)
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


## ------------------------------------------------------------------------------------------
#grid_arrange_shared_legend(line1.1fdr,line1.1tdp,nrow=1,ncol=2,position = "top")


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

ppvec = seq(0.1,0.4,length.out=11)

for (k in 1:length(ppvec) ){
  pp = ppvec[k]
  
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
    theta1 = rbinom(m1,1,0.2)
    theta2 = rbinom(m2,1,pp)
    theta = c(theta1,theta2)
    
    x1 = (1-theta1)*rnorm(m1,f0[1],f0[2]) + theta1*rnorm(m1,2,1)
    x2 = (1-theta2)*rnorm(m2,f0[1],f0[2]) + theta2*rnorm(m2,-4,1)
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
    
    #separate BH
    de.sbh1 = bh(x1,al)
    de.sbh2 = bh(x2,al)
    de.sbh = c(de.sbh1,de.sbh2)
    subfdr.sbh[l] = sum(de.sbh*(1-theta)) / max(1, sum(de.sbh) )
    subtdp.sbh[l] = sum(de.sbh*theta) / max(1, sum(theta) )
    
    #pooled Adadetect
    score.pool = group.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.pad = Adadetect(pools,poolst,al)
    subfdr.pad[l] = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
    subtdp.pad[l] = sum(de.pad*theta) / max(1, sum(theta) )
    
    #separate Adadetect
    de.sad1 = Adadetect(s1,st1,al)
    de.sad2 = Adadetect(s2,st2,al)
    de.sad = c(de.sad1,de.sad2)
    subfdr.sad[l] = sum(de.sad*(1-theta)) / max(1, sum(de.sad) )
    subtdp.sad[l] = sum(de.sad*theta) / max(1, sum(theta) )
    
    #IHW
    pv = 2*(1-pnorm(abs(x),0,1))
    ihw_fdr = ihw(pv, covariates=c(rep(1,m1),rep(2,m2)), alpha=al, nbins=2)
    ihw.de = (adj_pvalues(ihw_fdr)<=al)
    subfdr.ihw[l] = sum(ihw.de*(1-theta))/ max(sum(ihw.de),1)
    subtdp.ihw[l] = sum(ihw.de*theta)/ max(sum(theta),1)
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
ex1.5fdr = data.frame(ppvec, CLAW = fdr.new, PooledBH=fdr.pbh, SeparateBH=fdr.sbh,
                      PooledAD = fdr.pad, SeparateAD = fdr.sad, IHW = fdr.ihw)
ex1.5fdr = melt(ex1.5fdr,id="ppvec")
colnames(ex1.5fdr) = c("p","METHOD","value")


ex1.5tdp = data.frame(ppvec, CLAW = tdp.new, PooledBH=tdp.pbh, SeparateBH=tdp.sbh,
                      PooledAD = tdp.pad, SeparateAD = tdp.sad, IHW = tdp.ihw)
ex1.5tdp = melt(ex1.5tdp,id="ppvec")
colnames(ex1.5tdp) = c("p","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.5fdr = ggplot(ex1.5fdr)+aes(x=p, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+xlab("pi")+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.5tdp = ggplot(ex1.5tdp)+aes(x=p, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("AP")+xlab("pi")+
  theme_bw()


## ------------------------------------------------------------------------------------------
print(line1.5fdr)
print(line1.5tdp)


## ----warning=FALSE-------------------------------------------------------------------------
#m = 3000

m2vec = seq(1000,2000,by=100)
m1 = 3000


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

for (k in 1:length(m2vec) ){
  m2 = m2vec[k]
  m = m1+m2
  
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
    theta1 = rbinom(m1,1,0.2)
    theta2 = rbinom(m2,1,0.1)
    theta = c(theta1,theta2)
    
    x1 = (1-theta1)*rnorm(m1,f0[1],f0[2]) + theta1*rnorm(m1,2,0.5)
    x2 = (1-theta2)*rnorm(m2,f0[1],f0[2]) + theta2*rnorm(m2,-4,1)
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
    
    #separate BH
    de.sbh1 = bh(x1,al)
    de.sbh2 = bh(x2,al)
    de.sbh = c(de.sbh1,de.sbh2)
    subfdr.sbh[l] = sum(de.sbh*(1-theta)) / max(1, sum(de.sbh) )
    subtdp.sbh[l] = sum(de.sbh*theta) / max(1, sum(theta) )
    
    #pooled Adadetect
    score.pool = group.func(x,y,f0)
    pools = score.pool$s
    poolst = score.pool$st
    de.pad = Adadetect(pools,poolst,al)
    subfdr.pad[l] = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
    subtdp.pad[l] = sum(de.pad*theta) / max(1, sum(theta) )
    
    #separate Adadetect
    de.sad1 = Adadetect(s1,st1,al)
    de.sad2 = Adadetect(s2,st2,al)
    de.sad = c(de.sad1,de.sad2)
    subfdr.sad[l] = sum(de.sad*(1-theta)) / max(1, sum(de.sad) )
    subtdp.sad[l] = sum(de.sad*theta) / max(1, sum(theta) )
    
    #IHW
    pv = 2*(1-pnorm(abs(x),0,1))
    ihw_fdr = ihw(pv, covariates=c(rep(1,m1),rep(2,m2)), alpha=al, nbins=2)
    ihw.de = (adj_pvalues(ihw_fdr)<=al)
    subfdr.ihw[l] = sum(ihw.de*(1-theta))/ max(sum(ihw.de),1)
    subtdp.ihw[l] = sum(ihw.de*theta)/ max(sum(theta),1)
    
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
ex1.6fdr = data.frame(m2vec, CLAW = fdr.new, PooledBH=fdr.pbh, SeparateBH=fdr.sbh,
                      PooledAD = fdr.pad, SeparateAD = fdr.sad, IHW=fdr.ihw)
ex1.6fdr = melt(ex1.6fdr,id="m2vec")
colnames(ex1.6fdr) = c("m2","METHOD","value")


ex1.6tdp = data.frame(m2vec, CLAW = tdp.new, PooledBH=tdp.pbh, SeparateBH=tdp.sbh,
                      PooledAD = tdp.pad, SeparateAD = tdp.sad, IHW=tdp.ihw)
ex1.6tdp = melt(ex1.6tdp,id="m2vec")
colnames(ex1.6tdp) = c("m2","METHOD","value")


## ------------------------------------------------------------------------------------------
line1.6fdr = ggplot(ex1.6fdr)+aes(x=m2, y=value, group=METHOD, color=METHOD, shape=METHOD)+
  geom_point(size=3)+
  geom_line(size=0.7)+ylab("FDR")+ylim(0, 0.1)+
  geom_hline(aes(yintercept=al),size=1)+
  theme_bw()

line1.6tdp = ggplot(ex1.6tdp)+aes(x=m2, y=value, group=METHOD, color=METHOD, shape=METHOD)+
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
library(mvtnorm)
library(Compositional)


## ----warning=FALSE-------------------------------------------------------------------------
set.seed(123456)

al = 0.05
d = 2
pp01 = 0.8
pp02 = 0.9

m1=1000
m2=2500
m=m1+m2

mu1 = rep(sqrt(2*log(d)),2)
mu2 = rep(2,d)

subfdr.1 = NA
subtdp.1 = NA
subfdr.new = NA
subtdp.new = NA

for (l in 1:200) {
  
theta1 = rbinom(m1,1,0.2)
theta2 = rbinom(m2,1,0.1)
theta = c(theta1,theta2)

x1 = (1-theta1)*rmvnorm(m1,rep(0,d),diag(d)) + theta1*rmvnorm(m1,mu1,diag(d))
x2 = (1-theta2)*rmvnorm(m2,rep(0,d),diag(d)) + theta2*rmvnorm(m2,mu2,diag(d))
x = rbind(x1,x2)
    
y1 = rmvnorm(m1,rep(0,d),diag(d))
y2 = rmvnorm(m2,rep(0,d),diag(d))
y = rbind(y1,y2)

clfdr1 = 0.8*dmvnorm(x1,rep(0,d),diag(d))/mkde(x1)
clfdr2 = 0.9*dmvnorm(x2,rep(0,d),diag(d))/mkde(x2)

clfdr = c(clfdr1,clfdr2)
rcl = sort(clfdr)
hk=1
while(hk<m && (1/hk)*sum(rcl[1:hk])<al){
  hk=hk+1
  }
hk = hk-1
rej = (clfdr<=rcl[hk])
subfdr.1[l] = sum(rej*(1-theta)) / max(1, sum(rej) )
subtdp.1[l] = sum(rej*theta) / max(1, sum(theta) )

#new mwthod
data1 = rbind(x1,y1)
g1ds = (1+0.8)/2 * dmvnorm(data1,rep(0,d),diag(d))/mkde(data1)
g1ds = pmin(g1ds,0.999)
g1ds = (0.5-0.1)/(1-0.1)*g1ds/(1-g1ds)
s1 = g1ds[1:m1]
st1 = g1ds[-(1:m1)]


data2 = rbind(x2,y2)
g2ds = (1+0.9)/2 * dmvnorm(data2,rep(0,d),diag(d))/mkde(data2)
g2ds = pmin(g2ds,0.999)
g2ds = (0.5-0.05)/(1-0.05)*g2ds/(1-g2ds)
s2 = g2ds[1:m2]
st2 = g2ds[-(1:m2)]

s = c(s1,s2)
st = c(st1,st2)

#new mwthod
    st = st/(1+st)
    s = s/(1+s)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
}


Dimension = rep("d1=2",400)
METHOD = c(rep("Clfdr",200), rep("CLAW",200))
FDR = c(subfdr.1,subfdr.new)
AP = c(subtdp.1,subtdp.new)
data.box.fdr = data.frame(Dimension, METHOD, FDR)
data.box.ap = data.frame(Dimension, METHOD, AP)


## ----warning=FALSE-------------------------------------------------------------------------
al = 0.05
d = 4
pp01 = 0.8
pp02 = 0.9

m1=1000
m2=2500
m=m1+m2

mu1 = rep(sqrt(2*log(d)),4)
mu2 = rep(2,d)

subfdr.1 = NA
subtdp.1 = NA
subfdr.new = NA
subtdp.new = NA

for (l in 1:200) {
  
theta1 = rbinom(m1,1,0.2)
theta2 = rbinom(m2,1,0.1)
theta = c(theta1,theta2)

x1 = (1-theta1)*rmvnorm(m1,rep(0,d),diag(d)) + theta1*rmvnorm(m1,mu1,diag(d))
x2 = (1-theta2)*rmvnorm(m2,rep(0,d),diag(d)) + theta2*rmvnorm(m2,mu2,diag(d))
x = rbind(x1,x2)
    
y1 = rmvnorm(m1,rep(0,d),diag(d))
y2 = rmvnorm(m2,rep(0,d),diag(d))
y = rbind(y1,y2)

clfdr1 = 0.8*dmvnorm(x1,rep(0,d),diag(d))/mkde(x1)
clfdr2 = 0.9*dmvnorm(x2,rep(0,d),diag(d))/mkde(x2)

clfdr = c(clfdr1,clfdr2)
rcl = sort(clfdr)
hk=1
while(hk<m && (1/hk)*sum(rcl[1:hk])<al){
  hk=hk+1
  }
hk = hk-1
rej = (clfdr<=rcl[hk])
subfdr.1[l] = sum(rej*(1-theta)) / max(1, sum(rej) )
subtdp.1[l] = sum(rej*theta) / max(1, sum(theta) )

#new mwthod
data1 = rbind(x1,y1)
g1ds = (1+0.8)/2 * dmvnorm(data1,rep(0,d),diag(d))/mkde(data1)
g1ds = pmin(g1ds,0.999)
g1ds = (0.5-0.1)/(1-0.1)*g1ds/(1-g1ds)
s1 = g1ds[1:m1]
st1 = g1ds[-(1:m1)]


data2 = rbind(x2,y2)
g2ds = (1+0.9)/2 * dmvnorm(data2,rep(0,d),diag(d))/mkde(data2)
g2ds = pmin(g2ds,0.999)
g2ds = (0.5-0.05)/(1-0.05)*g2ds/(1-g2ds)
s2 = g2ds[1:m2]
st2 = g2ds[-(1:m2)]


s = c(s1,s2)
st = c(st1,st2)

#new mwthod
    st = st/(1+st)
    s = s/(1+s)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
}

Dimension = rep("d2=4",400)
METHOD = c(rep("Clfdr",200), rep("CLAW",200))
FDR = c(subfdr.1,subfdr.new)
AP = c(subtdp.1,subtdp.new)
add.fdr = data.frame(Dimension, METHOD, FDR)
add.ap = data.frame(Dimension, METHOD, AP)

data.box.fdr = rbind(data.box.fdr,add.fdr)
data.box.ap = rbind(data.box.ap,add.ap)


## ----warning=FALSE-------------------------------------------------------------------------
al = 0.05
d = 8
pp01 = 0.8
pp02 = 0.9

m1=1000
m2=2500
m=m1+m2

mu1 = c(rep(sqrt(2*log(d)),4), rep(0,d-4))
mu2 = rep(2,d)

subfdr.1 = NA
subtdp.1 = NA
subfdr.new = NA
subtdp.new = NA

for (l in 1:200) {
  
theta1 = rbinom(m1,1,0.2)
theta2 = rbinom(m2,1,0.1)
theta = c(theta1,theta2)

x1 = (1-theta1)*rmvnorm(m1,rep(0,d),diag(d)) + theta1*rmvnorm(m1,mu1,diag(d))
x2 = (1-theta2)*rmvnorm(m2,rep(0,d),diag(d)) + theta2*rmvnorm(m2,mu2,diag(d))
x = rbind(x1,x2)
    
y1 = rmvnorm(m1,rep(0,d),diag(d))
y2 = rmvnorm(m2,rep(0,d),diag(d))
y = rbind(y1,y2)

clfdr1 = 0.8*dmvnorm(x1,rep(0,d),diag(d))/mkde(x1)
clfdr2 = 0.9*dmvnorm(x2,rep(0,d),diag(d))/mkde(x2)

clfdr = c(clfdr1,clfdr2)
rcl = sort(clfdr)
hk=1
while(hk<m && (1/hk)*sum(rcl[1:hk])<al){
  hk=hk+1
  }
hk = hk-1
rej = (clfdr<=rcl[hk])
subfdr.1[l] = sum(rej*(1-theta)) / max(1, sum(rej) )
subtdp.1[l] = sum(rej*theta) / max(1, sum(theta) )

#new mwthod
data1 = rbind(x1,y1)
g1ds = (1+0.8)/2 * dmvnorm(data1,rep(0,d),diag(d))/mkde(data1)
g1ds = pmin(g1ds,0.999)
g1ds = (0.5-0.1)/(1-0.1)*g1ds/(1-g1ds)
s1 = g1ds[1:m1]
st1 = g1ds[-(1:m1)]


data2 = rbind(x2,y2)
g2ds = (1+0.9)/2 * dmvnorm(data2,rep(0,d),diag(d))/mkde(data2)
g2ds = pmin(g2ds,0.999)
g2ds = (0.5-0.05)/(1-0.05)*g2ds/(1-g2ds)
s2 = g2ds[1:m2]
st2 = g2ds[-(1:m2)]

s = c(s1,s2)
st = c(st1,st2)

#new mwthod
    st = st/(1+st)
    s = s/(1+s)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
}

Dimension = rep("d3=8",400)
METHOD = c(rep("Clfdr",200), rep("CLAW",200))
FDR = c(subfdr.1,subfdr.new)
AP = c(subtdp.1,subtdp.new)
add.fdr = data.frame(Dimension, METHOD, FDR)
add.ap = data.frame(Dimension, METHOD, AP)

data.box.fdr = rbind(data.box.fdr,add.fdr)
data.box.ap = rbind(data.box.ap,add.ap)


## ----warning=FALSE-------------------------------------------------------------------------
al = 0.05
d = 16
pp01 = 0.8
pp02 = 0.9

m1=1000
m2=2500
m=m1+m2

mu1 = c(rep(sqrt(2*log(d)),4), rep(0,d-4))
mu2 = rep(2,d)

subfdr.1 = NA
subtdp.1 = NA
subfdr.new = NA
subtdp.new = NA

for (l in 1:200) {
  
theta1 = rbinom(m1,1,0.2)
theta2 = rbinom(m2,1,0.1)
theta = c(theta1,theta2)

x1 = (1-theta1)*rmvnorm(m1,rep(0,d),diag(d)) + theta1*rmvnorm(m1,mu1,diag(d))
x2 = (1-theta2)*rmvnorm(m2,rep(0,d),diag(d)) + theta2*rmvnorm(m2,mu2,diag(d))
x = rbind(x1,x2)
    
y1 = rmvnorm(m1,rep(0,d),diag(d))
y2 = rmvnorm(m2,rep(0,d),diag(d))
y = rbind(y1,y2)

clfdr1 = 0.8*dmvnorm(x1,rep(0,d),diag(d))/mkde(x1)
clfdr2 = 0.9*dmvnorm(x2,rep(0,d),diag(d))/mkde(x2)

clfdr = c(clfdr1,clfdr2)
rcl = sort(clfdr)
hk=1
while(hk<m && (1/hk)*sum(rcl[1:hk])<al){
  hk=hk+1
  }
hk = hk-1
rej = (clfdr<=rcl[hk])
subfdr.1[l] = sum(rej*(1-theta)) / max(1, sum(rej) )
subtdp.1[l] = sum(rej*theta) / max(1, sum(theta) )

#new mwthod
data1 = rbind(x1,y1)
g1ds = (1+0.8)/2 * dmvnorm(data1,rep(0,d),diag(d))/mkde(data1)
g1ds = pmin(g1ds,0.999)
g1ds = (0.5-0.1)/(1-0.1)*g1ds/(1-g1ds)
s1 = g1ds[1:m1]
st1 = g1ds[-(1:m1)]

data2 = rbind(x2,y2)
g2ds = (1+0.9)/2 * dmvnorm(data2,rep(0,d),diag(d))/mkde(data2)
g2ds = pmin(g2ds,0.999)
g2ds = (0.5-0.05)/(1-0.05)*g2ds/(1-g2ds)
s2 = g2ds[1:m2]
st2 = g2ds[-(1:m2)]

s = c(s1,s2)
st = c(st1,st2)

#new mwthod
    st = st/(1+st)
    s = s/(1+s)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
}

Dimension = rep("d4=16",400)
METHOD = c(rep("Clfdr",200), rep("CLAW",200))
FDR = c(subfdr.1,subfdr.new)
AP = c(subtdp.1,subtdp.new)
add.fdr = data.frame(Dimension, METHOD, FDR)
add.ap = data.frame(Dimension, METHOD, AP)

data.box.fdr = rbind(data.box.fdr,add.fdr)
data.box.ap = rbind(data.box.ap,add.ap)


## ----warning=FALSE-------------------------------------------------------------------------
al = 0.05
d = 24
pp01 = 0.8
pp02 = 0.9

m1=1000
m2=2500
m=m1+m2

mu1 = c(rep(sqrt(2*log(d)),4), rep(0,d-4))
mu2 = rep(2,d)

subfdr.1 = NA
subtdp.1 = NA
subfdr.new = NA
subtdp.new = NA

for (l in 1:200) {
  
theta1 = rbinom(m1,1,0.2)
theta2 = rbinom(m2,1,0.1)
theta = c(theta1,theta2)

x1 = (1-theta1)*rmvnorm(m1,rep(0,d),diag(d)) + theta1*rmvnorm(m1,mu1,diag(d))
x2 = (1-theta2)*rmvnorm(m2,rep(0,d),diag(d)) + theta2*rmvnorm(m2,mu2,diag(d))
x = rbind(x1,x2)
    
y1 = rmvnorm(m1,rep(0,d),diag(d))
y2 = rmvnorm(m2,rep(0,d),diag(d))
y = rbind(y1,y2)

clfdr1 = 0.8*dmvnorm(x1,rep(0,d),diag(d))/mkde(x1)
clfdr2 = 0.9*dmvnorm(x2,rep(0,d),diag(d))/mkde(x2)

clfdr = c(clfdr1,clfdr2)
rcl = sort(clfdr)
hk=1
while(hk<m && (1/hk)*sum(rcl[1:hk])<al){
  hk=hk+1
  }
hk = hk-1
rej = (clfdr<=rcl[hk])
subfdr.1[l] = sum(rej*(1-theta)) / max(1, sum(rej) )
subtdp.1[l] = sum(rej*theta) / max(1, sum(theta) )

#new mwthod
data1 = rbind(x1,y1)
g1ds = (1+0.8)/2 * dmvnorm(data1,rep(0,d),diag(d))/mkde(data1)
g1ds = pmin(g1ds,0.99999)
g1ds = (0.5-0.1)/(1-0.1)*g1ds/(1-g1ds)
s1 = g1ds[1:m1]
st1 = g1ds[-(1:m1)]

data2 = rbind(x2,y2)
g2ds = (1+0.9)/2 * dmvnorm(data2,rep(0,d),diag(d))/mkde(data2)
g2ds = pmin(g2ds,0.99999)
g2ds = (0.5-0.05)/(1-0.05)*g2ds/(1-g2ds)
s2 = g2ds[1:m2]
st2 = g2ds[-(1:m2)]

s = c(s1,s2)
st = c(st1,st2)

#new mwthod
    st = st/(1+st)
    s = s/(1+s)
    #new mwthod
    de.new = conf.q(s,st,al)
    subfdr.new[l] = sum(de.new*(1-theta)) / max(1, sum(de.new) )
    subtdp.new[l] = sum(de.new*theta) / max(1, sum(theta) )
}

Dimension = rep("d5=24",400)
METHOD = c(rep("Clfdr",200), rep("CLAW",200))
FDR = c(subfdr.1,subfdr.new)
AP = c(subtdp.1,subtdp.new)
add.fdr = data.frame(Dimension, METHOD, FDR)
add.ap = data.frame(Dimension, METHOD, AP)

data.box.fdr = rbind(data.box.fdr,add.fdr)
data.box.ap = rbind(data.box.ap,add.ap)


## ----message=FALSE, warning=FALSE----------------------------------------------------------
box.fdr = ggplot(data.box.fdr,aes(x=Dimension, y=FDR, fill=METHOD))+
  geom_hline(aes(yintercept=al),size=1)+geom_boxplot(width=0.5,outlier.shape = NA)+
  ylim(0,0.3)+theme_bw()
print(box.fdr)


## ------------------------------------------------------------------------------------------
box.ap = ggplot(data.box.ap,aes(x=Dimension, y=AP, fill=METHOD))+
  geom_boxplot(width=0.5,outlier.shape = NA)+theme_bw()
print(box.ap)


## ----message=FALSE, warning=FALSE----------------------------------------------------------
grid_arrange_shared_legend(box.fdr,box.ap,nrow=1,ncol=2,position = "top")


