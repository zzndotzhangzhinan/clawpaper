library(IHWpaper)

source("coordinate_func2.R")
library(ggplot2)
library(reshape2)
library(lemon)
library(doParallel)
library(foreach)
library(adaptMT)
library(splines)
library(IHW)
library(Compositional)


data = read.csv("yeast.csv")
ss = data$covariate
pv = data$pvalue
m = length(pv)
pv = pmin(pv,0.999999)
pv = pmax(pv,0.000001)

x = qnorm(1-pv/2)#qnorm(1-pv/2)#qnorm(1-pv/2)#(1-2*rbinom(m,1,0.5))*qnorm(1-pv/2)
y = rnorm(m,0,1)

al = 0.05
f0 = c(0,1)



#pv = 2*(1-pnorm(abs(x),0,1))
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
sum(de.new)

#BH
#pv = 2*(1-pnorm(abs(x),0,1))

de.bh = bh(pv,al)
sum(de.bh)

#Adadetect
score.pool = exscore.func(x,y,f0)
pools = score.pool$s
poolst = score.pool$st
de.ad = Adadetect(pools,poolst,al)
sum(de.ad)

#LAWS
#pis = pis_1D.func(pv,tau=0.5, log(ss))#pis_1D.laws(pv,tau=0.5, ss)

de.law = law.func(pv, pis, al)
sum(de.law)


#SABHA
de.sab = sab.func(pv, pis, al)
sum(de.sab)

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
sum(adapt.de)

#IHW
ihw_fdr = ihw(pv, covariates=ss, alpha=al, nbins=5)
ihw.de = (adj_pvalues(ihw_fdr)<=0.05)
sum(ihw.de)


alvec = c(0.045,0.05,0.055,0.06)

tdp.new = matrix(NA,length(alvec),100)
tdp.bh = NA
tdp.ad = matrix(NA,length(alvec),100)
tdp.law = NA
tdp.sab = NA
tdp.adapt = NA
tdp.ihw = NA

for (k in 1:length(alvec)) {
  al = alvec[k]
  
  #BH
  de.bh = bh(pv,al)
  tdp.bh[k] = sum(de.bh)
  
  #LAWS
  de.law = law.func(pv, pis, al)
  tdp.law[k] = sum(de.law)
  
  #SABHA
  de.sab = sab.func(pv, pis, al)
  tdp.sab[k] = sum(de.sab)
  
  #AdaPT
  adapt.de = (adapt.res$qvals<=al)
  tdp.adapt[k] = sum(adapt.de)
  
  #IHW
  ihw.de = (adj_pvalues(ihw_fdr)<=al)
  tdp.ihw[k] = sum(ihw.de)
}


set.seed(123456)

for (l in 1:100) {
  y = rnorm(m,0,1)
  
  pvy = 2*(1-pnorm(abs(y),0,1))
  pisy = pis_1D.func(pvy,tau=0.5,ss)
  piclaw = (pis+pisy)/2
  pdf1 = dnorm(c(x,y),0,1)
  sdensity = density( ss,from=min(ss)-10, to=max(ss)+10, n=1000)
  pdf2 = lin.itp(ss, sdensity$x, sdensity$y)
  mat = matrix(c(x,y,ss,ss), ncol=2, byrow = F)
  bidens = mkde(mat)
  s = (1-piclaw)*pdf1[1:m]*pdf2/bidens[1:m]
  st = (1-piclaw)*pdf1[-(1:m)]*pdf2/bidens[-(1:m)]
  
  
  #Adadetect
  score.pool = exscore.func(x,y,f0)
  pools = score.pool$s
  poolst = score.pool$st
  
  for (k in 1:length(alvec)) {
    al = alvec[k]
    
    de.new = conf.q(s,st,al)
    tdp.new[k,l] =  sum(de.new)
    
    de.ad = Adadetect(pools,poolst,al)
    tdp.ad[k,l] = sum(de.ad)
  }
}


ex1 = data.frame(alvec=c("0.045","0.05","0.055","0.06"), 
                 CLAW=as.vector(apply(tdp.new, 1, median)), BH=tdp.bh, #LAWS=tdp.law,
                 AdaDetect=as.vector(apply(tdp.ad, 1, median)), LAWS=tdp.law,
                 SABHA=tdp.sab, AdaPT=tdp.adapt, IHW=tdp.ihw)
ex1 = melt(ex1,id="alvec")
colnames(ex1) = c("alvec","METHOD","value")

line1 = ggplot(ex1, aes(x=alvec, y=value))+
  geom_bar(stat = "identity", aes(fill=METHOD), position = position_dodge(0.6),
           width=0.5,color="black")+ labs(x="Nominal FDR level",y="Number of discoveries")+
  theme_bw()

print(line1)




box1tdp= NULL
for (i in 1:100) {
  
  rawtdp.df = data.frame(alvec, CLAW=tdp.new[,i], AdaDetect=tdp.ad[,i])
  rawtdp.df = melt(rawtdp.df,id="alvec")
  colnames(rawtdp.df) = c("alvec","METHOD","value")
  rawtdp.df = cbind(PID=rep(i,8),rawtdp.df)
  box1tdp = rbind(box1tdp, rawtdp.df)
}

boxplot1 = ggplot() + geom_boxplot(data=box1tdp,mapping = aes(x=as.factor(alvec),y=value,fill=METHOD), width=0.3, position = position_dodge(width=0.4)) + 
  labs(x="Nominal FDR level",y="Number of discoveries") +
  theme_bw()

print(boxplot1)








