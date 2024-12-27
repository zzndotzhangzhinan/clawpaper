lin.itp<-function(x, X, Y){
  ## x: the coordinates of points where the density needs to be interpolated
  ## X: the coordinates of the estimated densities
  ## Y: the values of the estimated densities
  ## the output is the interpolated densities
  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else 
      y[k]<-Y[i]
  }
  return(y)
}


exscore.func = function(data1, data2, f0){
  m = length(data1)
  mix.dens = density(c(data1,data2),from=min(c(data1,data2))-10, 
                     to=max(c(data1,data2))+10, n=1000)
  drx = dnorm(data1,f0[1],f0[2])/ lin.itp(data1, mix.dens$x, mix.dens$y)
  drm = dnorm(data2,f0[1],f0[2])/ lin.itp(data2, mix.dens$x, mix.dens$y)
  s = drx
  st = drm
  
  return(list(s=s,st=st))
}

Adadetect = function(x,y,al){
  m = length(x)
  t = sort(unique(c(x, y)))
  fr = NA
  for (j in 1:length(t)){
    fr[j] = (sum(y<=t[j])+1)/(m+1) /(sum(x<=t[j])/m)
  }
  t.hat = max(na.omit(t[fr<=al]))
  decision = (x<=t.hat)
  return(decision)
}


cor_1D.func = function(x, y, tau=0.5, h){
  ## pis_est.func calculates the conditional proportions pis
  ## Arguments
  # x: z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: the conditional proportions
  pval = 2*pnorm(-abs(x),0,1)
  pval2 = 2*pnorm(-abs(y),0,1)
  m <- length(pval)
  s <- 1:m # auxiliary variable
  p.est <-rep(0, m)
  
  mix.dens = density(c(x,y),from=min(c(x,y))-10, 
                     to=max(c(x,y))+10, n=1000)
  bw = mix.dens$bw
  dx = rep(NA,m)
  dy = rep(NA,m)
  
  for (i in 1:m) { 
    kht = dnorm(s-i, 0, h)
    p.est[i] = (sum(kht[which(pval>=tau)])+ sum(kht[which(pval2>=tau)]))/(2*(1-tau)*sum(kht))
    dx[i] = sum(kht * (dnorm(x[i]-x,0,bw)+dnorm(x[i]-y,0,bw)) ) /(2*sum(kht))
    dy[i] = sum(kht * (dnorm(y[i]-x,0,bw)+dnorm(y[i]-y,0,bw)) ) /(2*sum(kht))
  }
  p.est[which(p.est>1)] = 1
  
  
  drx = p.est*dnorm(x,0,1)/dx
  dry = p.est*dnorm(y,0,1)/dy
  
  drx = pmin(drx,0.999)
  dry = pmin(dry,0.999)
  
  #drx = drx/(1-drx)*(p.est-0.5)/p.est
  #dry = dry/(1-dry)*(p.est-0.5)/p.est
  drx = drx/(1+drx)
  dry = dry/(1+dry)
  
  return(list(drx=drx,dry=dry))
}


conf.q = function(x,y,al){
  m = length(x)
  Wt = (2*as.numeric(y>=x)-1) * pmax((1-y),(1-x))* (y!=x)
  t = sort(unique( abs(Wt) ))
  fr = NA
  for (j in 1:length(t)){
    fr[j] = (sum(Wt<=(-t[j]))+1) /(sum(Wt>=t[j]))
  }
  t.hat = min(na.omit(t[fr<=al]))
  decision = (Wt>=t.hat)
  return(decision)
}


bh = function(pv,al){
  #pv = 2*(1-pnorm(abs(x),0,1))
  adpv = p.adjust( pv ,method = "BH")
  decision = (pv <= max( pv[adpv<=al] ))
  return(decision)
}


law.func<-function(pvs, pis, q){
  ## implementing "spatial multiple testing by locally adaptive weighting"
  ## Arguments
  # pvs: p-values
  # pis: conditional probabilities
  # q: FDR level
  ## Values
  # de: the decision
  # th: the threshold for weighted p-values
  m<-length(pvs)
  nu<-10e-5
  pis[which(pis<nu)]<-nu # stabilization
  pis[which(pis>1-nu)]<-1-nu # stabilization
  ws<-pis/(1-pis)
  pws<-pvs/ws
  st.pws<-sort(pws)
  fdps<-sum(pis)*st.pws/(1:m)
  de<-rep(0, m)
  if(sum(fdps<=q)==0)
  {
    k<-0
    pwk<-1
  }
  else
  {
    k<-max(which(fdps<=q))
    pwk<-st.pws[k]
    de[which(pws<=pwk)]<-1
  }
  #y<-list(nr=k, th=pwk, de=de)
  y = de
  return (y)
}


sab.func<-function(pvs, pis, q){
  ## implementing "SABHA" by Li and Barber
  ## Arguments
  # pvs: p-values
  # pis: conditional probabilities
  # q: FDR level
  ## Values
  # de: the decision
  # th: the threshold for weighted p-values
  m<-length(pvs)
  nu<-10e-5
  pis[which(pis>1-nu)]<-1-nu # stabilization
  pws<-pvs*(1-pis)
  st.pws<-sort(pws)
  
  pwi<-st.pws/1:m
  de<-rep(0, m)
  if (sum(pwi<=q/m)==0)
  {
    k<-0
    pk<-1
  }
  else
  {
    k<-max(which(pwi<=(q/m)))
    pk<-st.pws[k]
    de[which(pws<=pk)]<-1
  }
  #y<-list(nr=k, th=pk, de=de)
  y = de
  return (y)
}	 


####################################################
## ## ## weights estimation for LAWS and SABHA

pis_1D.func<- function(pval, tau=0.5, h){
  ## pis_est.func calculates the conditional proportions pis
  ## Arguments
  # x: z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth h=50 default
  ## Values
  # pis: the conditional proportions
  #pval <- 2*pnorm(-abs(x))
  m <- length(pval)
  s <- 1:m # auxiliary variable
  p.est <-rep(0, m)
  for (i in 1:m) { 
    kht<-dnorm(s-i, 0, h)
    p.est[i]<-sum(kht[which(pval>=tau)])/((1-tau)*sum(kht))
  }
  p.est[which(p.est>1)] <-1
  return(1-p.est)
}


####################################################
## ## ## weights estimation for 2D case
disvec.func<-function(dims, s)
{
  # disvec computes the distances of all points on a m1 times m2 spatial domain to a point s
  ## Arguments:
  # dims=c(d1, d2): the dimensions
  # s=c(s1, s2): a spatial point
  ## Values:
  # a vector of distances
  
  m<-dims[1]*dims[2]
  dis.vec<-rep(0, m)
  for(i in 1:dims[1])
  {
    dis.vec[((i-1)*dims[2]+1):(i*dims[2])]<-sqrt((i-s[1])^2+(1:dims[2]-s[2])^2)
  }
  return(dis.vec) 
}

pis_2D.func<- function(x, tau=0.5, h)
{
  ## pis_2D.func calculates the conditional proportions pis
  ## Arguments
  # x: a matrix of z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: conditional proportions
  
  dims<-dim(x)
  m<-dims[1]*dims[2]
  x.vec<-c(t(x))
  pv.vec<-2*pnorm(-abs(x.vec))
  scr.idx<-which(pv.vec>=tau)
  p.est<-matrix(rep(0, m), dims[1], dims[2])  
  
  for (i in 1:dims[1]) 
  {
    for (j in 1:dims[2]) 
    {
      s<-c(i, j)
      dis.vec<-disvec.func(dims, s)
      kht<-dnorm(dis.vec, 0, h)
      p.est[i,j]<-min(1-1e-5, sum(kht[scr.idx])/((1-tau)*sum(kht)))
    }
  }
  pis.est<-1-c(p.est)
  return(pis.est)
}


cor_2D.func = function(x, y, tau=0.5, h){
  ## pis_est.func calculates the conditional proportions pis
  ## Arguments
  # x: z-values
  # tau: the screening threshold, which can be prespecified or chosen adaptively
  # bdw: bandwidth
  ## Values
  # pis: the conditional proportions
  
  dims<-dim(x)
  m<-dims[1]*dims[2]
  x.vec = c(x)
  y.vec = c(y)
  pv.vec = 2*pnorm(-abs(x.vec))
  pv.y = 2*pnorm(-abs(y.vec))
  
  scr.idx<-which(pv.vec>=tau)
  scr.idy<-which(pv.y>=tau)
  
  p.est<-matrix(rep(0, m), dims[1], dims[2])
  dx<-matrix(rep(0, m), dims[1], dims[2])
  dy<-matrix(rep(0, m), dims[1], dims[2])
  
  mix.dens = density(c(x.vec,y.vec),from=min(c(x.vec,y.vec))-10, 
                     to=max(c(x.vec,y.vec))+10, n=1000)
  bw = mix.dens$bw
  
  
  for (i in 1:dims[1]) 
  {
    for (j in 1:dims[2]) 
    {
      s<-c(i, j)
      dis.vec<-disvec.func(dims, s)
      kht<-dnorm(dis.vec, 0, h)
      p.est[i,j]<-min(1-1e-5, (sum(kht[scr.idx])+sum(kht[scr.idy]))/(2*(1-tau)*sum(kht)))
      dx[i,j] = sum(kht * (dnorm(x[i,j]-x.vec,0,bw)+dnorm(x[i,j]-y.vec,0,bw)) ) /(2*sum(kht))
      dy[i,j] = sum(kht * (dnorm(y[i,j]-x.vec,0,bw)+dnorm(y[i,j]-y.vec,0,bw)) ) /(2*sum(kht))
    }
  }
  
  dx = c(dx)
  dy = c(dy)
  p.est = c(p.est)
  
  
  drx = p.est*dnorm(x.vec,0,1)/dx
  dry = p.est*dnorm(y.vec,0,1)/dy
  
  drx = pmin(drx,0.999)
  dry = pmin(dry,0.999)
  
  drx = drx/(1-drx)*(p.est-0.5)/p.est
  dry = dry/(1-dry)*(p.est-0.5)/p.est
  
  drx = drx/(1+drx)
  dry = dry/(1+dry)

  return(list(drx=drx,dry=dry))
}



