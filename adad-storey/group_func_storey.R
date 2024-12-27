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


group.func = function(data1, data2, f0){
  m = length(data1)
  mix.dens = density(c(data1,data2),from=min(c(data1,data2))-10, 
                     to=max(c(data1,data2))+10, n=1000)
  drx = dnorm(data1,f0[1],f0[2])/ lin.itp(data1, mix.dens$x, mix.dens$y)
  drm = dnorm(data2,f0[1],f0[2])/ lin.itp(data2, mix.dens$x, mix.dens$y)
  
  pvs = 2*(1-pnorm(abs(c(data1,data2)),f0[1],f0[2]))
  pp0 = sum(pvs>0.5)/(0.5*2*m)
  if (pp0==Inf | is.na(pp0)){
    pp0 = 1
  }
  
  s = pp0 * drx
  st = pp0 * drm
  
  s = pmin(s,0.999)
  st = pmin(st,0.999)
  pp = 1-pp0
  
  s = (0.5-pp)/pp0*s/(1-s)
  st = (0.5-pp)/pp0*st/(1-st)
  s = s/(1+s)
  st = st/(1+st)
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

conf.q = function(x,y,al){
  m = length(x)
  Wt = (2*as.numeric(y>=x)-1) * pmax((1-y),(1-x))
  t = sort(unique( abs(Wt) ))
  fr = NA
  for (j in 1:length(t)){
    fr[j] = (sum(Wt<=(-t[j]))+1) /(sum(Wt>=t[j]))
  }
  t.hat = min(na.omit(t[fr<=al]))
  decision = (Wt>=t.hat)
  return(decision)
}

bh = function(x,al){
  pv = 2*(1-pnorm(abs(x),0,1))
  adpv = p.adjust( pv ,method = "BH")
  decision = (pv <= max( pv[adpv<=al] ))
  return(decision)
}


build.conf.p = function(sx,sy){
  m = length(sx)
  n = length(sy)
  pv = NA
  for (i in 1:m) {
    pv[i] = (1+sum(sx[i]>sy))/(n+1)
  }
  return(pv)
}


pi0.storey = function(pv,lam){
  pi0 = (1+sum(pv>lam))/(length(pv)*(1-lam))
  return(pi0)
}


