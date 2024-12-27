
bwfw1.hmm1<-function(x, pii, A, f0, f1)#function to calculate $s^{X}$
{
  
  ###################################################################
  
  ## USAGE
  # bwfw.hmm(x, pii, A, f0, f1)
  
  ## ARGUMENTS
  # x=(x[1], ..., x[m]): the observed data
  # pii=(pii[0], pii[1]): the initial state distribution
  # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
  # f0=(mu, sigma): the parameters for null distribution
  # f1=(mu, sigma): the non-null distribution
  
  ## DETAILS
  # bwfw.hmm calculates values for backward, forward variables, probabilities of hidden states, 
  # --the lfdr variables and etc. 
  # --using the forward-backward procedure (Baum et al.) 
  # --based on a sequence of observations for a given hidden markov model M=(pii, A, f)
  # the underflow problem was fixed by using the rescaled forward and backward variables
  # -- see Stamp (2004) for a detailed instruction on the coding of this algorithm
  
  ## VALUES
  # alpha: rescaled backward variables
  # beta: rescaled forward variables
  # lfdr: lfdr variables
  # gamma: probabilities of hidden states
  # dgamma: rescaled transition variables
  # omega: rescaled weight variables
  
  ################################################################## 
  
  ## Initialize
  
  NUM<-length(x)
  
  ## Densities
  
  f0x<-dnorm(x, f0[1], f0[2])
  f1x<-dnorm(x, f1[1], f1[2])
  
  ## the backward-forward procedure
  
  # a. the backward variables
  # --rescaled 
  
  alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=T)
  # scaling variable c_0
  c0<-rep(0, NUM)
  
  alpha[1, 1]<-pii[1]*f0x[1]
  alpha[1, 2]<-pii[2]*f1x[1]
  # rescaling alpha
  c0[1]<-1/sum(alpha[1, ])
  alpha[1, ]<-c0[1]*alpha[1, ]
  
  for (k in 1:(NUM-1))
  { 
    alpha[k+1, 1]<-(alpha[k, 1]*A[1, 1]+alpha[k, 2]*A[2, 1])*f0x[k+1]
    alpha[k+1, 2]<-(alpha[k, 1]*A[1, 2]+alpha[k, 2]*A[2, 2])*f1x[k+1]
    # rescaling alpha
    c0[k+1]<-1/sum(alpha[k+1, ])
    alpha[k+1, ]<-c0[k+1]*alpha[k+1, ]
  }
  
  # b. the forward variables
  # --rescaled
  
  beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=T)
  
  beta[NUM, 1]<-c0[NUM]
  beta[NUM, 2]<-c0[NUM]
  
  for (k in (NUM-1):1)
  { 
    beta[k, 1]<-A[1, 1]*f0x[k+1]*beta[k+1, 1]+A[1, 2]*f1x[k+1]*beta[k+1, 2]
    beta[k, 2]<-A[2, 1]*f0x[k+1]*beta[k+1, 1]+A[2, 2]*f1x[k+1]*beta[k+1, 2]
    # rescaling beta
    # using the same scaling factors as alpha 
    beta[k, ]<-c0[k]*beta[k, ]
  }
  
  # c. lfdr variables
  # --original
  # --the same formulae hold for the rescaled alpha and beta
  
  lfdr<-rep(0, NUM)
  
  for (k in 1:NUM)
  { 
    q1<-alpha[k, 1]*beta[k, 1]
    q2<-alpha[k, 2]*beta[k, 2]
    lfdr[k]<-q1/(q1+q2)
  }
  
  bwfw.var<-list(bw=alpha, fw=beta, lf=lfdr)
  return(bwfw.var)
  
}



bwfw1.hmm.ad<-function(x, pii, A, f0, f1, indx)#function to calculate $s_{i}^{Y}$
{
  
  NUM<-length(x)
  
  ## Densities
  
  f0x<-dnorm(x, f0[1], f0[2])
  f1x<-dnorm(x, f1[1], f1[2])
  
  ## the backward-forward procedure
  
  # a. the backward variables
  # --rescaled 
  
  alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=T)
  # scaling variable c_0
  c0<-rep(0, NUM)
  
  alpha[1, 1]<-pii[1]*f0x[1]
  alpha[1, 2]<-pii[2]*f1x[1]
  # rescaling alpha
  c0[1]<-1/sum(alpha[1, ])
  alpha[1, ]<-c0[1]*alpha[1, ]
  
  for (k in 1:(NUM-1))
  { 
    alpha[k+1, 1]<-(alpha[k, 1]*A[1, 1]+alpha[k, 2]*A[2, 1])*f0x[k+1]
    alpha[k+1, 2]<-(alpha[k, 1]*A[1, 2]+alpha[k, 2]*A[2, 2])*f1x[k+1]
    # rescaling alpha
    c0[k+1]<-1/sum(alpha[k+1, ])
    alpha[k+1, ]<-c0[k+1]*alpha[k+1, ]
  }
  
  # b. the forward variables
  # --rescaled
  
  beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=T)
  
  beta[NUM, 1]<-c0[NUM]
  beta[NUM, 2]<-c0[NUM]
  
  for (k in (NUM-1):1)
  { 
    beta[k, 1]<-A[1, 1]*f0x[k+1]*beta[k+1, 1]+A[1, 2]*f1x[k+1]*beta[k+1, 2]
    beta[k, 2]<-A[2, 1]*f0x[k+1]*beta[k+1, 1]+A[2, 2]*f1x[k+1]*beta[k+1, 2]
    # rescaling beta
    # using the same scaling factors as alpha 
    beta[k, ]<-c0[k]*beta[k, ]
  }
  
  # c. lfdr variables
  # --original
  # --the same formulae hold for the rescaled alpha and beta
  
  lfdr<- NA
  
  q1<-alpha[indx, 1]*beta[indx, 1]
  q2<-alpha[indx, 2]*beta[indx, 2]
  lfdr<-q1/(q1+q2)
  
  
  return(lfdr)
}



rdata1.hmm<-function(NUM, pii, A, f0, f1)#function to generate HMM observations
{
  theta<-rep(0, NUM)
  x<-rep(0, NUM)
  
  ## generating the states
  # initial state
  theta[1]<-rbinom(1, 1, pii[2])
  # other states
  for (i in 2:NUM)
  {
    if (theta[i-1]==0)
      theta[i]<-rbinom(1, 1, A[1, 2])
    else
      theta[i]<-rbinom(1, 1, A[2, 2])
  }
  
  ## generating the observations
  for (i in 1:NUM)
  {
    if (theta[i]==0)
    {
      x[i]<-rnorm(1, mean=f0[1], sd=f0[2])
    }
    else
    { 
      x[i]<-rnorm(1, mean=f1[1], sd=f1[2])
    }
  }
  data<-list(s=theta, o=x)
  return (data)
}




bwfw1.hmm<-function(x, pii, A, f0, f1)#function to estimate HMM
{
  
  ###################################################################
  
  ## USAGE
  # bwfw.hmm(x, pii, A, f0, f1)
  
  ## ARGUMENTS
  # x=(x[1], ..., x[m]): the observed data
  # pii=(pii[0], pii[1]): the initial state distribution
  # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
  # f0=(mu, sigma): the parameters for null distribution
  # f1=(mu, sigma): the non-null distribution
  
  ## DETAILS
  # bwfw.hmm calculates values for backward, forward variables, probabilities of hidden states, 
  # --the lfdr variables and etc. 
  # --using the forward-backward procedure (Baum et al.) 
  # --based on a sequence of observations for a given hidden markov model M=(pii, A, f)
  # the underflow problem was fixed by using the rescaled forward and backward variables
  # -- see Stamp (2004) for a detailed instruction on the coding of this algorithm
  
  ## VALUES
  # alpha: rescaled backward variables
  # beta: rescaled forward variables
  # lfdr: lfdr variables
  # gamma: probabilities of hidden states
  # dgamma: rescaled transition variables
  # omega: rescaled weight variables
  
  ################################################################## 
  
  ## Initialize
  
  NUM<-length(x)
  
  ## Densities
  
  f0x<-dnorm(x, f0[1], f0[2])
  f1x<-dnorm(x, f1[1], f1[2])
  
  ## the backward-forward procedure
  
  # a. the backward variables
  # --rescaled 
  
  alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=T)
  # scaling variable c_0
  c0<-rep(0, NUM)
  
  alpha[1, 1]<-pii[1]*f0x[1]
  alpha[1, 2]<-pii[2]*f1x[1]
  # rescaling alpha
  c0[1]<-1/sum(alpha[1, ])
  alpha[1, ]<-c0[1]*alpha[1, ]
  
  for (k in 1:(NUM-1))
  { 
    alpha[k+1, 1]<-(alpha[k, 1]*A[1, 1]+alpha[k, 2]*A[2, 1])*f0x[k+1]
    alpha[k+1, 2]<-(alpha[k, 1]*A[1, 2]+alpha[k, 2]*A[2, 2])*f1x[k+1]
    # rescaling alpha
    c0[k+1]<-1/sum(alpha[k+1, ])
    alpha[k+1, ]<-c0[k+1]*alpha[k+1, ]
  }
  
  # b. the forward variables
  # --rescaled
  
  beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=T)
  
  beta[NUM, 1]<-c0[NUM]
  beta[NUM, 2]<-c0[NUM]
  
  for (k in (NUM-1):1)
  { 
    beta[k, 1]<-A[1, 1]*f0x[k+1]*beta[k+1, 1]+A[1, 2]*f1x[k+1]*beta[k+1, 2]
    beta[k, 2]<-A[2, 1]*f0x[k+1]*beta[k+1, 1]+A[2, 2]*f1x[k+1]*beta[k+1, 2]
    # rescaling beta
    # using the same scaling factors as alpha 
    beta[k, ]<-c0[k]*beta[k, ]
  }
  
  # c. lfdr variables
  # --original
  # --the same formulae hold for the rescaled alpha and beta
  
  lfdr<-rep(0, NUM)
  
  for (k in 1:NUM)
  { 
    q1<-alpha[k, 1]*beta[k, 1]
    q2<-alpha[k, 2]*beta[k, 2]
    lfdr[k]<-q1/(q1+q2)
  }
  
  # d. probabilities of hidden states
  # -- and transition variables
  # -- both are rescaled
  
  gamma<-matrix(1:(NUM*2), NUM, 2, byrow=T)
  # initialize gamma[NUM]
  gamma[NUM, ]<-c(lfdr[NUM], 1-lfdr[NUM])
  dgamma<-array(rep(0, (NUM-1)*4), c(2, 2, (NUM-1)))
  
  for (k in 1:(NUM-1))
  {
    denom<-0
    for (i in 0:1)
    {
      for (j in 0:1)
      { 
        fx<-(1-j)*f0x[k+1]+j*f1x[k+1]
        denom<-denom+alpha[k, i+1]*A[i+1, j+1]*fx*beta[k+1, j+1]
      }
    }
    for (i in 0:1)
    {
      gamma[k, i+1]<-0
      for (j in 0:1)
      { 
        fx<-(1-j)*f0x[k+1]+j*f1x[k+1]
        dgamma[i+1, j+1, k]<-alpha[k, i+1]*A[i+1, j+1]*fx*beta[k+1, j+1]/denom
        gamma[k, i+1]<-gamma[k, i+1]+dgamma[i+1, j+1, k]  
      }
    }
  }
  
  # f. return the results of the bwfw proc.
  
  bwfw.var<-list(bw=alpha, fw=beta, lf=lfdr, pr=gamma, ts=dgamma)
  return(bwfw.var)
  
}
em1.hmm<-function(x, maxiter=200)
{
  
  #####################################################################################
  
  ## USAGE
  # em.hmm(x, L)
  
  ## ARGUMENTS
  # x: the observed data
  # L: num. of components in the non-null mixture
  # maxiter: the maximum number of iterations
  
  ## DETAILS
  # em.hmm calculates the MLE for a HMM model with hidden states being 0/1.
  # the distribution of state 0 is assumed to be normal
  # the distribution of state 1 is assumed to be a normal mixture with L components
  
  ## VALUES
  # fuction 'em.hmm' gives the MLE of model parameters and Lfdr estimates for a HMM 
  # pii: the initial state distribution
  # A=(a00 a01\\ a10 a11): transition matrix
  # f0: the null distribution
  # pc, i from 1 to L: 
  # --probability weights of each component in the non-null mixture 
  # f1: an L by 2 matrix
  # --specifying the dist. of each component in the non-null mixture 
  # lfdr: the lfdr variables
  # niter: number of iterations
  
  ####################################################################################
  
  NUM<-length(x)
  # precision tolerance level
  ptol<-1e-4
  niter<-0
  
  ### initializing model parameters
  
  pii.new<-c(1, 0)
  A.new<-matrix(c(0.8, 0.2, 0.4, 0.6), 2, 2, byrow=T)
  f0<-c(0, 1)
  f1.new<-c(2, 1)
  
  diff<-1
  
  ### The E-M Algorithm
  
  while(diff>ptol && niter<maxiter)
  {
    
    niter<-niter+1
    
    pii.old<-pii.new
    A.old<-A.new
    f1.old<-f1.new
    
    ## updating the weights and probabilities of hidden states
    
    bwfw.res<-bwfw1.hmm(x, pii.old, A.old, f0, f1.old)
    
    # the hidden states probabilities
    gamma<-bwfw.res$pr
    # the transition variables
    dgamma<-bwfw.res$ts
    
    ## updating the parameter estimates
    
    # a. initial state distribution
    
    for (i in 0:1)
    {
      pii.new[i+1]<-gamma[1, i+1]
    }
    
    # b. transition matrix
    
    for (i in 0:1)
    {
      for (j in 0:1)
      { 
        q1<-sum(dgamma[i+1, j+1, ])
        q2<-sum(gamma[1:(NUM-1), i+1])
        A.new[i+1, j+1]<-q1/q2  
      }
    }
    
    # c. non-null distribution 
    
    q1<-sum(gamma[, 2])
    q2<-sum(gamma[, 2]*x)
    mu1<-q2/q1
    q3<-sum(gamma[, 2]*(x-mu1)*(x-mu1))
    sd1<-sqrt(q3/q1)
    f1.new<-c(mu1, sd1)
    
    df1<-abs(A.old-A.new)
    df2<-abs(f1.old-f1.new)
    diff<-max(df1, df2)
    
  }
  
  # f. the final local fdr statistic
  lfdr<-gamma[, 1]
  
  # g. return the results of the E-M algorithm
  
  em.var<-list(pii=pii.new, A=A.new, f1=f1.new, lf=lfdr, ni=niter)
  return (em.var)
  
}



#linear interpolation
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



Func2 = function(j, Z, X, Y)#function to calculate PLIS statistics for x and m
  { 
  xtil = Z
  xtil[j] = Y[j]
  #lism[j] = bwfw1.hmm.ad(xtil, pii, B, f0, f1, j)
  sublism = bwfw1.hmm.ad(xtil, pii, estA, f0, estf1, j)
  xtil[j] = X[j]
  #lisx[j] = bwfw1.hmm.ad(xtil, pii, B, f0, f1, j)
  sublisx = bwfw1.hmm.ad(xtil, pii, estA, f0, estf1, j)
  return(c(sublisx,sublism))
}

