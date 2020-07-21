#Simon Wood's answer to stack overflow question
#https://stats.stackexchange.com/questions/464961/issues-with-posterior-simulation-from-gam


## some useful densities (require mgcv::rmvn)...

rmvt <- function(n,mu,V,df) {
  ## simulate multivariate t variates  
  y <- rmvn(n,mu*0,V)
  v <- rchisq(n,df=df)
  t(mu + t(sqrt(df/v)*y))
}

dmvt <- function(x,mu,V,df,R=NULL) {
  ## multivariate t log density...
  p <- length(mu);
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R),x-mu)
  k <- - sum(log(diag(R))) - p*log(df*pi)/2 + lgamma((df+p)/2) - lgamma(df/2)
  k - if (is.matrix(z)) (df+p)*log1p(colSums(z^2)/df)/2 else (df+p)*log1p(sum(z^2)/df)/2
}

dmvn <- function(x,mu,V,R=NULL) {
  ## multivariate normal density mgcv:::rmvn can be used for generation 
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R),x-mu)
  -colSums(z^2)/2-sum(log(diag(R))) - log(2*pi)*length(mu)/2
}






#Now some code to extract the joint density of data and coefs from a gam...

bSb <- function(b,beta=coef(b)) {
  ## evaluate penalty for fitted gam, possibly with new beta
  bSb <- k <-  0
  for (i in 1:length(b$smooth)) {
    m <- length(b$smooth[[i]]$S)
    if (m) {
      ii <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      for (j in 1:m) {
        k <- k + 1
        bSb <- bSb + b$sp[k]*(t(beta[ii])%*%b$smooth[[i]]$S[[j]]%*%beta[ii])
      }
    }  
  }
  bSb
} ## bSb

devg <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## evaluate the deviance of a fitted gam given possibly new coefs, beta
  ## for general families this is simply -2*log.lik
  if (inherits(b$family,"general.family")) {
    -2*b$family$ll(b$y,X,beta,b$prior.weights,
                   b$family,offset=b$offset)$l
  } else { ## exp or extended family
    sum(b$family$dev.resids(b$y,b$family$linkinv(X%*%beta+b$offset),
                            b$prior.weights))
  }
} ## devg

lpl <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## log joint density for beta, to within uninteresting constants
  -(devg(b,beta,X)+bSb(b,beta))/2
}


#Finally a simple sampler. It alternates proposals based directly on the fixed approximate posterior (but with heavier tails) with random walk proposals based on a shrunk version of the posterior covariance matrix. The fixed part promotes big jumps, the random walk part avoids getting stuck.

gam.mh <- function(b, ns=10000, burn=4000, 
                   t.df=4, rw.scale=.25, thin=1) {
  ## generate posterior samples for fitted gam using Metroplois Hastings sampler
  ## alternating fixed proposal and random walk proposal, both based on Gaussian
  ## approximation to posterior...
  beta <- coef(b);Vb <- vcov(b)
  X <- model.matrix(b); burn <- max(0,burn)
  bp <- rmvt(ns+burn,beta,Vb,df=t.df) ## beta proposals
  bs[1,] <- beta
  lfp <- dmvt(t(bp),beta,Vb,df=t.df) ## log proposal density
  
  R <- chol(Vb) 
  step <- rmvn(ns+burn,beta*0,Vb*rw.scale) ## random walk steps (mgcv::rmvn)
  
  u <- runif(ns+burn);us <- runif(ns+burn) ## for acceptance check
  bs <- bp;j <- 1;accept <- rw.accept <- 0
  
  lpl0 <- lpl(b,bs[1,],X)
  for (i in 2:(ns+burn)) { ## MH loop
    ## first a static proposal...
    lpl1 <- lpl(b,bs[i,],X)
    if (u[i] < exp(lfp[j]-lfp[i]+lpl1-lpl0)) {
      lpl0 <- lpl1;accept <- accept + 1
      j <- i ## row of bs containing last accepted beta
    } else bs[i,] <- bs[i-1,]
    ## now a random walk proposal...
    lpl1 <- lpl(b,bs[i,]+step[i,],X)
    if (us[i] < exp(lpl1-lpl0)) { ## accept random walk step
      lpl0 <- lpl1;j <- i
      bs[i,] <- bs[i,] + step[i,]
      rw.accept <- rw.accept+1 
      lfp[i] <- dmvt(bs[i,],beta,Vb,df=4,R=R) ## have to update static proposal density
    }
    if (i==burn) accept <- rw.accept <- 0
    if (i%%2000==0) cat(".")
  } ## MH loop
  if (burn>0) bs <- bs[-(1:burn),]
  if (thin>1) bs <- bs[seq(1,ns,by=thin),]
  list(bs=bs,rw.accept = rw.accept/ns,accept=accept/ns)
} ## gam.mh
