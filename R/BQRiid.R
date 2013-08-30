
#######################################################
####     MCMC code to run the model of Reich,     #####
####     Bondell, and Wang assuming independent   #####
####     errors.  The main function is BayesQReg  #####
####                                              #####
####     An analysis of a sample data set is      #####
####     conducted at the bottom of this file.    #####
#######################################################:



BayesQReg<-function(y,X,tau=0.5,M=10,
            runs=5000,burn=1000,
            aD=100,bD=100,
            mx.sig=NA,prec.beta=0.01){   
  #y:        vector of observations
  #X:        matrix of predictors (first column must be all ones!)
  #tau:      quantile of interest
  #M:        Number of terms in the finite approx
  #aD,bD:    V_k~Beta(1,D);D~dgamma(aD,bD)
  #siglim:   sigma~Unif(0,sigmlim)
  #prec.beta beta[j]~norm(0,1/sqrt(prec.beta))
  #runs:     number of MCMC iterations
  #burn:     number of MCMC iterations to discard as burn-in

  n<-length(y)
  p<-ncol(X)

  #Generate initial values
  D<-1
  mn.sig<-0.2
  if(is.na(mx.sig)){mx.sig<-2*sd(y)}
  beta<-rnorm(p)
  gamma<-rep(0,p);gamma[1]<-1
  sigh<-as.vector(X%*%gamma)
  lambda<-10

  q<-sig1<-sig2<-mu1<-mu2<-rep(2,M)
  for(j in 1:M){while(q[j]<0 | q[j]>1){
     mu1[j]<-rASL(1,lambda,tau)
     mu2[j]<-rASL(1,lambda,tau)
     sig1[j]<-runif(1,mn.sig,mx.sig)
     sig2[j]<-runif(1,mn.sig,mx.sig)
     q[j]<-make.q(mu1[j],mu2[j],sig1[j],sig2[j],tau)
  }}

  v<-rbeta(M,1,D);v[M]<-1
  probs<-makeprobs(v)
  g<-sample(1:M,n,replace=T,prob=probs)
  h<-rbinom(n,1,q[g])    #h=I(mu<0) P(h=1)=q

  #keep track of stuff
  keepbeta<-keepgamma<-matrix(0,runs,p)
  keepheter <- matrix(0, runs,p)
  TruncProb<-rep(0,n)
  x.grid<-seq(-20,20,0.01)
  sumdense<-sumdense2<-0*x.grid
  acc<-att<-rep(0,6);can<-rep(.25,6)
  count<-afterburn<-0

  #START THE MCMC SAMPLING:
  for(i in 1:runs){

   #update beta:
    sd<-sigh*(h*sig1[g]+(1-h)*sig2[g])
    var<-solve(t(X)%*%diag(1/sd^2)%*%X+prec.beta*diag(p))
    mn<-var%*%t(X)%*%((y-sigh*(h*mu1[g]+(1-h)*mu2[g]))/sd^2)
    beta<-mn+t(chol(var))%*%rnorm(p)
    resids<-y-X%*%beta

    #update gamma:
    mmm<-(1-h)*mu2[g]+h*mu1[g]
    sss<-h*sig1[g]+(1-h)*sig2[g]
    if(p>1){for(s in 2:p){
      att[6]<-att[6]+1
      cangamma<-gamma;cangamma[s]<-rnorm(1,gamma[s],2*can[6])
      cansigh<-as.vector(X%*%cangamma)
      keepheter[i, s] <- min(cansigh)
      if(min(cansigh)>0){
        MHrate<-sum(dnorm(resids,cansigh*mmm,cansigh*sss,log=T)-
                    dnorm(resids,sigh*mmm,sigh*sss,log=T))+
                dnorm(cangamma[s],0,1/sqrt(prec.beta),log=T)-dnorm(gamma[s],0,1/sqrt(prec.beta),log=T)
        if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
          gamma<-cangamma;sigh<-cansigh;acc[6]<-acc[6]+1}} 
      }
    }}
   gamma[1]<-1
   sigh<-as.vector(X%*%gamma)

    #update g,h
    for(s in 1:n){
       if(h[s]==0){lll<-log(1-q)+dnorm(resids[s],sigh[s]*mu2,sigh[s]*sig2,log=T)+log(probs)}
       if(h[s]==1){lll<-log(q)+dnorm(resids[s],sigh[s]*mu1,sigh[s]*sig1,log=T)+log(probs)}
       g[s]<-sample(1:M,1,prob=exp(lll-max(lll)))
    }
    p1<-q[g]*dnorm(resids,sigh*mu1[g],sigh*sig1[g])
    p0<-(1-q[g])*dnorm(resids,sigh*mu2[g],sigh*sig2[g])
    h<-rbinom(n,1,p1/(p0+p1))

    for(j in 1:M){
 
      #update mu1
       att[1]<-att[1]+(probs[j]>.2)
       canmu1<-mu1;canmu1[j]<-rnorm(1,mu1[j],can[1])
       canq<-q;canq[j]<-make.q(canmu1[j],mu2[j],sig1[j],sig2[j],tau)
       if(canq[j]>0 & canq[j]<1){
          MHrate<-dASL(canmu1[j],0,lambda,tau)-dASL(mu1[j],0,lambda,tau)+
                  sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                  sum(h*log(q[g])+(1-h)*log(1-q[g]))
          if(sum(g==j & h==1)>0){
              MHrate<-MHrate+
                      sum(dnorm(resids[g==j & h==1],sigh[g==j & h==1]*canmu1[j],sigh[g==j & h==1]*sig1[j],log=T)-
                          dnorm(resids[g==j & h==1],sigh[g==j & h==1]*mu1[j],sigh[g==j & h==1]*sig1[j],log=T))
          }
          if(runif(1)<exp(MHrate)){acc[1]<-acc[1]+(probs[j]>.2);q<-canq;mu1<-canmu1}
        }

        #update mu2
        att[2]<-att[2]+(probs[j]>.2)
        canmu2<-mu2;canmu2[j]<-rnorm(1,mu2[j],can[2])
        canq<-q;canq[j]<-make.q(mu1[j],canmu2[j],sig1[j],sig2[j],tau)
        if(canq[j]>0 & canq[j]<1){
          MHrate<-dASL(canmu2[j],0,lambda,tau)-dASL(mu2[j],0,lambda,tau)+
                  sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                  sum(h*log(q[g])+(1-h)*log(1-q[g]))
          if(sum(g==j & h==0)>0){
              MHrate<-MHrate+
                      sum(dnorm(resids[g==j & h==0],sigh[g==j & h==0]*canmu2[j],sigh[g==j & h==0]*sig2[j],log=T)-
                          dnorm(resids[g==j & h==0],sigh[g==j & h==0]*mu2[j],sigh[g==j & h==0]*sig2[j],log=T))
          }
          if(runif(1)<exp(MHrate)){acc[2]<-acc[2]+(probs[j]>.2);q<-canq;mu2<-canmu2}
        }

        r<-resids-sigh*(h*mu1[g]+(1-h)*mu2[g])
        #update sigma1
        att[3]<-att[3]+(probs[j]>.2)
        cansig1<-sig1;cansig1[j]<-rnorm(1,sig1[j],can[3])
        if(cansig1[j]>mn.sig & cansig1[j]<mx.sig){
          canq<-q;canq[j]<-make.q(mu1[j],mu2[j],cansig1[j],sig2[j],tau)
          if(canq[j]>0 & canq[j]<1){
            MHrate<-sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                    sum(h*log(q[g])+(1-h)*log(1-q[g]))
            if(sum(g==j & h==1)>0){
                MHrate<-MHrate+
                        sum(dnorm(r[g==j & h==1],0,sigh[g==j & h==1]*cansig1[j],log=T)-
                            dnorm(r[g==j & h==1],0,sigh[g==j & h==1]*sig1[j],log=T))
             }
            if(runif(1)<exp(MHrate)){acc[3]<-acc[3]+(probs[j]>.2);q<-canq;sig1<-cansig1}
           }
        }

        #update sigma2
        att[3]<-att[3]+(probs[j]>.2)
        cansig2<-sig2;cansig2[j]<-rnorm(1,sig2[j],can[3])
        if(cansig2[j]>mn.sig & cansig2[j]<mx.sig){
          canq<-q;canq[j]<-make.q(mu1[j],mu2[j],sig1[j],cansig2[j],tau)
          if(canq[j]>0 & canq[j]<1){
            MHrate<-sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                    sum(h*log(q[g])+(1-h)*log(1-q[g]))
            if(sum(g==j  & h==0)>0){
                MHrate<-MHrate+
                        sum(dnorm(r[g==j & h==0],0,sigh[g==j & h==0]*cansig2[j],log=T)-
                            dnorm(r[g==j & h==0],0,sigh[g==j & h==0]*sig2[j],log=T))
            }
           if(runif(1)<exp(MHrate)){acc[3]<-acc[3]+(probs[j]>.2);q<-canq;sig2<-cansig2}
          }
        }
    }
  
    #update lambda
    SSS<-sum(abs(mu1)*ifelse(mu1<0,1-tau,tau))+
         sum(abs(mu2)*ifelse(mu2<0,1-tau,tau))
    lambda<-rgamma(1,2*M+0.01,rate=SSS+0.01)

    #new v
    for(j in 1:(M-1)){v[j]<-rbeta(1,sum(g==j)+1,sum(g>j)+D)}
    v[M]<-1
    probs<-makeprobs(v)

    #new D:
    D<-rgamma(1,M-1+aD,bD-sum(log(1-v[-M])))
    if(is.na(D)){D<-1}
    if(D<0.001){D<-0.001}

    #keep track of stuff
    keepbeta[i,]<-beta
    keepgamma[i,]<-gamma
    TruncProb[i]<-probs[M]

    for(j in 1:6){if(att[j]>50 & 2*i<burn){
      can[j]<-can[j]*ifelse(acc[j]/att[j]<0.3,0.5,1)*
                     ifelse(acc[j]/att[j]>0.5,1.5,1)
      acc[j]<-att[j]<-1
    }}


    if(i>burn){
      ddd<-0*x.grid
      for(j in 1:M){
          ddd<-ddd+q[j]*probs[j]*dnorm(x.grid,mu1[j],sig1[j])+
                   (1-q[j])*probs[j]*dnorm(x.grid,mu2[j],sig2[j])
          sumdense<-sumdense+ddd/(runs-burn)
          sumdense2<-sumdense2+ddd*ddd/(runs-burn)
      }
    }

#    if(i%%1000==0){print(i)}
  }

list(beta=keepbeta,gamma=keepgamma,
    TruncProb=TruncProb,
    x.grid=x.grid,dense.mn=sumdense,dense.var=sumdense2-sumdense^2,
     heter=keepheter, runs=runs, burn=burn)
}


#define other functions called by BayesQReg:

#compute the stick-breaking weights
makeprobs<-function(v){ 
   N<-length(v)
   probs<-v
   probs[2:N]<-probs[2:N]*cumprod(1-v[2:N-1])
probs}

make.q<-function(mu1,mu2,sig1,sig2,tau){
  (tau-pnorm(0,mu2,sig2))/(pnorm(0,mu1,sig1)-pnorm(0,mu2,sig2))}


#Asymmetric Laplace density:
dASL<-function(y,mu,lambda,tau){
   length(y)*log(lambda) - sum(ifelse(y>mu,tau,1-tau)*abs(y-mu)*lambda)}

#generate and asymmetric rv:
rASL<-function(n,lambda,tau){
   posneg<-rbinom(n,1,tau)
(1-posneg)*rexp(n,tau*lambda)-posneg*rexp(n,(1-tau)*lambda)}





##############################################################:
##############       NUMERICAL EXAMPLE             ###########:
##############################################################:

## x<-matrix(rnorm(500),100,5)
## x[,1]<-1
## y<-x[,1]+2*x[,2]+rgamma(100,2,2)

## fit<-BayesQReg(y,x,tau=0.9)

## print("Regression coefficients (beta)")
## print(round(apply(fit$beta[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))

## print("Coefficients in the standard deviation (gamma)")
## print(round(apply(fit$gamma[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))

## print("Stick-breaking truncation probability")
## print(round(quantile(fit$TruncProb[1000:5000],c(0.025,0.5,0.975)),3))


## plot(fit$x.grid,fit$dense.mn,type="l",
##      xlab="x",ylab="residual density",xlim=c(-5,5))


############  my simulation 01/12 #############

## x<-matrix(rnorm(500),100,5)
## x[,1]<-1
## y<-x[,1]+2*x[,2]+rnorm(100, 100, 0.1)

## begin <- proc.time()
## fit<-BayesQReg(y,x,tau=0.5)
## proc.time()-begin

## print("Regression coefficients (beta)")
## print(round(apply(fit$beta[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))

## y = b0 + b1 *xi + ei , ei ~ N(0,1)
## tao = 0.5
## n = 1000 take long time
## true: 1 + 2*x2 + e
## > print("Regression coefficients (beta)")
## [1] "Regression coefficients (beta)"
## > print(round(apply(fit$beta[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))
##        [,1]  [,2]   [,3]   [,4]   [,5]
## 2.5%  1.026 1.908 -0.035 -0.038 -0.027
## 50%   1.090 1.968  0.028  0.024  0.036
## 97.5% 1.155 2.030  0.090  0.089  0.101

# tao = 0.25
# n = 100
# true: 1+ 2*x2 -0.67  -> 0.33+2*x2
## > print(round(apply(fit$beta[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))
##         [,1]  [,2]   [,3]   [,4]   [,5]
## 2.5%  -0.009 1.656 -0.188 -0.476 -0.114
## 50%    0.292 1.890  0.008 -0.135  0.119
## 97.5%  0.526 2.117  0.216  0.153  0.372

## tao = 0.75
##        [,1]  [,2]   [,3]   [,4]   [,5]
## 2.5%  1.479 1.615 -0.129 -0.420 -0.123
## 50%   1.728 1.828  0.071 -0.188  0.150
## 97.5% 2.014 2.049  0.276  0.082  0.402


## y<-x[,1]+2*x[,2]+rnorm(100, 100, 0.1)
## print(round(apply(fit$beta[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))
##          [,1]  [,2]   [,3]   [,4]   [,5]
## 2.5%  100.955 1.951 -0.035 -0.033 -0.019
## 50%   101.000 1.983  0.003  0.001  0.016
## 97.5% 101.042 2.018  0.040  0.032  0.045
## > 
##  proc.time()-begin
##    user  system elapsed 
## 251.664   0.392 253.093 

## begin <- proc.time()
## for (i in 1:1000000)  1
## proc.time()-begin
