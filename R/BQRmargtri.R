
##############################################################
####     MCMC code to run the marginal model of Reich,   #####
####     Bondell, and Wang for clustered data.           #####
####     The main function is BayesQReg                  #####
####                                                     #####
####     An analysis of a sample data set is             #####
####     conducted at the bottom of this file.           #####
##############################################################:


BayesQReg<-function(y,X,subject,tau=0.5,
           runs=5000,burn=1000,
           aD=1,bD=1,M=25,prec.beta=0.01,mx.sig=NA){   

  #y:          vector of observations
  #X:          matrix of predictors (includes the intercept!)
  #subject:    vector of subject labels
  #tau:        quantile of interest
  #siglim:     sigma~Unif(0,sigmlim)
  #aD,bD:      V_k~Beta(1,D);D~gamma(aD,bD)
  #prec.beta   beta[j]~N(0,1/sqrt(prec.beta)
  #runs:       number of MCMC iterations
  #burn:       number of MCMC iterations to discard as burn-in
  #M:          Number of terms in the finite approx

  n<-length(y)
  p<-ncol(X)
  nsub<-max(subject)

  #Generate initial values
  D<-1
  mn.sig<-0.02
  if(is.na(mx.sig)){mx.sig<-2*sd(y)}
  
  beta<-rnorm(p)
  alpha<-rep(0,nsub)
  siga<-.01
  gamma<-rep(0,p);gamma[1]<-1; gamma[4] <- 1; gamma[7] <- 1;
  sigh<-as.vector(X%*%gamma)
  lambda<-10
  
  q<-sig1<-sig2<-mu1<-mu2<-rep(2,M)
  for(j in 1:M){while(q[j]<0 | q[j]>1){
     mu1[j]<-rASL(1,lambda,tau)
     mu2[j]<-rASL(1,lambda,tau)
     sig1[j]<-runif(1,mn.sig,mx.sig)
     sig2[j]<-runif(1,mn.sig,mx.sig)
     q[j]<-make.q(mu1[j],mu2[j],sig1[j],sig2[j],siga,tau)
  }}

  v<-rbeta(M,1,D);v[M]<-1
  probs<-makeprobs(v)
#  g<-rcat(n,probs)
  g<-c(rcat(n,probs))
  h<-rbinom(n,1,q[g])    #h=I(mu<0) P(h=1)=q

  #keep track of stuff
  keepbeta<-keepgamma<-matrix(0,runs,p)
  keepalpha<-matrix(0,runs,nsub)
  keeplambda<-keepD<-keepsiga<-keepICC<-rep(0,runs)
  TruncProb<-rep(0,n)
  x.grid<-seq(-10,5,0.01)
  sumdense<-sumdense2<-0*x.grid
  acc<-att<-rep(0,7);can<-rep(.25,7)
  count<-afterburn<-0


  #LET'S ROLL!
  for(i in 1:runs){

   #update beta:
    sd<-sigh*(h*sig1[g]+(1-h)*sig2[g])
    taue<-1/sd^2
    var<-solve(t(X)%*%diag(taue)%*%X+prec.beta*diag(p))
    mn<-var%*%t(X)%*%(taue*(y-sigh*(alpha[subject]+h*mu1[g]+(1-h)*mu2[g])))
    beta<-mn+t(chol(var))%*%rnorm(p)
    resids<-y-X%*%beta-sigh*(h*mu1[g]+(1-h)*mu2[g])

    #Update alpha
    taue<-1/(h*sig1[g]+(1-h)*sig2[g])^2
    for(j in 1:nsub){
       VAR<-1/(sum(taue[subject==j])+1/siga^2)
       MN<-sum(taue[subject==j]*resids[subject==j]/sigh[subject==j]) 
       alpha[j]<-rnorm(1,MN*VAR,sqrt(VAR))
    }

    cansiga<-rnorm(1,siga,.25*can[7])
    att[7]<-att[7]+1
    if(cansiga>0 & cansiga<mx.sig){
        canq<-q
        for(j in 1:M){canq[j]<-make.q(mu1[j],mu2[j],sig1[j],sig2[j],cansiga,tau)
        if(min(canq)>0 & max(canq)<1){
            MHrate<-sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                    sum(h*log(q[g])+(1-h)*log(1-q[g]))+
                    sum(dnorm(alpha,0,cansiga,log=T))-
                    sum(dnorm(alpha,0,siga,log=T))
            if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
                 q<-canq;siga<-cansiga;acc[7]<-acc[7]+1}}
          }
        }
    }
    taua<-1/siga^2
    resids<-y-X%*%beta

    #update gamma:
    mmm<-(1-h)*mu2[g]+h*mu1[g]+alpha[subject]
    sss<-h*sig1[g]+(1-h)*sig2[g]
    if(p>1){for(s in 2:p){
      if (s!=4 && s!=7) {
        att[6]<-att[6]+1
        cangamma<-gamma;cangamma[s]<-rnorm(1,gamma[s],2*can[6])
        cansigh<-as.vector(X%*%cangamma)
        if(min(cansigh)>0){
          MHrate<-sum(dnorm(resids,cansigh*mmm,cansigh*sss,log=T)-
                      dnorm(resids,sigh*mmm,sigh*sss,log=T))+
                        dnorm(cangamma[s],0,1,log=T)-dnorm(gamma[s],0,1,log=T)
          if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
            gamma<-cangamma;sigh<-cansigh;acc[6]<-acc[6]+1}} 
        }
      }
    }}
    gamma[1]<-1
    gamma[4] <- 1
    gamma[7] <- 1
   sigh<-as.vector(X%*%gamma)
   resids<-y-X%*%beta-sigh*alpha[subject]

    if (sum(is.nan(sigh)) !=0 ) print('nan, sigh')
    
    #update g,h
    for(s in 1:n){
       if(h[s]==0){lll<-log(1-q)+dnorm(resids[s],sigh[s]*mu2,sigh[s]*sig2,log=T)+log(probs)}
       if(h[s]==1){lll<-log(q)+dnorm(resids[s],sigh[s]*mu1,sigh[s]*sig1,log=T)+log(probs)}
       g[s]<-rcat(1,exp(lll-max(lll)))
    }
    p1<-q[g]*dnorm(resids,sigh*mu1[g],sigh*sig1[g])
    p0<-(1-q[g])*dnorm(resids,sigh*mu2[g],sigh*sig2[g])
    h<-rbinom(n,1,p1/(p0+p1))

    for(j in 1:M){
 
      #update mu1
       att[1]<-att[1]+(probs[j]>.2)
       canmu1<-mu1;canmu1[j]<-rnorm(1,mu1[j],can[1])
       canq<-q;canq[j]<-make.q(canmu1[j],mu2[j],sig1[j],sig2[j],siga,tau)
       if(canq[j]>0 & canq[j]<1){
          MHrate<-dASL(canmu1[j],0,lambda,tau)-dASL(mu1[j],0,lambda,tau)+
                  sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                  sum(h*log(q[g])+(1-h)*log(1-q[g]))
          if(sum(g==j & h==1)>0){
              MHrate<-MHrate+
                      sum(dnorm(resids[g==j & h==1],sigh[g==j & h==1]*canmu1[j],sigh[g==j & h==1]*sig1[j],log=T)-
                          dnorm(resids[g==j & h==1],sigh[g==j & h==1]*mu1[j],sigh[g==j & h==1]*sig1[j],log=T))
          }
          if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
            acc[1]<-acc[1]+(probs[j]>.2);q<-canq;mu1<-canmu1}}
        }

        #update mu2
        att[2]<-att[2]+(probs[j]>.2)
        canmu2<-mu2;canmu2[j]<-rnorm(1,mu2[j],can[2])
        canq<-q;canq[j]<-make.q(mu1[j],canmu2[j],sig1[j],sig2[j],siga,tau)
        if(canq[j]>0 & canq[j]<1){
          MHrate<-dASL(canmu2[j],0,lambda,tau)-dASL(mu2[j],0,lambda,tau)+
                  sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                  sum(h*log(q[g])+(1-h)*log(1-q[g]))
          if(sum(g==j & h==0)>0){
              MHrate<-MHrate+
                      sum(dnorm(resids[g==j & h==0],sigh[g==j & h==0]*canmu2[j],sigh[g==j & h==0]*sig2[j],log=T)-
                          dnorm(resids[g==j & h==0],sigh[g==j & h==0]*mu2[j],sigh[g==j & h==0]*sig2[j],log=T))
          }
          if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
             acc[2]<-acc[2]+(probs[j]>.2);q<-canq;mu2<-canmu2}}
        }

        r<-resids-sigh*(h*mu1[g]+(1-h)*mu2[g])
        #update sigma1
        att[3]<-att[3]+(probs[j]>.2)
        cansig1<-sig1;cansig1[j]<-rnorm(1,sig1[j],can[3])
        if(cansig1[j]>mn.sig & cansig1[j]<mx.sig){
          canq<-q;canq[j]<-make.q(mu1[j],mu2[j],cansig1[j],sig2[j],siga,tau)
          if(canq[j]>0 & canq[j]<1){
            MHrate<-sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                    sum(h*log(q[g])+(1-h)*log(1-q[g]))
            if(sum(g==j & h==1)>0){
                MHrate<-MHrate+
                        sum(dnorm(r[g==j & h==1],0,sigh[g==j & h==1]*cansig1[j],log=T)-
                            dnorm(r[g==j & h==1],0,sigh[g==j & h==1]*sig1[j],log=T))
             }
            if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
              acc[3]<-acc[3]+(probs[j]>.2);q<-canq;sig1<-cansig1}}
           }
        }

        #update sigma2
        att[3]<-att[3]+(probs[j]>.2)
        cansig2<-sig2;cansig2[j]<-rnorm(1,sig2[j],can[3])
        if(cansig2[j]>mn.sig & cansig2[j]<mx.sig){
          canq<-q;canq[j]<-make.q(mu1[j],mu2[j],sig1[j],cansig2[j],siga,tau)
          if(canq[j]>0 & canq[j]<1){
            MHrate<-sum(h*log(canq[g])+(1-h)*log(1-canq[g]))-
                    sum(h*log(q[g])+(1-h)*log(1-q[g]))
            if(sum(g==j  & h==0)>0){
                MHrate<-MHrate+
                        sum(dnorm(r[g==j & h==0],0,sigh[g==j & h==0]*cansig2[j],log=T)-
                            dnorm(r[g==j & h==0],0,sigh[g==j & h==0]*sig2[j],log=T))
            }
           if(!is.na(MHrate)){if(runif(1)<exp(MHrate)){
             acc[3]<-acc[3]+(probs[j]>.2);q<-canq;sig2<-cansig2}}
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
    keeplambda[i]<-lambda
    keepalpha[i,]<-alpha
    keepsiga[i]<-siga
    keepD[i]<-D
    TruncProb[i]<-probs[M]
    E1<-sum(probs*q*mu1+probs*(1-q)*mu2)
    E2<-sum(probs*q*(mu1^2+sig1^2)+probs*(1-q)*(mu2^2+sig2^2))
    sig2eps<-E2-E1^2
    keepICC[i]<-siga*siga/(siga*siga+sig2eps)

    for(j in 1:length(can)){if(att[j]>50 & 2*i<burn){
      can[j]<-can[j]*ifelse(acc[j]/att[j]<0.3,0.5,1)*
                     ifelse(acc[j]/att[j]>0.5,1.5,1)
      acc[j]<-att[j]<-1
    }}

    if(i>=burn){
      ddd<-0*x.grid
      for(j in 1:M){
          ddd<-ddd+q[j]*probs[j]*dnorm(x.grid,mu1[j],sig1[j])+
                   (1-q[j])*probs[j]*dnorm(x.grid,mu2[j],sig2[j])
          sumdense<-sumdense+ddd/(runs-burn)
          sumdense2<-sumdense2+ddd*ddd/(runs-burn)
      }
    }

    if(i%%100==0){
#      print(paste("Done with",i,"iterations"))
    }
  }

list(beta=keepbeta,gamma=keepgamma,alpha=keepalpha,D=keepD,lambda=keeplambda,
    MHrate=acc/att,TruncProb=TruncProb,ICC=keepICC,siga=keepsiga,
    x.grid=x.grid,dense.mn=sumdense,dense.sd=sqrt(abs(sumdense2-sumdense^2)), runs=runs, burn=burn)}

#define other functions called by BayesQReg:

#compute the stick-breaking weights
makeprobs<-function(v){ 
   #compute the stick-breaking weights
   N<-length(v)
   probs<-v
   probs[2:N]<-probs[2:N]*cumprod(1-v[2:N-1])
   probs}

make.q<-function(mu1,mu2,sig1,sig2,siga,tau){
  (tau-pnorm(0,mu2,sqrt(sig2^2+siga^2)))/
  (pnorm(0,mu1,sqrt(sig1^2+siga^2))-pnorm(0,mu2,sqrt(sig2^2+siga^2)))}

rcat<-function(ndraws,prob){
  #generate categorical variables:
  prob[is.na(prob)]<-0
  if(sum(prob)==0){prob[1]<-1}
  (1:length(prob))%*%rmultinom(ndraws,1,prob)
}

#Asymmetric Laplace density:
dASL<-function(y,mu,lambda,tau){
   length(y)*log(lambda) - sum(ifelse(y>mu,tau,1-tau)*abs(y-mu)*lambda)}

#generate and asymmetric rv:
rASL<-function(n,lambda,tau){
   posneg<-rbinom(n,1,tau)
   (1-posneg)*rexp(n,tau*lambda)-posneg*rexp(n,(1-tau)*lambda)}


BQR.margtri.Summary <- function(foo, truebetatau, comp){

  runs <- dim(foo$beta)[1]
  beta.coef <- apply(foo$beta[foo$burn:foo$runs,],2,mean)
#  lbd <- apply(foo$beta[1000:5000,], 2, function(x) quantile(x, 0.025))
#  ubd <- apply(foo$beta[1000:5000,], 2, function(x) quantile(x, 0.975))
#  len <- ubd-lbd
  if (comp==1) beta <- beta.coef[1:3]
  if (comp==2) beta <- beta.coef[4:6]
  if (comp==3) beta <- beta.coef[7:9]
#  mse <- mean((beta[-1]-truebetatau[-1])^2)
  mse <- (beta-truebetatau)^2
#  p <- length(truebetatau)
#  for (j in 2:p)   cover <- prod((truebetatau[j]>lbd[j] && truebetatau[j]<ubd[j]))

#  return(list(beta=beta.coef, lbd=lbd, ubd=ubd, length=len, mse=mse, cover=cover))
  return(list(mse=mse))
}


##############################################################:
##############       NUMERICAL EXAMPLE             ###########:
##############################################################:

## x<-matrix(rnorm(500),100,5)
## subject<-sample(1:20,100,replace=T)
## alpha<-rnorm(20)
## x[,1]<-1
## y<-alpha[subject]+x[,1]+2*x[,2]+rgamma(100,2,2)

## fit<-BayesQReg(y,x,subject,tau=0.9)

## print("Regression coefficients (beta)")
## print(round(apply(fit$beta[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))

## print("Coefficients in the standard deviation (gamma)")
## print(round(apply(fit$gamma[1000:5000,],2,quantile,c(0.025,0.5,0.975)),3))

## print("Stick-breaking truncation probability")
## print(round(quantile(fit$TruncProb[1000:5000],c(0.025,0.5,0.975)),3))


## plot(fit$x.grid,fit$dense.mn,type="l",
##      xlab="x",ylab="residual density",xlim=c(-5,5))


## ## my simulation and comparison
## subject <- sample(1:50,100,replace=T)

## begin <- proc.time()
## fit <- reichBayesQReg(y,x,subject, tau= 0.5)
## proc.time()-begin

