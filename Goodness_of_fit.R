source("./Tulap.R")
library(rmutil)
library(dgof)
library(ROptEstOld)

KS <- function(x,cdf,ep,Kuiper = FALSE){
  n = length(x)
  ks = ks.test(x,cdf,alternative = "two.sided")$statistic
  if(Kuiper == "Kuiper"){
    Dplus = ks.test(x,cdf,alternative="greater")$statistic
    Dminus = ks.test(x,cdf,alternative="less")$statistic
    ks = Dplus + Dminus
  }
  N = (1/n)*rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
  ks = ks+N
  return(ks)
}

Unknown_GOF<-function(x,ep,Kuiper = FALSE){
  #the minimum KS parameter estimates, estimate mean and sd of normal based on the sample x
  KSEstimator <- function(x){
    KSDis <- function(param){
      ks_min = ks.test(x,"pnorm",mean =param[1],sd = param[2],alternative="two.sided")$statistic
      if(Kuiper == "Kuiper"){
        Dplus = ks.test(x,"pnorm",mean =param[1],sd = param[2],alternative="greater")$statistic
        Dminus = ks.test(x,"pnorm",mean =param[1],sd = param[2],alternative="less")$statistic
        ks_min = Dplus+Dminus
      }
      return(ks_min)
    }
    res <- optim(c(0,1), fn = KSDis)$par
    return(list(mean = res[1], sd = res[2]))
  }
  ks_param = KSEstimator(x=x)
  fitted_cdf = function(t){
    return(pnorm(t,mean = ks_param$mean, sd = ks_param$sd))
  }
  #Then to evaluate the distance we can just use the above implementation of KS:
  ks = KS(x,fitted_cdf,ep = ep,Kuiper = Kuiper)
  return(ks)
}

Cramer <- function(x,cdf,ep) {
  U = sort(x)
  n = length(x)
  rank = seq_len(n)
  ###JA "U" should be replaced with cdf(U)
  omega2 = 1/(12 * n) + sum(((2*rank - 1)/(2*n)-cdf(U))^2)
  omega2 = sqrt(omega2/n)
  ### Add (1/n)*rlaplace(s=1/ep)
  omega2 = omega2 + (1/n)*rlaplace(s=1/ep)
  return(omega2)
}

cdf <- function(q){
  return (pnorm(q,m=0,s=1))
}

pdf<- function(x){
  return(dnorm(x,0,1))
}

# Wasserstein <- function(x,cdf,pdf,ep){
#   # note that x is not an input to Integrand, 
#   # but because it is defined inside of Wasserstein, it will know what x is
#   Integrand = function(t){
#     x_ecdf = ecdf(x) 
#     return(abs(cdf(t)-x_ecdf(t))*pdf(t))
#   }
#   n = length(x)
#   w = integrate(Integrand,lower = -Inf,upper = Inf,subdivisions=2000)
#   w = w$value + (1/n)*rlaplace(s=1/ep)
#   return(w)
# }

reference = function(n,ep,reps,type,Kuiper = FALSE){
  x = rnorm(n*reps,0,1)
  x_mat = matrix(x,nrow=reps,ncol=n)
  if(type == "cramer")
    return(apply(X=x_mat,MARGIN=1, FUN=Cramer,cdf = pnorm,ep=ep))
  if(type == "ks")
    return(apply(X=x_mat,MARGIN=1, FUN=KS,cdf = pnorm,ep=ep,Kuiper = Kuiper))
  # if(type == "Wasserstein")
  #   return(apply(X=x_mat,MARGIN=1, FUN=Wasserstein,cdf = pnorm, pdf = dnorm, ep=ep))
  if(type == "Unknown")
    return(apply(X=x_mat,MARGIN=1, FUN=Unknown_GOF, ep=ep,Kuiper = Kuiper))

}
###JA: Anderson-Darling test may not be possible, since sensitivity seems to be unbounded

ep= .1
reps = 10
al = 0.05
nVec = c(50,100,200,400,800,1600)
#nVec = c(25,50,100,200,400)
p_GOF_KS = p_Cramer =p_KS_V = p_WA = p_Unknown = p_Unknown_Kuiper = rep(0,reps)
power_GOF_KS = power_Cramer = power_ksv = power_wa = power_Unknown = power_Unknown_Kuiper= rep(0,length(nVec))
for (k in 1:length(nVec)){
  print(k/length(nVec))
  n = nVec[k]
  null_cramer = ecdf(reference(n=n,ep=ep,reps=1000,type = "cramer"))
  null_ks= ecdf(reference(n=n,ep=ep,reps=1000,type = "ks"))
  null_ksv= ecdf(reference(n=n,ep=ep,reps=1000,type = "ks",Kuiper = "Kuiper"))
  #null_wa= ecdf(reference(n=n,ep=ep,reps=10000,type = "Wasserstein"))
  null_Unknown= ecdf(reference(n=n,ep=ep,reps=1000,type = "Unknown"))
  null_Unknown_Kuiper= ecdf(reference(n=n,ep=ep,reps=1000,type = "Unknown",Kuiper = "Kuiper"))
  for(i in 1:reps){
    x = rnorm(n,mean = 1,sd = 1)
    #x = rcauchy(n,0,1)
    #x = rlaplace(n,0,1)
    #x = rlogis(n,0,1)
    cramer = Cramer(x,cdf,ep) 
    ks = KS(x,cdf,ep,Kuiper = FALSE)
    ksv = KS(x,cdf,ep,Kuiper = "Kuiper")
    #wa = Wasserstein(x,cdf,pdf,ep)
    unknown = Unknown_GOF(x,ep,Kuiper = FALSE)
    unknown_kuiper = Unknown_GOF(x,ep,Kuiper = "Kuiper")
    p_Cramer[i] = 1-null_cramer(cramer)
    p_GOF_KS[i] = 1-null_ks(ks)
    p_KS_V[i] = 1-null_ksv(ksv)
    #p_WA[i] = 1-null_wa(wa)
    p_Unknown[i] = 1-null_Unknown(unknown)
    p_Unknown_Kuiper[i] = 1-null_Unknown_Kuiper(unknown_kuiper)
  }
  power_GOF_KS[k]=mean(p_GOF_KS < al)
  power_Cramer[k]=mean(p_Cramer < al)
  power_ksv[k]=mean(p_KS_V < al)
  #power_wa[k]=mean(p_WA < al)
  power_Unknown[k]=mean(p_Unknown < al)
  power_Unknown_Kuiper[k]=mean(p_Unknown_Kuiper < al)
}

power_Cramer
power_GOF_KS
power_ksv
power_Unknown
power_Unknown_Kuiper



plot(log(nVec),power_Cramer,col="black",ylim=c(0,1),type="l",xaxt="n",ylab = "power",xlab="n",lty=1,lwd = 3)
lines(log(nVec),power_GOF_KS,col="red",lty=2,lwd = 3)
lines(log(nVec),power_ksv,col="green",lty=3,lwd = 3)
#lines(log(nVec),power_wa,col="blue",lty=6,lwd = 1.5)
lines(log(nVec),power_Unknown,col="orange",lty=4,lwd = 3)
lines(log(nVec),power_Unknown_Kuiper,col="blue",lty=5,lwd = 3)
axis(at=log(nVec),labels=nVec,side=1)
legend("bottomright",c("Cramer","KS","Kuiper","KS_Unknown","Kuiper_Unknown"),col=c("black","red","green","orange","blue"),lty=c(1,2,3,4,5),bty="o",lwd = 3)
legend("bottomright",c("Cramer","KS","Kuiper"),col=c("black","red","green"),lty=c(1,2,3),bty="o",lwd = 3)
