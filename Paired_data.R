source("./Tulap.R")
library(dplyr)


sign_test_DP = function(x,y,ep){
  n = length(x)
  S = sum(y>x)
  if(ep<Inf){
  # generate random sample from Tulap
  N = rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
  T = S+N
  values = seq(0,n)
  p_val = ptulap(values-T,median=0,lambda=exp(-ep),cut=0)%*%dbinom(values,size=n,prob=1/2)
  }
  # else if(ep==Inf){
  # N = runif(n=1,-1/2,1/2)
  # #rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
  # T = S+N
  # values = seq(0,n)
  # p_val = punif(values-T,min=-1/2,max=1/2)%*%dbinom(values,size=n,prob=1/2)
  # }
  return(p_val)
}


median_test_DP = function(x,y,ep){
  n=length(x)
  z = c(x,y)
  T = sum(rank(z)[1:n]>n)
  ###  test
  if(ep<Inf){
    N = rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
    T = S+N
    values = seq(0,n)
    p_val = ptulap(values-T,median=0,lambda=exp(-ep),cut=0)%*%dhyper(m=n,k=n)
  }
  else if(ep==Inf){
    N = runif(n=1,-1/2,1/2)
    T = S+N
    values = seq(0,n)
    p_val = punif(values-T,min=-1/2,max=1/2)%*%dhyper(m=n,k=n)
  }
  return(p_val)
  ###   more to do
}

# `wc_new` is the Pratt variant of the test, which our test uses
# this variant essentially keeps all of the zeroes in the ranking step, and 
# then replaces their rankings with zero after the ranking takes place.
# the input, `x`, is a vector of d_i's

wc_new <- function(x,y,ep) {
  z = y-x
  n=length(z)
  ranks <- z %>% abs() %>% rank(ties.method = "average")
  ranks[z == 0] <- 0
  W <- sum(ranks * sign(z))
  N = rlaplace(n=1,m=0,s=2*n/ep)
  return(W+N)
}
# `gen_null` generates a null distribution of test statistics with a given 
# sample size `n` and privacy parameter `epsilon`. reps is the
# number of simulations to carry out.
gen_null <- function(n, epsilon, reps) {
  s <- sqrt(n*(n+1)*(2*n+1)/6)
  Z <- rnorm(n = reps,
             mean = 0,
             sd = s)
  noise <- rlaplace(n = reps,m = 0,s = 2*n/(epsilon))  
  return(noise + Z)
}

#modified kolmogorov-smirnov test
symmetrized_KS = function(x,y,ep,Kuiper = FALSE){
  z = y-x
  z_neg = -z
  ks = ks.test(z,z_neg,alternative = "two.sided")$statistic
  if(Kuiper == "Kuiper"){
    Dplus = ks.test(z,z_neg,alternative="greater")$statistic
    Dminus = ks.test(z,z_neg,alternative="less")$statistic
    ks = Dplus + Dminus
  }
  N = (2/n)*rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
  ks = ks+N
  return(ks)
}

reference_KS = function(n,ep,reps,Kuiper = FALSE){
  z = rnorm(n*reps,0,1)
  z_mat = matrix(z,nrow=reps,ncol=n)
  # margin 1对应每一列，2对应每一行
  return(apply(X=z_mat,MARGIN=1, FUN=symmetrized_KS,x=0,ep=ep,Kuiper = Kuiper))
}

# power calculation

ep= 10
reps = 10000
nVec = c(50,100,200,400,800,1600)
al=.05
power_signDP = power_wcDP  = power_ks = power_kuiper  = rep(0,length(nVec))
p_signDP = p_wcDP  = p_ks =  p_kuiper =  rep(0,reps)
# p_sign = p_wc = p_ks = rep(0,reps)

for (k in 1:length(nVec)){
  print(k/length(nVec))
  n = nVec[k]
  null_wcDP = ecdf(gen_null(n=n,epsilon=ep,reps=10000))
  #null_wc = function(q){return(pnorm(q,m=0,s=sqrt(n*(n+1)*(2*n+1)/6)))}
  null_ksDP = ecdf(reference_KS(n=n,ep=ep,reps=10000))
  null_kuiper = ecdf(reference_KS(n=n,ep=ep,reps=10000,Kuiper = "Kuiper"))

  # calculate the p-value
  for(i in 1:reps){
    #print(i/reps)
    x = rep(0,n)# rnorm(n,0,5)###   this doesnt matter
    y = rexp(n,1)-log(2)/1
    #y = x+rcauchy(n,.2,1)
    #y=x+rnorm(n,.2,1)
    #y = x+rlaplace(n,.2,1)
    #y = x+rlogis(n,.2,1)
    sign = sign_test_DP(x,y,ep)
    p_signDP[i] = 2*min(sign,1-sign)
    wc = wc_new(x,y,ep)
    p_wcDP[i] = 2*(1-null_wcDP(abs(wc)))
    ks = symmetrized_KS(x,y,ep)
    p_ks[i] = 1-null_ksDP(ks)
    kuiper = symmetrized_KS(x,y,ep,Kuiper = "Kuiper")
    p_kuiper[i] = 1-null_kuiper(kuiper)
     #sign = sign_test_DP(x,y,ep=Inf)
    #p_sign[i] = 2*min(sign,1-sign)
    #wc = wc_new(x,y,ep=Inf)
    #p_wc[i] = 2*(1-null_wc(abs(wc)))
    #ks = symmetrized_KS(x,y,ep=Inf)
    #p_ks[i] = 1-null_ks(ks)
  }
  power_signDP[k] = mean(p_signDP < 0.05)
  power_wcDP[k] = mean(p_wcDP < 0.05)
  power_ks[k] = mean(p_ks < 0.05)
  power_kuiper[k] = mean(p_kuiper < 0.05)

}
plot(log(nVec),power_signDP,type="l",col="black",ylim=c(0,1),lty=1,xaxt="n",ylab="Power",xlab="n",lwd = 1.7)
lines(log(nVec),power_wcDP,col="red",lty=3,lwd = 1.7)
lines(log(nVec),power_ks,col="green",lty=5,lwd = 1.7)
lines(log(nVec),power_kuiper,col="blue",lty=6,lwd = 1.7)
axis(at=log(nVec),labels=nVec,side=1)
legend("topleft",c("Sign","Wilcoxon","Ks","Kuiper"),col=c("black","red","green","blue"),lty=c(1,3,5,6),bty="o",lwd = 1.7)

