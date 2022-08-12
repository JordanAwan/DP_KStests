source("./Tulap.R")
library(rmutil)
library(dgof)

Kruskal_Wallis <- function(x = NULL,y = NULL,z = NULL,ep) {
  # old kw
  # n1 = length(x)
  # n2 = length(y)
  # n = n1+n2
  # z = c(x,y)
  # term_1 = 12/(n*(n + 1))
  # l = n1+1
  # term_2 = (sum(rank(z)[1:n1])/n1)^2 + sum(rank(z)[l:length(z)]/n2)^2
  # term_3 = 3*(n + 1)
  # H  = term_1 * term_2 - term_3
  # H = H +87*rlaplace(s=1/ep)
  # couch version kw
  if (!is.null(z)) {
    x = z[1:(length(z)/2)]
    y = z[(length(z)/2+1):(length(z))]
  }
  n1 = length(x)
  n2 = length(y)
  z = c(x,y)
  l = n1+1
  n = n1+n2
  if ((n %% 2) == 0) {
    # even numbers
    H = 4*(n - 1) / (n^2) * (n1*abs((sum(rank(z)[1:n1])/n1)-(n+1)/2) + n2* abs(sum(rank(z)[l:length(z)]/n2)-(n+1)/2))
  }
  else {
    H = 4 / (n + 1) * (n1*abs((sum(rank(z)[1:n1])/n1)-(n+1)/2) + n2* abs(sum(rank(z)[l:length(z)]/n2)-(n+1)/2))
  }
  H = H +8*rlaplace(s=1/ep)
  return(H)
}

Cramer <- function(x,cdf,ep) {
  U = sort(x)
  n = length(x)
  rank = seq_len(n)
  ###JA "U" should be replaced with cdf(U)
  omega2 = 1/(12 * n) + sum(((2*rank - 1)/(2*n)-cdf(U))^2)
  ### Add (1/n)*rlaplace(s=1/ep)
  omega2 = omega2 + (1/n)*rlaplace(s=1/ep)
  return(omega2)
}

Mann_Whitney <- function(x = NULL,y=NULL,z = NULL, ep){
  if (!is.null(z)) {
    x = z[1:(length(z)/2)]
    y = z[(length(z)/2+1):(length(z))]
  }
  n1 = length(x)
  n2 = length(y)
  z = c(x,y)
  l = n1+1
  U1 = sum(rank(z)[1:n1],na.rm = TRUE)
  U2 = sum(rank(z)[l:length(z)],na.rm = TRUE)
  # U1 = n1*n2-sum(rank(z)[1:n1],na.rm = TRUE) + (n1*(n1 + 1)/2)
  # U2 = n1*n2-sum(rank(z)[l:length(z)],na.rm = TRUE) + (n2*(n2 + 1)/2)
  U = min(U1,U2)
  # U1 = (n1*(n1 + 1)/2) - sum(rank(z)[1:n1],na.rm = TRUE)
  # U2 = (n2*(n2 + 1)/2) - sum(rank(z)[l:length(z)],na.rm = TRUE)
  # U = min(U1, U2)
  #U = abs(U1)
  # m = min(n1,n2)
  # m_t = m+rlaplace(n=n,s=1/ep_m)
  # c = -ln(2*delta)/ep_m
  # m_star = max(ceiling(m_t-c),0)
  u_t = U + (n1)*rlaplace(n=1,m=0,s= 1/ep)### check this
  return(u_t)
}

median_test_DP <- function(x = NULL,y = NULL,z = NULL,ep=ep){
  if (!is.null(z)) {
    x = z[1:(length(z)/2)]
    y = z[(length(z)/2+1):(length(z))]
  }
  n=length(x)
  z = c(x,y)
  S = sum(rank(z)[1:n]>n)
  ###  test
  if(ep<Inf){
    N = rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
    T = S+N
    return(abs(T-n/2))
    # p_val = ptulap(T,median=0,lambda=exp(-ep),cut=0)*dhyper(x=n,m=n,n=n,k=n)
    # return(p_val)
  }
  # else if(ep==Inf){
  #   N = runif(n=1,-1/2,1/2)
  #   T = S+N
  #   values = seq(0,n)
  #   p_val = punif(values-T,min=-1/2,max=1/2)*dhyper(x=n,m=n,n=n,k=n)
  # }
  # return(p_val)
}

# Two_Sample_KS <- function(x = NULL,y = NULL,z = NULL,ep,Kuiper = FALSE){
#   # if z is not null
#   if (!is.null(z)) {
#     x = z[1:(length(z)/2)]
#     y = z[(length(z)/2+1):(length(z))]
#   }
#   x_ecdf = ecdf(x)
#   y_ecdf = ecdf(y)
#   z2 = c(pmin(x,0),pmin(y,0))
#   Dplus=max(c(x_ecdf(z2)-y_ecdf(z2),0))
#   Dminus=max(c(y_ecdf(z2)-x_ecdf(z2),0))
#   ### the old KS test statistic (Bn)
#   KS=max(Dplus,Dminus)
#   ### the new one (Vn)
#   if(Kuiper == "Vn")
#     KS = Dplus + Dminus
#   if(Kuiper == "Un")
#     KS=min(Dplus,Dminus)
#   if(ep == Inf)
#     return(KS+runif(n=1,-1/(2*n),1/(2*n)))
#   N = (1/n)*rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
#   return(KS+N)
# }

KS <- function(x = NULL,y = NULL,z = NULL,ep,Kuiper = FALSE){
  if (!is.null(z)) {
    x = z[1:(length(z)/2)]
    y = z[(length(z)/2+1):(length(z))]
  }
  KS = ks.test(x,y,alternative = "two.sided")$statistic
  Dplus = ks.test(x,y,alternative="greater")$statistic
  Dminus = ks.test(x,y,alternative="less")$statistic
  if(Kuiper == "Kuiper")
    KS = Dplus + Dminus
  N = (1/n)*rtulap(n=1, median = 0, lambda = exp(-ep), cut=0)
  return(KS+N)
  
}
gen_null <- function(n, epsilon, reps,laplace=TRUE) {
  s <- sqrt(n*(n+1)*(2*n+1)/6)
  Z <- rnorm(n = reps,
             mean = 0,
             sd = s)
  if (is.na(epsilon)) {
    # no epsilon return True
    return(Z)
  } else {
    if(laplace==TRUE)
      noise <- rlaplace(n = reps,
                        m = 0,
                        s = 2*n/(epsilon))  
    else if(laplace==FALSE)
      noise= 2*n*rtulap(n=reps,median=0,lambda=exp(-epsilon),cut=0)
    return(noise + Z)
  }
}

reference = function(n,ep,reps,type,Kuiper = FALSE){
  x = rnorm(n*reps,0,1)
  x_mat = matrix(x,nrow=reps,ncol=n)
  y = rnorm(n*reps,0,1)
  y_mat = matrix(y,nrow=reps,ncol=n)
  z_mat = cbind(x_mat,y_mat)
  # margin 1 for each column
  if(type == "ks")
    return(apply(X=z_mat,MARGIN=1, FUN=KS,x = NULL,y = NULL,ep=ep,Kuiper = Kuiper))
  if(type == "mn")
    return(apply(X=z_mat,MARGIN=1, FUN=Mann_Whitney,x = NULL,y = NULL,ep=ep))
  if(type == "kw")
    return(apply(X=z_mat,MARGIN=1, FUN=Kruskal_Wallis,x = NULL,y = NULL,ep=ep))
  if(type == "median")
    return(apply(X=z_mat,MARGIN=1, FUN=median_test_DP,x = NULL,y = NULL,ep=ep))
}

ep= .1
reps = 10000
al = 0.05
#nVec = c(50,100,200,400,800,1600)
nVec = c(25,50,100,200,400)
p_kw = p_mw = p_median  = p_ks = p_ksv =  rep(0,reps)
power_kw = power_mw = power_median = power_ks =power_ksv = rep(0,length(nVec))
for (k in 1:length(nVec)){
  print(k/length(nVec))
  n = nVec[k]
  null_kw = ecdf(reference(n=n,ep=ep,reps=10000,type = "kw"))
  null_MN = ecdf(reference(n=n,ep=ep,reps=10000,type = "mn"))
  null_median = ecdf(reference(n=n,ep=ep,reps=10000,type = "median"))
  null_ks = ecdf(reference(n=n,ep=ep,reps=10000,type = "ks"))
  null_ksv = ecdf(reference(n=n,ep=ep,reps=10000,type = "ks",Kuiper = "Kuiper"))
  for(i in 1:reps){
    #x = rnorm(n,mean = 0, sd= 1)
    #y = rnorm(n,mean = 1, sd= 1)
    x = rcauchy(n,1,1)
    y = rcauchy(n,0,1)
    #x = rlaplace(n,1,1)
    #y = rlaplace(n,0,1)
    #x = rlogis(n,1,1)
    #y = rlogis(n,0,1)
    # x = rexp(n,rate = 1)-log(2)/1
    # y = rexp(n,rate = 1)-log(2)/1
    kw = Kruskal_Wallis(x = x,y = y,z = NULL,ep = ep)
    p_kw[i] = 1 - null_kw(kw)
    u_t = Mann_Whitney(x = x,y =y,z=NULL,ep = ep)
    p_mw[i] = null_MN(u_t)
    median = median_test_DP(x = x,y =y,z=NULL,ep = ep)
    p_median[i] = 1-null_median(median)
    ks = KS(x = x,y = y,z = NULL,ep = ep,Kuiper = FALSE)
    p_ks[i] = 1-null_ks(ks)
    ksv = KS(x = x,y = y,z = NULL,ep = ep,Kuiper = "Kuiper")
    p_ksv[i] = 1-null_ksv(ksv)
  }
  power_kw[k] = mean(p_kw < al)
  power_mw[k] = mean(p_mw < al)
  power_median[k] = mean(p_median < al)
  power_ks[k] = mean(p_ks < al)
  power_ksv[k] = mean(p_ksv < al)
}
plot(log(nVec),power_kw,type="l",col="black",ylim=c(0,1),xaxt="n",ylab="Power",xlab="n",lty=1,lwd = 1.7)
lines(log(nVec),power_mw,col="red",lty=3,lwd = 1.7)
lines(log(nVec),power_median,col="yellow",lty=5,lwd = 1.7)
lines(log(nVec),power_ks,col="green",lty=7,lwd = 1.7)
lines(log(nVec),power_ksv,col="blue",lty=9,lwd = 1.7)
axis(at=log(nVec),labels=nVec,side=1)
legend("topleft",c("Kruskal Wallis","Mann Whitney","Median","KS","Kuiper"),col=c("black","red","yellow","green","blue"),lty=c(1,3,5,7,9),lwd = 1.7,bty="o")
bottomright
topleft

power_kw
power_mw
power_ks
power_median
power_ksv
# different epsilon.1,.5,1, different distribution, rlaplace, cauchy, logistics, normal 0,1;0,0.1;
# two sided
