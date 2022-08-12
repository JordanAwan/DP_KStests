# CDF function for tulap
#   t - domain of CDF
#   median - median
#   lambda - Laplace parameter, as seen on Wikipedia
#   lcut - cut the leftmost lcut amount
#   rcut - cut the rightmost rcut amount
ptulap <- function (t, median = 0, lambda = 0, cut=0) {
    lcut=cut/2
    rcut=cut/2
  # Normalize
  t <- t - median

  # Split the positive and negsative t calculations, and factor out stuff
  r <- round(t)
  g <- -log(lambda)
  l <- log(1 + lambda)
  k <- 1 - lambda
  negs <- exp((r * g) - l + log(lambda + ((t - r + (1/2)) * k)))
  poss <- 1 - exp((r * (-g)) - l + log(lambda + ((r - t + (1/2)) * k)))

  # Check for infinities
  negs[is.infinite(negs)] <- 0
  poss[is.infinite(poss)] <- 0

  # Truncate wrt the indicator on t's positivity
  is.leq0 <- t <= 0
  trunc <- (is.leq0 * negs) + ((1 - is.leq0) * poss)

  # Handle the cut adjustment and scaling
  cut <- lcut + rcut
  is.mid <- (lcut <= trunc) & (trunc <= (1 - rcut))
  is.rhs <- (1 - rcut) < trunc
  return (((trunc - lcut) / (1 - cut)) * is.mid + is.rhs)
}



#####
# Sample from tulap, where m is median and lambda is Laplace parameter.
#   n - number of samples
#   median - median
#   lambda - Laplace parameter (cf Wikipedia)
#   lcut - cut the LHS q/2
#   rcut - cut the RHS q/2
rtulap <- function (n, median = 0, lambda = 0, cut) {
    if(cut>=0){
        alpha=.95
        lcut=cut/2
        rcut=cut/2

                                        # Approximates M such that given n Bernoulli trials with success rate prob,
                                        # is such that alpha of the times, there are at least n successes among M.
                                        #   n - the number of trials that we want more than success of
                                        #   prob - Bernoulli probability success of each iid trial
                                        #   alpha - probability that among the M trials, at least n successes.
        approx.trials <- function (n, prob = 1, alpha = 0) {
                                        # Solve a quadratic form for this:
            a <- prob^2
            b <- -((2 * n * prob) + ((qnorm(alpha)^2) * prob * (1 - prob)))
            c <- n^2
            return ((-b + sqrt(b^2 - (4 * a * c))) / (2 * a))
        }

                                        # Calculate actual amount needed
        cut <- lcut + rcut
        n2 <- approx.trials(n, prob=(1 - cut), alpha=alpha)

                                        # Sample from the original Tulambda distribution
        geos1 <- rgeom(n2, prob=(1 - lambda))
        geos2 <- rgeom(n2, prob=(1 - lambda))
        unifs <- runif(n2, min=(-1/2), max=(1/2))
        samples <- median + geos1 - geos2 + unifs

                                        # Cut the tails based on the untampered CDF (ie no cuts)
        probs <- ptulap(samples, median=median, lambda=lambda)
        is.mid <- (lcut <= probs) & (probs <= (1 - rcut))
        
                                        # Abuse the NA property of R wrt arithmetics
        mids <- samples[is.mid]
        while ({len <- length(mids); len} < n) {
            diff <- n - len
            mids <- c(mids, rtulap(diff, median=median, lambda=lambda,cut=cut))
        }
        return (mids[1:n])
    }
    geos1 <- rgeom(n2, prob=(1 - lambda))
    geos2 <- rgeom(n2, prob=(1 - lambda))
    unifs <- runif(n2, min=(-1/2), max=(1/2))
    samples <- median + geos1 - geos2 + unifs
    return(samples)
    
}

#####
# Test the tulap implementations.
# This function generates plots.
tulap.test <- function () {
  start.time <- Sys.time()

  n <- 1000000
  xmin <- -10
  xmax <- 10
  x <- runif(n, min=-10, max=10)
  e <- 0.1

  res1 <- rtulap(n, lambda=exp(-e), cut=.1)
  print(length(res1))
  print(var(res1))
  print(2 / e)
  hist(res1, breaks=1e3, xlim=c(xmin, xmax))

  end.time <- Sys.time()
  print(paste("runtime (secs):",
              toString(as.numeric(difftime(end.time, start.time)))))
}

tulap.test2 <- function () {
  values = seq(-10, 10, by=0.01)
  e <- 1
  res <- ptulap(values, lambda = exp(-e),cut=.1)
  plot(res, type="l")
}

# tulap.test()
