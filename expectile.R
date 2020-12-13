cterSimData <- function(n, bet0, tau, modtype = 1, errtype = 1)
{
  
  ### calculate the true expectile
  pefun <- function(t, errtype){
    
    if (errtype==1){
      ## example1: normal distribution
      F <- pnorm(t, 0, 1)
      
      integrand <- function(x) {x*dnorm(x, 0, 1)}
      
    }
    G <- integrate(integrand, lower = -(1e+3), upper = t)$value
    gmean <-  integrate(integrand, lower = -(1e+3), upper = 1e+3)$value
    u <- G -t * F
    asy <- u/(2*u + t-gmean)
    
    return(asy)
  }
  
  efun <- function (tau, errtype){
    tau[tau > 1 | tau < 0] = NA
    zz = 0 * tau
    lower = rep(-10, length(tau))
    upper = rep(10, length(tau))
    diff = 1
    index = 1
    while (diff > 1e-10 && index < 1000) {
      root = pefun(zz, errtype) - tau
      root[is.na(root)] = 0
      lower[root < 0] = zz[root < 0]
      upper[root > 0] = zz[root > 0]
      zz = (upper + lower)/2
      diff = max(abs(root), na.rm = T)
      index = index + 1
    }
    zz[is.na(tau)] = NA
    return(zz)
  }
  
  
  ### parameter settings
  x <- runif(n, -2, 4) 
  X <- cbind(1, x)
  
  
  if (errtype==1){
    ## example1: normal distribution
    err <- rnorm(n, 0, 1)    # normal distribution
  } 
  err0 <- err - efun(tau, errtype)   # the tau-th expectile of err is zero
  
  if (modtype == 1){
    Y <- X %*% bet0 + err0
  }
  
  dat <- cbind(Y, X)
  
  return(dat)
}

expTLlaws <- function(X, y, tau,  max.iter, tol)
{
  # X is the model matrix
  # y is the response vector of observed proportion
  # maxIter is the maximum number of iterations
  # tol is a convergence criterion
  #X <- cbind(1, X) # add constant
  b <- bLast <- rep(0, ncol(X)) # initialize
  it <- 1 # iteration index
  while (it <= max.iter){
    
    ypred <- c(X %*% b)
    w <- as.vector(tau *(y>= ypred) + (1-tau)* (y<ypred))
    b <- lsfit(X, y, w, intercept=FALSE)$coef
    if (max(abs(b - bLast)/(abs(bLast) + 0.01*tol)) < tol) break
    bLast <- b
    it <- it + 1 # increment index
  }
  if (it > max.iter) warning('maximum iterations exceeded')
  
  ## the loss function
  loss <- sum(w*(y-c(X%*%b))^2)
  
  ## the variance
  Am <- t(X) %*% diag(w) %*% X
  if(min(eigen(Am)$values)<1e-8){Am=as.matrix(nearPD(Am)$mat)}
  
  A <- solve(Am)
  H <- w * diag(X %*% A %*% t(X))
  B <- t(X) %*% diag(w^2*(y-ypred)^2/(1-H)) %*% X
  Vb <- A %*% B %*% A  # variance function
  list(coefficients = b, variance = Vb, loss = loss, it = it)
}



## simulation
set.seed(123)
n <- 500
theta0 <- c(1, 5)  # true value
tau <- 0.7
NS <- 1000

p <- length(theta0)
theta.est <- matrix(0, NS, p)
theta.ese <- matrix(0, NS, p)
for (kk in 1:NS){
  
  ## generated data
  d <- cterSimData(n,theta0,tau = 0.7)
  
  fit <- expTLlaws(X=d[,2:3],y=d[,1],tau = 0.7,max.iter = 100,tol = 1)
  
  theta.est[kk, ] <- fit$coefficients
  
  theta.ese[kk, ] <- diag(sqrt(fit$variance))
}

### 
bias <- apply(theta.est, 2, mean) - theta0
sd <- apply(theta.est, 2, sd)
ese <- apply(theta.ese, 2, mean)
alp <- 0.05
c.alp <- qnorm(1-alp/2, 0, 1)
theta000 <- matrix( rep(theta0, NS), nrow=NS, ncol=2, byrow = TRUE)
coverage <-  (theta.est-c.alp* theta.ese <= theta000) *
  (theta000 <= theta.est + c.alp* theta.ese )
CP <- apply(coverage, 2, mean)  


out <- rbind(bias, sd, ese, CP)
colnames(out) <- c("bet0", "bet1")

out
