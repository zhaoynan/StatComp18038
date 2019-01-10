## ----eval=FALSE----------------------------------------------------------
#  Cau <- function(theta, eta, y)
#  {
#    .dcauchy <- function(theta, eta, x)
#    {
#      #compute the density function of cauchy distribution
#      if(theta <= 0) return(NA)
#  
#      r <- theta * pi * (1 + ((x-eta)/theta)^2)
#      return(1 / r)
#    }
#  
#    n <- length(y)
#    res <- numeric(n)
#    for(i in 1:n)
#    {
#      res[i] <- integrate(.dcauchy, lower = -Inf, upper = y[i],
#                          rel.tol = .Machine$double.eps^0.25,
#                          theta = theta, eta = eta)$value
#    }
#    return(res)
#  }
#  

## ------------------------------------------------------------------------
library(StatComp18038)
y <- seq(-10, 10, 1)
theta <- 2
eta <- 3
res1 <- Cau(theta, eta, y)
res2 <- pcauchy(y, eta, theta)
round(rbind(res1, res2), 5)

## ------------------------------------------------------------------------
MCe <- function(x, n)
{
  y <- runif(n)
  cdf <- numeric(length(x))

  for (i in 1:length(x))
  {
      g <- (x[i])^3*(y*(1-y*x[i]))^2
      cdf[i] <- mean(g)/beta(3,3)
  }
   return(cdf)
}

## ------------------------------------------------------------------------
library(StatComp18038)
x <- seq(0.1, 0.9, 0.1)
n <- 10000
res1 <- MCe(x, n)
res2 <- pbeta(x,3,3)
round(rbind(res1, res2), 3)

## ----eval=FALSE----------------------------------------------------------
#  MC.Phi <- function(x, sigma, r, antithetic = TRUE)
#  {
#    y <- runif(r)
#    if(!antithetic) v <- runif(r)
#    else
#      v <- 1-y
#    y <- c(y, v)
#    cdf <- numeric(length(x))
#    for(i in 1:length(x))
#    {
#      g <- (y*x[i]^2)*exp(-(y*x[i]^2)/(2*(sigma^2)))
#      cdf[i] <-  1/(sigma^2)*mean(g)
#    }
#    cdf
#  }

## ----eval=TRUE-----------------------------------------------------------
library(StatComp18038)
x <- seq(0.1, 0.9, 0.1)
sigma <- 2
r <- 10000
res1 <- MC.Phi(x, sigma, r, antithetic = TRUE)
res2 <- 1-exp(-((x)^2)/(2*(sigma^2))) 
round(rbind(res1, res2), 3)

## ------------------------------------------------------------------------
ID <- c(1, 2, 3, 4)
age <- c(24, 25, 50, 30)
diabetes <- c("type1", "type2", "type1", "type1")
status <- c("poor", "improved", "excellent", "poor")
diabetes <- factor(diabetes)
status <- factor(status, order = TRUE)
pdata <- data.frame(ID, age, diabetes, status)
pdata
table(pdata$diabetes, pdata$status)
summary(pdata)

## ------------------------------------------------------------------------
x <- c(1:10)
y <- x
z <- 10/x
opar <- par(no.readonly = TRUE)
par(mar = c(5, 4, 4, 8)+0.1)   
plot(x, y, type = "b", pch = 22, col = "red", 
     yaxt = "n", lty = 3, ann = FALSE)
lines(x, z, type = "b", pch = 25, col = "blue", lty = 2) 
axis(side = 2, at = x, labels = x, col.axis = "red", las = 2) 
axis(side = 4, at = z, labels = round(z, digits = 2),
     col.axis = "blue", las = 2, cex.axis = 0.7, tck = -0.1)
mtext("y=1/x", side = 4, line = 3, cex.lab = 1, las = 2, col= "blue") 
title(main = "An example of creative axes", col.main = "lightgreen", xlab = "X values", ylab = "Y=X")

## ------------------------------------------------------------------------
x <- c(0, 1, 2, 3, 4)    # the number of xi for i=1,2,3,4,5
p <- c(0.1, 0.2, 0.2, 0.2, 0.3)   # the probability of each xi 
cp <- cumsum(p)  # caculate the cumulative distribution function 
m <- 1000     
r <- numeric(m)
set.seed(1234)  
r <- x[findInterval(runif(m), cp)+1]
table(r)
ct <- as.vector(table(r))
q <- ct/sum(ct)/p  #the specific values
cbind(x = x, p = p, freq = q)

## ------------------------------------------------------------------------
set.seed(1234)
r1 <- sample(c(0,1,2,3,4), m, replace = T, c(0.1,rep(0.2,3),0.3))
ct1 <- as.vector(table(r1))
q1 <- ct1/sum(ct1)/p
cbind(x = x, p = p, freq = q1)

## ------------------------------------------------------------------------
n <- 1000
j <- 0  # counter for accepted
k <- 0  #  iterations
y <- numeric(n)
while(k < n)
{
u <- runif(1)
j <- j+1
x <- runif(1)   # random variable form g(x)
if(4*x^2*(1-x) > u)
{  # we accept x
  k <- k+1
  y[k] <- x
}
}
j 
hist(y, prob = TRUE, main = expression(f(x)==12*x^2*(1-x)))
p <- seq(0, 1, 0.01)
lines(p, 12*p^2*(1-p))

## ------------------------------------------------------------------------
n <- 1000
r <- 4
beta <- 2
lamba <- rgamma(n, r, beta)         # lamba is generate randomly
y <- rexp(n, lamba)                 # generate 1000 random observations
plot(y, ylim = c(0,6), main = "The random observations")


## ------------------------------------------------------------------------
x <- seq(0.1, 0.9, length = 9)
n <- 10000
y <- runif(n)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
  g <- (x[i])^3*(y*(1-y*x[i]))^2
  cdf[i] <- mean(g)/beta(3,3)
}
phi <- pbeta(x,3,3) # gengerate the values by the pbeta function
print(round(rbind(x, cdf, phi), 3)) # Compare the estimates with the values returned by the pbeta function


## ------------------------------------------------------------------------
MC.Phi <- function(x, r = 10000, antithetic = TRUE){
  y <- runif(r/2)
  if(!antithetic) v <- runif(r/2)
  else
    v <- 1-y
  y <- c(y, v)
  cdf <- numeric(length(x))
  for(i in 1:length(x)){
    g <- (y*x[i]^2)*exp(-(y*x[i]^2)/(2*(2^2)))
    cdf[i] <-  1/(2^2)*mean(g) 
  }
  cdf
} # use antithetic varibules to implement a function to generate samples from a Rayleigh(2) distribution

x <- seq(0.1, 2.5, length = 10)
set.seed(12345)
Phi <- 1-exp(-((x)^2)/(2*(2^2)))  # Results generated from Rayleigh(2) 
MC1 <- MC.Phi(x, antithetic = FALSE)
MC2 <- MC.Phi(x)
print(round(rbind(x,Phi, MC1, MC2), 5)) # A comparison of estimates obtained from a single Monte Carlo experiment

## ------------------------------------------------------------------------
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.65
for (i in 1:m) {
  MC1[i] <- MC.Phi(x, r = 1000, antithetic = FALSE)
  MC2[i] <- MC.Phi(x, r = 1000)
}
print(sd(MC1))
print(sd(MC2))
print((var(MC1)-var(MC2))/var(MC1))

## ------------------------------------------------------------------------
x <- seq(1, 10, 0.1)
g <- (1/(sqrt(2*pi)))*x^2 * exp(-x^2/2)
f1 <- dnorm(x, 2, 1.1) 
f2 <- dgamma(x, 3, 2)
# importance functions: f1,f2 with g(x) and ratios g(x)/f(x) 
# figure 1
plot(x, g, type = "l", main = "", ylab = "", ylim = c(0, 0.4), lwd = 2)
lines(x, f1, lty = 2, lwd = 2, col = "red")
lines(x, f2, lty = 3, lwd = 2, col = "purple")
legend("topright", legend = c("g", "f1", "f2"), lty = 1:3, col = c("black", "red", "purple"), lwd = 2, inset = 0.02)
# figure 2
plot(x, g, type = "l", main = "", ylab = "", ylim = c(0, 3), lwd = 2)
lines(x, g/f1, lty = 2, lwd = 2, col = "red")
lines(x, g/f2, lty = 3, lwd = 2, col = "purple")
legend("topright", legend = c("g", "f1", "f2"), lty = 1:3, col = c("black", "red", "purple"), lwd = 2, inset = 0.02)

m <-10000 
se <- numeric(2)
g <- function(x){
  (1/(sqrt(2*pi)))*x^2 * exp(-x^2/2) * (x > 1)
}
set.seed(123)
x <- rnorm(m, 2, 1.1)
fg <- g(x) / dnorm(x, 2, 1.1)
se[1] <- sd(fg)

x <- rgamma(x, 3,2)
fg <- g(x) /  dgamma(x, 3, 2)
se[2] <- sd(fg)
print(round(se, 5)) #the corresponding standard errors

## ------------------------------------------------------------------------
m <-10000 
theta.hat  <- numeric(2)
g <- function(x){
  (1/(sqrt(2*pi)))*x^2 * exp(-x^2/2) * (x > 1)
}
set.seed(123)
x <- rnorm(m, 2, 1.1)
fg <- g(x) / dnorm(x, 2, 1.1)
theta.hat[1] <- mean(fg)

x <- rgamma(x, 3,2)
fg <- g(x) /  dgamma(x, 3, 2)
theta.hat[2] <- mean(fg)
phi <- (1/sqrt(2*pi))*exp(-1/2)+1-pnorm(1) # the real value of the function 
print(round(c(theta.hat, phi), 5))

## ------------------------------------------------------------------------
f <- function( n, x, mu){
  y <- numeric(n)
  for(i in 1:n){
    y[i] <- (2*i-n-1) * x[i]
  }
  return(1/(n^2*mu)*sum(y))
}
set.seed(123)
#the mean,median and deciles of G.hat with X is standard lognormal
ucl1 <- replicate(1000, expr = f( n <- 20, x <- sort(rlnorm(n)), mu <- mean(rlnorm(n))))
print(mean(ucl1))
print(median(ucl1))
quantile(ucl1, seq(0,1,0.1))
hist(ucl1,probability = TRUE, main = "hist of standard lognormal") # histogram of the replicates with x is standard lognormal

## ------------------------------------------------------------------------
#the mean,median and deciles of G.hat with X is uniform distribution
ucl2 <- replicate(1000, expr = f( n <- 20, x <- sort(runif(n)),mu <- mean(runif(n))))
print(mean(ucl2))
print(median(ucl2))
quantile(ucl2, seq(0,1,0.1))
hist(ucl2, probability = TRUE,main = "hist of uniform distribution") # histogram of the replicates with x is uniform distribution

## ------------------------------------------------------------------------
#the mean,median and deciles of G.hat with X is Bernoulli(0.1)
ucl3 <- replicate(1000, expr = f( n <- 20, x <- sort(rbinom(n,20,prob = 0.1)),mu <- mean(rbinom(n,20,prob = 0.1))))
print(mean(ucl3))
print(median(ucl3))
quantile(ucl3, seq(0,1,0.1))
hist(ucl3, probability = TRUE,main = "hist of Bernoulli") # histogram of the replicates with x is Bernoulli(0.1)

## ------------------------------------------------------------------------
sigma <- 1
f <- function(n, y){# function to construct confident interval for sigma 
  s2 <- (1/(n-1))*sum((y-mean(y))^2)
  c((n-1)*s2/qchisq(0.975, n-1), (n-1)*s2/qchisq(0.025, n-1))
}
set.seed(123)
result <- replicate(1000, f(n <- 20, y <- log(rlnorm(n))))
r.true <- 2*pnorm(sigma/sqrt(2))-1 # the true value of r
r.hat <- 2*pnorm(sqrt(result/2))-1
c.r <- cbind(r.hat[1,], r.hat[2,])
judge <- numeric(1000)
for(i in 1:1000)
{
  judge[i] <- I(c.r[i,1] <= r.true & c.r[i,2] >= r.true)
}
sum(judge)
mean(judge)


## ------------------------------------------------------------------------
library(MASS)
n <- 20
m <- 1000
c <- c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))
l <- length(c)
power.p <- numeric(l)
library(MASS)
set.seed(123)
for(i in 1:l)
{
  p.p <- replicate(m, expr = {
    x <- mvrnorm(n, c(0,1), matrix(c(1, c[i], c[i], 1), 2, 2))
     cortest <- cor.test(x[,1], x[,2], alternative = "two.sided", method = "pearson")
         cortest$p.value } )
  power.p[i] <- mean(  p.p <= 0.05 )
}
power.p 

## ------------------------------------------------------------------------
n <- 20
m <- 1000
c <- c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))
l <- length(c)
power.s <- numeric(l)
library(MASS)
set.seed(123)
for(i in 1:l)
{
  p.s <- replicate(m, expr = {
    x <- mvrnorm(n, c(0,1), matrix(c(1, c[i], c[i], 1), 2, 2))
    cortest <- cor.test(x[,1], x[,2], alternative = "two.sided", method = "spearman")
    cortest$p.value } )
  power.s[i] <- mean(  p.s <= 0.05 )
}
power.s

## ------------------------------------------------------------------------
n <- 20
m <- 1000
c <- c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))
l <- length(c)
power.k <- numeric(l)
library(MASS)
set.seed(123)
for(i in 1:l)
{
  p.k <- replicate(m, expr = {
    x <- mvrnorm(n, c(0,1), matrix(c(1, c[i], c[i], 1), 2, 2))
    cortest <- cor.test(x[,1], x[,2], alternative = "two.sided", method = "kendall")
    cortest$p.value } )
  power.k[i] <- mean(  p.k <= 0.05 )
}
power.k

## ------------------------------------------------------------------------
power <- cbind(power.p, power.s, power.k)
rownames(power) <- c(seq(-0.9, -0.1, 0.1), seq(0.1, 0.9, 0.1))
power

## ------------------------------------------------------------------------
n <- 20
m <- 1000
set.seed(123)
p.p1 <- replicate(m, expr = {
  x <- rexp(n, 4)
  y <- runif(n, 0, 1)
  cortest <- cor.test(x, y, alternative = "two.sided", method = "pearson")
  cortest$p.value } )
power.p1 <- mean(  p.p1 <= 0.05 )
power.p1

## ------------------------------------------------------------------------
n <- 20
m <- 1000
set.seed(123)
p.s1 <- replicate(m, expr = {
  x <- rexp(n, 4)
  y <- runif(n, 0, 1)
  cortest <- cor.test(x, y, alternative = "two.sided", method = "spearman")
  cortest$p.value } )
power.s1 <- mean(  p.s1 <= 0.05 )
power.s1

## ------------------------------------------------------------------------
library(bootstrap)
n <- nrow(law)
LSAT <- law$LSAT
GPA <- law$GPA
theta.hat <- cor(LSAT, GPA)

# compute the jackknife replictes,leave-one-out estimates
theta.jack <- numeric(n)
for(i in 1:n)
  theta.jack[i] <- cor(LSAT[-i], GPA[-i]) 
bias <- (n-1) * (mean(theta.jack) - theta.hat)   #jackknife estimate of bias
se <- sqrt((n-1) * mean((theta.jack- mean(theta.jack))^2))   #jackknife estimate of standard eror
  print(round(c(bias, se), 5)) 

## ------------------------------------------------------------------------
library(boot)
  data("aircondit", package = "boot")

    # make the set of aircondit a vector
  aircondit <- unlist(aircondit)
  aircondit <- as.vector(aircondit)
  lambda <- mean(aircondit)
  
  lambda.boot <- function(x,i)
  {
    # function to compute the statistic
    mean(x[i])
  }
  x <- aircondit
boot.obj <- boot(x, statistic = lambda.boot, R = 2000 )
print(boot.ci(boot.obj, type = c("norm", "basic", "perc","bca"))) 

## ------------------------------------------------------------------------

data(scor, package = "bootstrap")

# turn scor into a matrix x
names(scor) <- NULL
x <- as.matrix(scor)

n <- nrow(x)
m <- ncol(x)
x.bar <- numeric(m)

#compute the theta
theta <- function(cov)
{
  eigen.cov <- eigen(cov)   
  eigen.v <- sort(eigen.cov$values, decreasing = TRUE)
  eigen.v[1]/sum(eigen.v)
}

# compute the MLE of covariance matrix
cov.hat <- function(j)
{
  x.bar <- colMeans(j)
  I <- matrix(rep(1, nrow(j)), nrow = nrow(j), ncol = 1)
  x.bar <- I %*% x.bar
  crossprod(j - x.bar, j - x.bar)/nrow(j)
}
theta.hat <- theta(cov.hat(x))

#  the jackknife estimates of bias and standard error 
cov.jack <- matrix(0, nrow = m, ncol = m)
theta.jack <- numeric(n)
for(i in 1:n)
{
  cov.jack <- cov.hat(x[-i, ])
  theta.jack[i] <- theta(cov.jack)
}

bias <- (n-1) * (mean(theta.jack) - theta.hat)   #jackknife estimate of bias
se <- sqrt((n-1) * mean((theta.jack- mean(theta.jack))^2))   #jackknife estimate of standard eror
print(round(c(bias, se), 5)) 


## ------------------------------------------------------------------------
library(lattice)
library(DAAG)
attach(ironslag)
n <- length(magnetic) 

# the function of compute the change of R-square
shrinkage <- function(fit, k= n/2)
{
  require(bootstrap)
  
  theta.fit <- function(x,y)
  {
    lsfit(x,y)
  }
  
  theta.predict <- function(fit,x)
  {
    cbind(1,x) %*% fit$coef
  }
  
  x <- fit$model[,2:ncol(fit$model)]
  y <- fit$model[,1]
  
  results <- crossval(x, y, theta.fit, theta.predict, ngroup = k)
  
  r2 <- cor(y, fit$fitted.values)^2
  r2cv <- cor(y, results$cv.fit)^2
  return(r2-r2cv)   # the change of R-square
}

fit1 <- lm(magnetic ~ chemical)

fit2 <- lm(magnetic ~ chemical + I(chemical ^ 2))

fit3 <- lm(log(magnetic) ~ chemical)

fit4 <- lm(log(magnetic) ~ log(chemical))

print(c(shrinkage(fit1), shrinkage(fit2), shrinkage(fit3), shrinkage(fit4)))

## ------------------------------------------------------------------------
w2 <- function(x, y)
{
  n <- length(x)
  m <- length(y)
  
  # compute the ecdf of x and y 
  Fn <- ecdf(x)
  Gm <- ecdf(y)
  
  d <- (m * n) / (m + n)^2
  m1 <- sum((Fn(x) - Gm(x))^2)
  m2 <- sum((Fn(y) - Gm(y))^2)
  w2 <- d * (m1 + m2)
  return(w2)
}

attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

R = 999
n <- length(x)
m <- length(y)
K <- 1:(m + n)
z <- c(x, y)
w0 <- w2(x, y)
w <- numeric(R)
set.seed(12345)

for(i in 1:R)
{
  k <- sample(K, size = n, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  w[i] <- w2(x1, y1)
}
p <- mean(c(w0, w) >= w0)
p

hist(w, main = "", breaks = "scott", freq = FALSE,
     xlab = "W2(p = 0.407)")
points(w0, 0, cex = 1, pch = 16)

## ------------------------------------------------------------------------
f <- function(x)
{
  1 / (pi * (1 + x^2))
}

m <- 1e4
sigma <- 1 # Set standard deviation of normal proposal density
x <- numeric(m)
x[1] <- rnorm(1, mean = 0, sd = sigma)
k <- 0
u <- runif(m)
for(i in 2:m)
{
  xt <- x[i-1]
  y <- rnorm(1, mean = xt, sd = sigma)
  num <- f(y) * dnorm(xt, mean = y, sd = sigma)
  den <- f(xt) * dnorm(y, mean = xt, sd = sigma)
  
  if(u[i] <= num / den)
    x[i] <- y
  else
  {
    x[i] <- xt
    k <- k + 1
  }
}
print(k)


## ------------------------------------------------------------------------
index <- 5001:5500
y1 <- x[index]
plot(index, y1, type = "l", main = "", ylab = "X")

## ------------------------------------------------------------------------
b <- 1001 # discard the burnin sample
y <- x[b:m]
#a <- ppoints(100)
a <- seq(0.1, 0.9, 0.1)
QR <- qcauchy(a)
Q <- quantile(y, a)
print(round(cbind(QR, Q), 3))


## ------------------------------------------------------------------------
k <- 100
a <- seq(0, sqrt(k), 0.1)
y <- numeric(length(a))
f <- function(a)
{
  pt(sqrt(k * a^2 / (k+1-a^2)),k ) -
      pt(sqrt((k-1) * a^2 / (k-a^2)), k-1)
}

for(i in 1:length(a))
{
  y[i] <- f(a[i])
}
print(y)
plot(a, y, type = "l", xlab = "a", ylab = "y", main = "f(a)")
abline(h = 0, lty = 2, col = "red")

## ------------------------------------------------------------------------
k <- c(4:25, 100, 500, 1000)
n <-length(k)

out <- numeric(n)

for(i in 1:n)
{
  
  f <- function(a)
  {
    
    pt(sqrt(k[i] * a^2 / (k[i]+1-a^2)), k[i]) -
      pt(sqrt((k[i]-1) * a^2 / (k[i]-a^2)), k[i]-1)
  }
 
  out[i] <- uniroot(f, c(1, 2))$root
}
print(out)

## ------------------------------------------------------------------------
.dcauchy <- function(theta, eta, x)
  {
    #compute the density function of cauchy distribution
    if(theta <= 0) return(NA)
    
    r <- theta * pi * (1 + ((x-eta)/theta)^2)
    return(1 / r)
  }

.pcauchy <- function(theta, eta, y)
{
   integrate(.dcauchy, lower = -Inf, upper = y,
                   rel.tol = .Machine$double.eps^0.25,
                   theta = theta, eta = eta)$value
}

y <- matrix(seq(-10, 10, 1))

res11 <- apply(y, 1, .pcauchy, eta = 0, theta = 0.5)
res12 <-  pcauchy(y,  0, 0.5)



plot(y, res11, type = "l", col = "red", lwd=2, ylab = "cdf", ylim = c(0, 1))
lines(y, res12, type = "p")



## ------------------------------------------------------------------------
nA. = 28
nB. = 24
nOO = 41
nAB = 70

m <- 10 #max.number of interations
theta0 <- c(0.3, 0.2)  #initial est. for p,q,r
theta <- matrix(0, m, 3)
theta[1,1] <- theta0[1]
theta[1,2] <- theta0[2]
L <- numeric(m)
#tol <- .Machine$double.eps^0.5
for(i in 1:(m-1))
{
  ll <- function(theta1)
  {
    p.hat <- theta[i,1]
    q.hat <- theta[i,2]
    p <- theta1[1]
    q <- theta1[2]
    r <- 1-p-q
    r.hat <- 1-p.hat-q.hat
    f1 <- (nA.* 2*(p.hat + r.hat)/(p.hat + 2*r.hat) + nAB)*log(p)
    f2 <- (nB.* 2*(q.hat + r.hat)/(q.hat + 2*r.hat) + nAB)*log(q)
    f3 <- 2*(nOO + nA.* r.hat/(p.hat + 2*r.hat)+ nB.* r.hat/(q.hat + 2*r.hat))*log(r)
    return(- (f1+f2+f3))
  }
  theta[i+1, 1] <- optim(c(0.4, 0.2), ll)$par[1]
  theta[i+1, 2] <- optim(c(0.4, 0.2), ll)$par[2]
  theta[i+1, 3] <- 1-theta[i+1, 2]-theta[i+1, 1]
  L[i] <- ll(c(0.4, 0.2))
}
print(colMeans(theta))
x <- seq(1, m-1 , 1)
print(L[-m])
plot(x, L[-m])


## ------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# with loop
lm.mtcars <- vector("list", length(formulas))

for(i in  seq_along(formulas))
{
  lm.mtcars[[i]] <- lm(formulas[[i]], data = mtcars)$coefficients
}
lm.mtcars

# with lapply()
lapply(formulas, function(x) lm(x, data = mtcars)$coefficients)



## ------------------------------------------------------------------------
# with loop
set.seed(123)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

 # with loop
cor.mtcars <- vector("list", length(bootstraps))

for(i in seq_along(bootstraps))
{
  cor.mtcars[[i]] <- lm(mpg ~ disp, 
                        data = bootstraps[[i]])$coefficients  # coefficient
}

# with lapply()
y <- lapply(bootstraps, function(x) lm(mpg ~ disp, data = x)$coefficients)

cor.mtcars <- as.matrix(unlist(cor.mtcars))
y <- as.matrix(unlist(y))
z <- cbind(cor.mtcars, unname(y))  # turn to matrix
z

## ------------------------------------------------------------------------
# R2 with question 3
rsq <- function(mod) summary(mod)$r.squared
r1 <- numeric(length(formulas))
for(i in  seq_along(formulas))
{
  lm.mtcars[[i]] <- lm(formulas[[i]], data = mtcars)
  r1[i] <- rsq(lm.mtcars[[i]])
}
r1

# R2 with question 4
r2 <- numeric(length(bootstraps))


# with loop
cor.mtcars <- vector("list", length(bootstraps))

for(i in seq_along(bootstraps))
{
  cor.mtcars[[i]] <- lm(mpg ~ disp, 
                        data = bootstraps[[i]])# coefficient
  r2[i] <- rsq(cor.mtcars[[i]])
}
r2


## ------------------------------------------------------------------------
set.seed(123)
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
p.value <- function(test) test$p.value

round(sapply(trials, p.value), 5)

# extra challenge
x <- vector("list", length(trials))
for(i in seq_along(trials))
{
  x[i] <- trials[[i]]$p.value
}
round(unlist(x), 5)

## ------------------------------------------------------------------------
set.seed(123)
df <- data.frame(replicate(5, runif(5)))
summary <- function(x)
{
  funs <- c(mean, median, sd, mad, IQR)
  sapply(funs, function(f) f(x))
}
x <- vapply(df, summary, c(mean = 0, median = 0, sd = 0, mad = 0, IQR = 0))
x

## ------------------------------------------------------------------------
chisq.test1 <- function (x)
{
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  
    if (length(x[1,]) != length(x[2,])) 
      stop("'x'  must have the same length")
  if ((n <- sum(x)) == 0) 
    stop("at least one entry of 'x' must be positive")
  
  if (is.matrix(x)) {
    METHOD <- "Pearson's Chi-squared test"
    nr <- as.integer(nrow(x))
    nc <- as.integer(ncol(x))
    
    sr <- rowSums(x)
    sc <- colSums(x)
    E <- outer(sr, sc, "*")
    statistic <- (sum((x*x)/E)-1) * n
  }
  return(statistic)
}
M <- rbind(c(762, 327, 468, 1000, 956, 859), c(484, 239, 477, 560, 458, 589))

library(microbenchmark)
microbenchmark(t <- chisq.test1(M), Xsq <- chisq.test(M))
t
Xsq


## ------------------------------------------------------------------------
table1 <- function(x)
{
  
   z <- array(tabulate(x))
  z
}
x <- 1:10
microbenchmark(t <- table1(x), xsq <- table(x))
t
xsq

