# APPLICATION to  REAL DATA
library(mix)
data(stlouis)
stlouis
?stlouis
y <- stlouis[,4]
x <- stlouis[,1]
data <- cbind(x,y)
index <- complete.cases(data)
data <- data[index,]
data[,1]<- 4-data[,1] # reordering the ratings
x <- data[,1]
y <- data[,2]
cor(x,y)

# scatter-plot
op<-par()
par(mai=c(0.75,0.75,0.25,0.1), mgp=c(2.5,1,0))
plot(y,x,axes=FALSE,ylab=expression(X[2*scriptstyle(D)]),xlab=expression(X[1]), pch=20)
# Calculate frequency of each pair
pairs_freq <- table(paste(y, x))
# Adjust point size based on frequency
for (pair in names(pairs_freq)) {
  pair_values <- as.numeric(strsplit(pair, " ")[[1]])
  points(pair_values[1], pair_values[2], pch=20, cex=sqrt(pairs_freq[pair]))
}
axis(2,at=c(1,2,3))
axis(1)
box()
#  --> Figure 8
par(op)
# 1.
# estimation of the copula using BiCopSelect and the ranks for both variables
aggregate(data[,2] ~ data[,1], data, mean)
cor(sort(x), sort(y))
library(VineCopula)
u <- pobs(cbind(x,y), ties.method="first")
BiCopSelect(u[,1],u[,2])
library(copula)
cop<-frankCopula(dim=2, param=2.48) # or 2.700245
cop<-frankCopula(dim=2, param=2.700245) # or 2.700245
# estimation of the margins
table(x)
hist(y)
library(moments)
skewness(y)
kurtosis(y)
library(fitdistrplus)
fit.p <- fitdist(y,"weibull")
sh <- fit.p$estimate[1] # shape parameter
sc <- fit.p$estimate[2] # scale parameters
mu <- sc*gamma(1+1/sh)  # expectation of a Weibull rv
sigma2 <- sc^2*gamma(1+2/sh)-mu^2 # variance of a Weibull rv
# constructing the bivariate katent distribution with Weibull margins and Frank copula
dist <- mvdc(cop, c("weibull","weibull"), param=list(list(scale=sc, shape=sh),list(scale=sc, shape=sh)))
library(cubature)
integrand <- function(arg)
{
  x1 <- arg[1]
  x2 <- arg[2]
  x1*x2*dMvdc(c(x1,x2), dist)
}
# computing the correlation before discretization (polyserial correlation)
res   <- cuhre(f = integrand,
               lowerLimit = c(0,0),
               upperLimit = c(Inf, Inf),
               relTol = 1e-4, absTol= 1e-12)
rho <- (res$integral-mu*mu)/sqrt(sigma2*sigma2)
rho
# correlation after discretization (point-polyserial correlation)
integrand.d <- function(arg)
{
  x1 <- arg[1]
  x2 <- arg[2]
  x1*dMvdc(c(x1,x2), dist)
}
p <- (table(x)/sum(table(x)))
P <- cumsum(p)
res.1   <- cuhre(f = integrand.d,
               lowerLimit = c(0,0),
               upperLimit = c(Inf, qweibull(P[1],sh,sc)),
               relTol = 1e-4, absTol= 1e-12)
res.2   <- cuhre(f = integrand.d,
                 lowerLimit = c(0,qweibull(P[1],sh,sc)),
                 upperLimit = c(Inf, qweibull(P[2],sh,sc)),
                 relTol = 1e-4, absTol= 1e-12)
res.3   <- cuhre(f = integrand.d,
                 lowerLimit = c(0,qweibull(P[2],sh,sc)),
                 upperLimit = c(Inf, Inf),
                 relTol = 1e-4, absTol= 1e-12)
mu.d <- sum(1:3*p)
sigma2.d <- sum((1:3)^2*p)-mu.d^2
rho.d <- (res.1$integral+2*res.2$integral+3*res.3$integral-mu*mu.d)/sqrt(sigma2*sigma2.d)
rho.d

# GRAPH

x.x <- seq(50,150,0.1)
y.y <- seq(50,150,0.1)
z <- dMvdc(cbind(x.x,y.y), dist)

op<-par()
layout(matrix(c(1, 2), ncol = 1), heights = c(2, 5))
par(mar = c(0, 4, 0, 2) + 0.1, mgp=c(1,1,0))
plot(function(x) dweibull(x, shape=sh, scale=sc),xaxt = "n",yaxt="n", xlim=c(70,140),xlab="", ylab="pdf")
par(mar = c(4, 4, 0, 2) + 0.1, mgp=c(2.5,1,0))
cp <- contour(dist, dMvdc, n.grid=36, xlim=c(70,140), ylim=c(70,140),xlab=expression(X[1]),ylab=expression(X[2]))
layout(1)
#  --> Figure 9
par(op)

# a piece of simulation for checking the results
set.seed(0123)
s <- rMvdc(100000,dist)
cor(s)
q <- cumsum(table(x)/sum(table(x)))
q <- c(1,2,3)/3
s[s[,1]<quantile(s[,1],q[1]),1] <- 1
s[s[,1]>quantile(s[,1],q[1]) & s[,1]<quantile(s[,1],q[2]),1] <- 2
s[s[,1]>quantile(s[,1],q[2]),1] <- 3
cor(s)

# new part
# 2. taking into account the discrete nature of the X variable
library(rvinecopulib)
#
empirical_cdf <- function(x)
{
  x_sorted <- sort(x)
  ranks <- rank(x_sorted, ties.method = "max")  # Compute ranks with handling ties
  cdf_values <- ranks / length(x)  # Normalize to get empirical CDF values
  return(cdf_values[order(order(x))])  # Return values in original order
}
#
empirical_cdf_left <- function(x) {
  n <- length(x)
  x_sorted <- sort(x)  # Sort data
  unique_x <- unique(x_sorted)  # Unique values
  counts <- table(x_sorted)  # Frequency table
  cum_counts <- cumsum(as.numeric(counts)) / n  # CDF values
  cdf_left <- c(0, cum_counts[-length(cum_counts)])  # Left limit F(x-)
  
  # Map each x value to its corresponding F(x-), avoiding names
  cdf_left_vector <- as.numeric(cdf_left[match(x, unique_x)])
  return(cdf_left_vector)
}
# running the function bicop in rvinecopulib
# which accommodates for ties for discrete ("d") variables
# we restrict ourselves to one-parameter copula families
data.1  <- empirical_cdf(x)      # cdf
data.11 <- empirical_cdf_left(x) # left limit of the cdf
data.2 <- pobs(y) 
data    <- cbind(data.1,data.2,data.11)
res <- bicop(data, var_types = c("d", "c"), family_set = "onepar")
summary(res)
# still a Frank copula, with a slightly different parameter value than before

# 3.
# let us follow a full-likelihood approach
# assuming Weibull margins and Frank copula
# for the UNDERLYING bivariate distribution
# for the mixed discrete-continuous sample the log-likelihood function is as follows:
logliklat <- function(par)
{
  print(par)
  qx1<-par[1] # 1st threshold
  qx2<-par[2] # 2nd threshold
  lambda<-par[3]
  kappa<-par[4]
  theta<-par[5] # copula parameter
  qx <-c(0,qx1,qx2,Inf)
  cop <- frankCopula(param=theta, dim=2)
  joint <- mvdc(copula=cop, margins="weibull", paramMargins=list(list(shape=as.numeric(kappa),scale=as.numeric(lambda))),marginsIdentical=TRUE)
  n <- length(x)
  S <- 0
  for(i in 1:n)
  {
    fz <- Vectorize(function(z) dMvdc(c(z,y[i]),joint))
    S <- S +  log(integrate(function(z) fz(z),lower=qx[x[i]],upper=qx[x[i]+1])$value)
  }
  return(-S)
}
library(Rsolnp)
# start values for the parameters
theta.1 <- qweibull(mean(x==1), fit.p$estimate[1], fit.p$estimate[2])
theta.2 <- qweibull(mean(x<=2), fit.p$estimate[1], fit.p$estimate[2])
par <- c(theta.1, theta.2, fit.p$estimate[2],fit.p$estimate[1], 2.48)
# maximizing the log-likelihood function, we obtain:
sol <- solnp(pars=par, fun=logliklat, LB=c(0,0,0,0,-Inf), UB=rep(Inf,5))
# $pars
# [1]  98.181284 114.851221 114.545798   8.251416   2.623862
round(sol$pars,2)
