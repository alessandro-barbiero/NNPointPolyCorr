# PLOTS of PEARSON'S RHO for a bivariate RV with GAUSS COPULA
# and
# NORMAL
# or
# UNIFORM
# or
# EXPONENTIAL MARGINS

n <- 1000
rho.v <- (1:n)/n
# vector or rho values between 0.001 and 1 with constant step=0.001
rho.exp  <- numeric(n)
#
for(i in 1:n)
{
  Cop   <- normalCopula(param=rho.v[i], dim=2)
  # ASSIGN THE MARGINS
  model <- mvdc(Cop,  c("exp","exp"), list(1,1))
  # CONSTRUCT the function for the MIXED MOMENT of the MODEL
  integrand0 <- function(arg)
  {
    x1 <- arg[1]
    x2 <- arg[2]
    x1*x2*dMvdc(c(x1,x2),model)
  }
  res0 <- cuhre(f = integrand0,
                lowerLimit = c(0, 0),
                upperLimit = c(Inf, Inf),
                relTol = 1e-4, absTol= 1e-12)
  rho.exp[i]  <- (res0$integral-1)
  print(i)
}
# save(list=ls(),file="rhoExpGauss.Rdata")

par(mai=c(.8,.8,.1,.1), mgp=c(2.5,1,0))
plot(function(rho) 6/pi*asin(rho/2), xlim=c(0,1), xlab=expression(rho), ylab=expression(rho[X[1]*X[2]]), lwd=1.5, col="blue")
points(rho.v, rho.exp, col="red", cex=.025, pch=12, type="p")

points(rho.v[400], rho.exp[400], col="red", pch=12)
points(rho.v[500], 6/pi*asin(rho.v[500]/2), pch=8,col="blue")
points(0.6,0.6,pch=1)

legend(x="bottomright", col=c("black","blue","red"), text.col=c("black","blue","red"), legend=c("normal", "uniform","exp"), pch=c(1,8,12))

abline(0,1,col="black")


# PLOT OF SPEARMAN'S RHO and PEARSON'S RHO as FUNCTIONS of the COPULA PARAMETER
# GIVEN THE MARGINS

# Gumbel copula and exponential margins
library(copula)
library(cubature)

theta.v <- seq(1,5,0.02)
n <- length(theta.v)
rhoS.exp <- numeric(n)
rho.exp <- numeric(n)

for(i in 1:n)
{
  Cop   <- gumbelCopula(param=theta.v[i], dim=2)
  rhoS.exp[i] <- rho(Cop)
  # ASSIGN THE MARGINS
  model <- mvdc(Cop,  c("exp","exp"), list(1,1))
  # CONSTRUCT the function for the MIXED MOMENT of the MODEL
    integrand0 <- function(arg)
    {
      x1 <- arg[1]
      x2 <- arg[2]
      x1*x2*dMvdc(c(x1,x2),model)
    }
    res0 <- cuhre(f = integrand0,
                  lowerLimit = c(0, 0),
                  upperLimit = c(Inf, Inf),
                  relTol = 1e-4, absTol= 1e-12)
   rho.exp[i]  <- (res0$integral-1)
   print(i)
}

par(mai=c(.8,.8,.1,.1), mgp=c(2.5,1,0),cex=1.15)

plot(theta.v, rho.exp, xlab=expression(theta), ylab=expression(paste(rho[S],",",rho)), lwd=1.5, col="blue",xlim=c(1,5),ylim=c(0,1),,pch=0)
points(theta.v, rhoS.exp, col="red", pch=10, type="b")

legend(x="bottomright", col=c("blue","red"), text.col=c("blue","red"), legend=c("Pearson's rho", "Spearman's rho"), pch=c(0,10))

# save(list=ls(),file="rhorhoSGumbel.Rdata")


# Clayton copula and exponential margins
library(copula)
library(cubature)

theta.v <- seq(0,10,0.02)
n <- length(theta.v)
rhoS.exp <- numeric(n)
rho.exp <- numeric(n)

for(i in 1:n)
{
  Cop   <- claytonCopula(param=theta.v[i], dim=2)
  rhoS.exp[i] <- rho(Cop)
  # ASSIGN THE MARGINS
  model <- mvdc(Cop,  c("exp","exp"), list(1,1))
  # CONSTRUCT the function for the MIXED MOMENT of the MODEL
  integrand0 <- function(arg)
  {
    x1 <- arg[1]
    x2 <- arg[2]
    x1*x2*dMvdc(c(x1,x2),model)
  }
  res0 <- cuhre(f = integrand0,
                lowerLimit = c(0, 0),
                upperLimit = c(Inf, Inf),
                relTol = 1e-4, absTol= 1e-12)
  rho.exp[i]  <- (res0$integral-1)
  print(i)
}

par(mai=c(.8,.8,.1,.1), mgp=c(2.5,1,0),cex=1.15)

plot(theta.v, rho.exp, xlab=expression(theta), ylab=expression(paste(rho[S],",",rho)), lwd=1.5, col="blue",xlim=c(0,10),ylim=c(0,1),pch=0)
points(theta.v, rhoS.exp, col="red", pch=10, type="b")

legend(x="bottomright", col=c("blue","red"), text.col=c("blue","red"), legend=c("Pearson's rho", "Spearman's rho"), pch=c(0,10))

# save(list=ls(),file="rhorhoSClayton.Rdata")


