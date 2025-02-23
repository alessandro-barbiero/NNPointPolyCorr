# the main code for computing the ratio between point-polyserial and polyserial correlations
# for several combinations of identical margins and copula
rm(list=ls())
library(copula)  # for constructing bivariate distributions combining margins and copula
library(cubature)# for computing the bivariate integral involved in the mixed moment definition
library(VGAM)
library(copBasic)# for calibrating the copula parameter
###########################
dependence <- "gauss"
#######################################
# mean and variance of the WEIBULL rv #
#######################################
Ew <- function(scale=2, shape=1.5) scale*gamma(1+1/shape)
Vw <- function(scale=2, shape=1.5) scale^2*(gamma(1+2/shape)-(gamma(1+1/shape))^2)
########################
# utility function for switching between different copulas
f.cop <- function(copula)
{
  switch(copula,
         "gumbel"  = gumbelCopula(dim=2),
         "clayton" = claytonCopula(dim=2),
         "frank"   = frankCopula(dim=2),
         "gauss"   = normalCopula(dim=2)
         )
}
#######################
# utility function for switching between different calibrating functions
# one for each copula examined
find.theta.cop <- function(copula, rho)
{
  switch(copula,
         "gumbel"  = find.theta.gumbel(rho),
         "clayton" = find.theta.clayton(rho),
         "frank"   = find.theta.frank(rho),
         "gauss"   = find.theta.gauss(rho),
  )
}
#######################
# *calibrating functions*
# for each copula, these functions calculate
# the value of the copula parameter
# inducing the assigned Spearman's correlation rho
# between two continuous margins

find.theta.clayton <- function(rho)
{
  if(rho>.9989851) theta<-78*(1-.9989851)/(1-rho) # added to handle high values of correlation
  else
  {
    res<-uniroot(function(theta) rhoCOP(CLcop,theta)-rho, interval=c(0,100))
    theta<-res$root
  }
  return(theta)
}
find.theta.gumbel <- function(rho)
{
    uniroot(function(theta) rhoCOP(GHcop,theta)-rho, interval=c(1,20), extendInt="yes")$root
}
find.theta.frank <- function(rho)
{
  iRho(frankCopula(dim=2), rho)
}
find.theta.gauss <- function(rho)
{
  iRho(normalCopula(dim=2), rho)
}
########################
# FUNCTION FIND.THETA  #
# EXPONENTIAL or ..... #
########################
# This is the central function,
# which calculates the value of the copula parameter
# inducing a target Pearson's correlation rho,
# given the margins
# It incorporates an updating rule for 
# Spearman's correlation,
# which is used as a proxy
# for Pearson's correlation
find.theta <- function(rho.i, copula, margins, par)
{
  if(copula=="clayton")   {
    cop <- claytonCopula
  }  else if(copula=="gumbel")    {
    cop <- gumbelCopula
  }  else if(copula=="frank")  {
    cop <- frankCopula
  }  else if(copula=="gauss") {
    cop <- normalCopula
  }
  rhot <- rho.i
  rhos <- rho.i
  rho  <- rho.i
  iter <- 0
  t1   <- Sys.time()
  while(abs(rho-rhot)/rhot >.0002 | iter==0) # convergence condition
  {
    # ASSIGN THE COPULA 
    Cop <- cop(dim=2)
    # find the value of theta inducing the assigned Spearman's rho
    # special case to be handled with care: Clayton copula with Spearman's rho >.995
    if(copula=="clayton" & rhos>0.995) {
    theta <- find.theta.clayton(rhos)    # if rhos close to 1, it may not work!
    } else  theta <- iRho(Cop,rhos) 
    # set the copula parameter theta
    Cop   <- cop(param=theta, dim=2)
    # ASSIGN THE MARGINS
    model <- mvdc(Cop,  c(margins,margins), list(par,par))
    # CONSTRUCT the function for the MIXED MOMENT of the MODEL
    integrand0 <- function(arg)
    {
      x1 <- arg[1]
      x2 <- arg[2]
      x1*x2*dMvdc(c(x1,x2),model)
    }
    # COMPUTE THE MIXED MOMENT for the MODEL
    if(margins=="exp" | margins=="weibull"){
    res0 <- cuhre(f = integrand0,
                  lowerLimit = c(0, 0),
                  upperLimit = c(Inf, Inf),
                  relTol = 1e-4, absTol= 1e-12)
    } else if (margins=="norm"){
    res0 <- cuhre(f = integrand0,
                    lowerLimit = c(-Inf, -Inf),
                    upperLimit = c(Inf, Inf),
                    relTol = 1e-4, absTol= 1e-12) 
    }
    # COMPUTE PEARSON'S LINEAR CORRELATION (EVALUATE THE MARGINAL MEANS AND VARIANCES!)
    if(margins=="exp") {
      rho  <- (res0$integral-1)           # for EXP (lambda=1)
    } else if (margins=="weibull") {      # for WEIBULL (scale=2, shape=1.5)  
      rho  <- (res0$integral-Ew()^2)/Vw()
    } else if (margins=="norm") {         # standard normal
      rho <- res0$integral
    }
    iter <- iter+1
    # ADJUST SPEARMAN'S RHO
    rhos <- rhot^(log(rhos)/log(rho)) 
    # RANDOM ADJUSTMENT AFTER THE 7th ITERATION
    a <- ifelse(iter>7, runif(1,-min(1e-5,1-rhos),min(1e-5,1-rhos)),0)
    rhos <- rhos + a
    }
  return(theta)
}
# examples
# find.theta(.9, "gumbel","exp",1)
# find.theta(.2, "clayton","exp",1)
# find.theta(.2, "frank","exp",1)
# find.theta(.2, "gauss","exp",1)

##########################################################################################
##################################### UNIFORM MARGINS ####################################
##########################################################################################
set.seed(12345)
# GUMBEL/CLAYTON/FRANK/GAUSSIAN COPULA
# we consider values of linear correlation rho from 0.01 to 0.99 with constant step=0.01
rho.v   <- (1:99)/100                 # biserial correlation, between 0.01 and 0.99
ratio <- numeric(length(rho.v))
rhopp <- numeric(length(rho.v))
ratio.a <- numeric(length(rho.v))     # asymmetric
rhopp.a <- numeric(length(rho.v))     # asymmetric
ratio.tri.u <- numeric(length(rho.v)) # triserial uniform
rhopp.tri.u <- numeric(length(rho.v)) # triserial uniform
ratio.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
rhopp.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
ratio.tri.a <- numeric(length(rho.v)) # triserial asymmetric
rhopp.tri.a <- numeric(length(rho.v)) # triserial asymmetric
theta <- numeric(length(rho.v))

copula <- dependence

for(i in 1:length(rho.v))
{
  # find the copula parameter inducing rho.v[i] for the copula
  cop <- f.cop(copula)
  theta[i] <- find.theta.cop(copula, rho.v[i])
  cop@parameters <- theta[i]
  # this is the integrand function we have to use in order to compute
  # the point-biserial correlation
  integrand <- function(arg)
  {
    u1 <- arg[1]
    u2 <- arg[2]
    u1*dCopula(c(u1,u2), cop)
  }
  # biserial
  # symmetric case: p_1=p_2=1/2  
  res.1 <- cuhre(f = integrand,
                 lowerLimit = c(0, 0),
                 upperLimit = c(1, 1/2),
                 relTol = 1e-4, absTol= 1e-12)
  #
  res.2 <- cuhre(f = integrand,
                 lowerLimit = c(0, 1/2),
                 upperLimit = c(1, 1),
                 relTol = 1e-4, absTol= 1e-12)
  # asymmetric case: p_1=2/3, p_2=1/3
  res.1a <- cuhre(f = integrand,
                  lowerLimit = c(0,0),
                  upperLimit = c(1,2/3),
                  relTol = 1e-4, absTol= 1e-12)
  #
  res.2a <- cuhre(f = integrand,
                  lowerLimit = c(0, 2/3),
                  upperLimit = c(1, 1),
                  relTol = 1e-4, absTol= 1e-12)
  # triserial
  # discrete uniform
  res.1.tri.u <- cuhre(f = integrand,
                 lowerLimit = c(0, 0),
                 upperLimit = c(1, 1/3),
                 relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.u <- cuhre(f = integrand,
                 lowerLimit = c(0, 1/3),
                 upperLimit = c(1, 2/3),
                 relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.u <- cuhre(f = integrand,
                 lowerLimit = c(0, 2/3),
                 upperLimit = c(1, 1),
                 relTol = 1e-4, absTol= 1e-12)
  # symmetrical unimodal
  res.1.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(1, 1/6),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, 1/6),
                       upperLimit = c(1, 5/6),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, 5/6),
                       upperLimit = c(1, 1),
                       relTol = 1e-4, absTol= 1e-12)
  # asymmetrical
  res.1.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(1, 1/2),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, 1/2),
                       upperLimit = c(1, 5/6),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, 5/6),
                       upperLimit = c(1, 1),
                       relTol = 1e-4, absTol= 1e-12)
  # results
  # biserial
  # symm
  rhopp[i] <- (res.1$integral + 2*res.2$integral - 3/4)/sqrt(1/12*1/4)
  ratio[i] <- (rhopp[i]/rho.v[i])
  # asymm
  rhopp.a[i] <-(res.1a$integral + 2*res.2a$integral - 2/3)/sqrt(1/12*2/9)
  ratio.a[i]<-(rhopp.a[i]/rho.v[i])
  # triserial u
  rhopp.tri.u[i] <-(res.1.tri.u$integral + 2*res.2.tri.u$integral + 3*res.3.tri.u$integral - 2*1/2)/sqrt(1/12*2/3)
  ratio.tri.u[i]<-(rhopp.tri.u[i]/rho.v[i])
  # triserial s
  rhopp.tri.s[i] <-(res.1.tri.s$integral + 2*res.2.tri.s$integral + 3*res.3.tri.s$integral - 2*1/2)/sqrt(1/12*1/3)
  ratio.tri.s[i]<-(rhopp.tri.s[i]/rho.v[i])
  # triserial a
  rhopp.tri.a[i] <-(res.1.tri.a$integral + 2*res.2.tri.a$integral + 3*res.3.tri.a$integral - 5/3*1/2)/sqrt(1/12*5/9)
  ratio.tri.a[i]<-(rhopp.tri.a[i]/rho.v[i])
  
  print(i)
}
# saving the results
name <- paste("uniform",copula,"Rdata",sep=".")
save.image(file=name)

######################################################################################################
######################################### EXPONENTIAL MARGINS ########################################
######################################################################################################
# rm(list=ls())
# library(copula)
# library(cubature)
# set.seed(12345)

copula <- dependence

rho.v <- (1:99)/100 # biserial correlation between 0.01 and 0.99
ratio <- numeric(length(rho.v))
rhopp <- numeric(length(rho.v))
ratio.a <- numeric(length(rho.v))     # asymmetric
rhopp.a <- numeric(length(rho.v))     # asymmetric
ratio.tri.u <- numeric(length(rho.v)) # triserial uniform
rhopp.tri.u <- numeric(length(rho.v)) # triserial uniform
ratio.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
rhopp.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
ratio.tri.a <- numeric(length(rho.v)) # triserial asymmetric
rhopp.tri.a <- numeric(length(rho.v)) # triserial asymmetric
theta <- numeric(length(rho.v))

for(i in 1:length(rho.v))
{
  # find the copula parameter inducing rho.v[i] for the bivariate rv 
  theta[i] <- find.theta(rho.v[i], copula=copula, margins="exp", par=1)
  cop   <- f.cop(copula)
  cop@parameters <- theta[i]
  model <- mvdc(cop,  c("exp","exp"), list(rate=1,rate=1))
  # this is the integrand function we have to use in order to compute
  # the point-biserial correlation
  integrand <- function(arg)
  {
    x1 <- arg[1]
    x2 <- arg[2]
    x1*dMvdc(c(x1,x2), model)
  }
  res.1 <- cuhre(f = integrand,
                 lowerLimit = c(0, 0),
                 upperLimit = c(Inf, qexp(1/2,1)),
                 relTol = 1e-4, absTol= 1e-12)
  
  res.2 <- cuhre(f = integrand,
                 lowerLimit = c(0, qexp(1/2,1)),
                 upperLimit = c(Inf, Inf),
                 relTol = 1e-4, absTol= 1e-12)
  rhopp[i] <- (res.1$integral + 2*res.2$integral - 3/2)/sqrt(1*1/4)
  ratio[i] <- rhopp[i]/rho.v[i]
  # asymmetric
  res.1a <- cuhre(f = integrand,
                 lowerLimit = c(0, 0),
                 upperLimit = c(Inf, qexp(2/3,1)),
                 relTol = 1e-4, absTol= 1e-12)
  
  res.2a <- cuhre(f = integrand,
                 lowerLimit = c(0, qexp(2/3,1)),
                 upperLimit = c(Inf, Inf),
                 relTol = 1e-4, absTol= 1e-12)
  # triserial
  # discrete uniform
  res.1.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(Inf, qexp(1/3)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(0, qexp(1/3)),
                       upperLimit = c(Inf, qexp(2/3)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(0, qexp(2/3)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  # symmetrical unimodal
  res.1.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(Inf, qexp(1/6)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, qexp(1/6)),
                       upperLimit = c(Inf, qexp(5/6)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, qexp(5/6)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  # asymmetrical
  res.1.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(Inf, qexp(1/2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, qexp(1/2)),
                       upperLimit = c(Inf, qexp(5/6)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, qexp(5/6)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  
  # asymm
  rhopp.a[i] <- (res.1a$integral + 2*res.2a$integral - 4/3)/sqrt(1*2/9)
  ratio.a[i] <- rhopp.a[i]/rho.v[i]
  # triserial u
  rhopp.tri.u[i] <-(res.1.tri.u$integral + 2*res.2.tri.u$integral + 3*res.3.tri.u$integral - 2*1)/sqrt(1*2/3)
  ratio.tri.u[i]<-(rhopp.tri.u[i]/rho.v[i])
  # triserial s
  rhopp.tri.s[i] <-(res.1.tri.s$integral + 2*res.2.tri.s$integral + 3*res.3.tri.s$integral - 2*1)/sqrt(1*1/3)
  ratio.tri.s[i]<-(rhopp.tri.s[i]/rho.v[i])
  # triserial a
  rhopp.tri.a[i] <-(res.1.tri.a$integral + 2*res.2.tri.a$integral + 3*res.3.tri.a$integral - 5/3*1)/sqrt(1*5/9)
  ratio.tri.a[i]<-(rhopp.tri.a[i]/rho.v[i])
  print(i)
}
# saving the results
name <- paste("exponential",copula,"Rdata",sep=".")
save.image(file=name)

########################################################################################
#################################### WEIBULL MARGINS ###################################
########################################################################################

# rm(list=ls())
# library(copula)
# library(cubature)
# set.seed(12345)

# GUMBEL/CLAYTON/FRANK/GAUSS COPULA
rho.v <- (1:99)/100 # biserial correlation between 0.01 and 0.99
ratio       <- numeric(length(rho.v))
rhopp       <- numeric(length(rho.v))
ratio.a     <- numeric(length(rho.v)) # asymmetric
rhopp.a     <- numeric(length(rho.v)) # asymmetric
ratio.tri.u <- numeric(length(rho.v)) # triserial uniform
rhopp.tri.u <- numeric(length(rho.v)) # triserial uniform
ratio.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
rhopp.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
ratio.tri.a <- numeric(length(rho.v)) # triserial asymmetric
rhopp.tri.a <- numeric(length(rho.v)) # triserial asymmetric
theta <- numeric(length(rho.v))

copula <- dependence

for(i in 1:length(rho.v))
{
  # find the copula parameter inducing rho.v[i] for the bivariate rv 
  theta[i] <- find.theta(rho.v[i], copula=copula, margins="weibull", par=list(scale=2,shape=1.5))
  cop <- f.cop(copula)
  cop@parameters <- theta[i]
  model <- mvdc(cop,  c("weibull","weibull"), list(list(scale=2,shape=1.5), list(scale=2,shape=1.5)))
  # this is the integrand function we have to use in order to compute
  # the point-biserial correlation
  integrand <- function(arg)
  {
    x1 <- arg[1]
    x2 <- arg[2]
    x1*dMvdc(c(x1,x2), model)
  }
  res.1 <- cuhre(f = integrand,
                 lowerLimit = c(0,0),
                 upperLimit = c(Inf, qweibull(1/2,shape=1.5,scale=2)),
                 relTol = 1e-4, absTol= 1e-12)
  
  res.2 <- cuhre(f = integrand,
                 lowerLimit = c(0, qweibull(1/2,shape=1.5,scale=2)),
                 upperLimit = c(Inf, Inf),
                 relTol = 1e-4, absTol= 1e-12)
  # asymmetric case: p_1=2/3, p_2=1/3
  res.1a <- cuhre(f = integrand,
                 lowerLimit = c(0,0),
                 upperLimit = c(Inf, qweibull(2/3,shape=1.5,scale=2)),
                 relTol = 1e-4, absTol= 1e-12)
  
  res.2a <- cuhre(f = integrand,
                 lowerLimit = c(0, qweibull(2/3,shape=1.5,scale=2)),
                 upperLimit = c(Inf, Inf),
                 relTol = 1e-4, absTol= 1e-12)
  
    
  rhopp[i] <-(res.1$integral + 2*res.2$integral-3/2*Ew())/sqrt(Vw()*1/4)
  ratio[i]<-(rhopp[i]/rho.v[i])
  
  rhopp.a[i] <- (res.1a$integral + 2*res.2a$integral-4/3*Ew())/sqrt(Vw()*2/9)
  ratio.a[i] <- (rhopp.a[i]/rho.v[i])

  # polyserial
  # discrete uniform
  res.1.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(Inf, qweibull(1/3,shape=1.5,scale=2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(0, qweibull(1/3,shape=1.5,scale=2)),
                       upperLimit = c(Inf, qweibull(2/3,shape=1.5,scale=2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(0, qweibull(2/3,shape=1.5,scale=2)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  # symmetrical uniform
  res.1.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(Inf, qweibull(1/6,shape=1.5,scale=2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, qweibull(1/6,shape=1.5,scale=2)),
                       upperLimit = c(Inf, qweibull(5/6,shape=1.5,scale=2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(0, qweibull(5/6,shape=1.5,scale=2)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  # symmetrical uniform
  res.1.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, 0),
                       upperLimit = c(Inf, qweibull(1/2,shape=1.5,scale=2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, qweibull(1/2,shape=1.5,scale=2)),
                       upperLimit = c(Inf, qweibull(5/6,shape=1.5,scale=2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(0, qweibull(5/6,shape=1.5,scale=2)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  
  # triserial u
  rhopp.tri.u[i] <-(res.1.tri.u$integral + 2*res.2.tri.u$integral + 3*res.3.tri.u$integral - 2*Ew())/sqrt(Vw()*2/3)
  ratio.tri.u[i]<-(rhopp.tri.u[i]/rho.v[i])
  # triserial s
  rhopp.tri.s[i] <-(res.1.tri.s$integral + 2*res.2.tri.s$integral + 3*res.3.tri.s$integral - 2*Ew())/sqrt(Vw()*1/3)
  ratio.tri.s[i]<-(rhopp.tri.s[i]/rho.v[i])
  # triserial a
  rhopp.tri.a[i] <-(res.1.tri.a$integral + 2*res.2.tri.a$integral + 3*res.3.tri.a$integral - 5/3*Ew())/sqrt(Vw()*5/9)
  ratio.tri.a[i]<-(rhopp.tri.a[i]/rho.v[i])
  print(i)
  }
# saving the results
name <- paste("Weibull", copula, "Rdata", sep=".")
save.image(file=name)

########################################################################################
#################################### NORMAL MARGINS ###################################
########################################################################################
# rm(list=ls())
# library(copula)
# library(cubature)
# set.seed(12345)

# GUMBEL/CLAYTON/FRANK/GAUSS COPULA
rho.v <- (1:99)/100 # biserial correlation between 0.01 and 0.99
ratio       <- numeric(length(rho.v))
rhopp       <- numeric(length(rho.v))
ratio.a     <- numeric(length(rho.v)) # asymmetric
rhopp.a     <- numeric(length(rho.v)) # asymmetric
ratio.tri.u <- numeric(length(rho.v)) # triserial uniform
rhopp.tri.u <- numeric(length(rho.v)) # triserial uniform
ratio.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
rhopp.tri.s <- numeric(length(rho.v)) # triserial symmetric unimodal
ratio.tri.a <- numeric(length(rho.v)) # triserial asymmetric
rhopp.tri.a <- numeric(length(rho.v)) # triserial asymmetric
theta <- numeric(length(rho.v))

copula <- dependence

for(i in 1:length(rho.v))
{
  # find the copula parameter inducing rho.v[i] for the bivariate rv 
  theta[i] <- find.theta(rho.v[i], copula=copula, margins="norm", par=0)
  cop <- f.cop(copula)
  cop@parameters <- theta[i]
  model <- mvdc(cop,  c("norm","norm"), list(list(mean=0,sd=1), list(mean=0,sd=1)))
  # this is the integrand function we have to use in order to compute
  # the point-biserial correlation
  integrand <- function(arg)
  {
    x1 <- arg[1]
    x2 <- arg[2]
    x1*dMvdc(c(x1,x2), model)
  }
  res.1 <- cuhre(f = integrand,
                 lowerLimit = c(-Inf,-Inf),
                 upperLimit = c(Inf, qnorm(1/2)),
                 relTol = 1e-4, absTol= 1e-12)
  
  res.2 <- cuhre(f = integrand,
                 lowerLimit = c(-Inf, qnorm(1/2)),
                 upperLimit = c(Inf, Inf),
                 relTol = 1e-4, absTol= 1e-12)
  # asymmetric case: p_1=2/3, p_2=1/3
  res.1a <- cuhre(f = integrand,
                  lowerLimit = c(-Inf,-Inf),
                  upperLimit = c(Inf, qnorm(2/3)),
                  relTol = 1e-4, absTol= 1e-12)
  
  res.2a <- cuhre(f = integrand,
                  lowerLimit = c(-Inf, qnorm(2/3)),
                  upperLimit = c(Inf, Inf),
                  relTol = 1e-4, absTol= 1e-12)
  
  
  rhopp[i] <-(res.1$integral + 2*res.2$integral)/sqrt(1*1/4)
  ratio[i]<-(rhopp[i]/rho.v[i])
  
  rhopp.a[i] <- (res.1a$integral + 2*res.2a$integral)/sqrt(1*2/9)
  ratio.a[i] <- (rhopp.a[i]/rho.v[i])
  
  # triserial
  # discrete uniform
  res.1.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, -Inf),
                       upperLimit = c(Inf, qnorm(1/3)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, qnorm(1/3)),
                       upperLimit = c(Inf, qnorm(2/3)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.u <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, qnorm(2/3)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  # symmetrical unimodal
  res.1.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, -Inf),
                       upperLimit = c(Inf, qnorm(1/6)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, qnorm(1/6)),
                       upperLimit = c(Inf, qnorm(5/6)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.s <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, qnorm(5/6)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  # symmetrical uniform
  res.1.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, -Inf),
                       upperLimit = c(Inf, qnorm(1/2)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.2.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, qnorm(1/2)),
                       upperLimit = c(Inf, qnorm(5/6)),
                       relTol = 1e-4, absTol= 1e-12)
  #
  res.3.tri.a <- cuhre(f = integrand,
                       lowerLimit = c(-Inf, qnorm(5/6)),
                       upperLimit = c(Inf, Inf),
                       relTol = 1e-4, absTol= 1e-12)
  
  # triserial u
  rhopp.tri.u[i] <-(res.1.tri.u$integral + 2*res.2.tri.u$integral + 3*res.3.tri.u$integral)/sqrt(1*2/3)
  ratio.tri.u[i]<-(rhopp.tri.u[i]/rho.v[i])
  # triserial s
  rhopp.tri.s[i] <-(res.1.tri.s$integral + 2*res.2.tri.s$integral + 3*res.3.tri.s$integral)/sqrt(1*1/3)
  ratio.tri.s[i]<-(rhopp.tri.s[i]/rho.v[i])
  # triserial a
  rhopp.tri.a[i] <-(res.1.tri.a$integral + 2*res.2.tri.a$integral + 3*res.3.tri.a$integral)/sqrt(1*5/9)
  ratio.tri.a[i]<-(rhopp.tri.a[i]/rho.v[i])
  print(i)
}
# saving the results
name <- paste("normal", copula, "Rdata", sep=".")
save.image(file=name)
