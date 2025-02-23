# Two examples where the discretization of one continuous component of a bivariate random vector
# leads to an inflation of Pearson's linear correlation
library(copula)
#
# Frank copula combined with Exponential margins
cop <- frankCopula(dim=2, param=1.25)
mod <- mvdc(copula=cop, margins = c("exp","exp"), paramMargins = list(1,1))
set.seed(01234)
x <- rMvdc(n=1000, mod)
plot(x)
rhox<-cor(x)[1,2]
rhox
y <- x
y[x[,2]<=qexp(1/2,1),2] <- 1
y[x[,2]>qexp(1/2,1) ,2] <- 2
rhoxd<-cor(y)[1,2]
rhoxd
rhoxd/rhox

op<-par()
par(mfrow=c(1,2), mai=c(0.8,0.8,0.4,0.1), mgp=c(2.5,1,0))
plot(x,xlab=expression(x),ylab=expression(y),main=bquote(rho[XY] == .(round(rhox, 3))))
points(x[y[,2]==1,1],x[y[,2]==1,2], col="red")
points(x[y[,2]==2,1],x[y[,2]==2,2], col="cyan")
plot(y,xlab=expression(x),ylab=expression(y[d]),main=bquote(rho[XY[d]] == .(round(rhoxd, 3))),axes=FALSE)
points(y[y[,2]==1,1],y[y[,2]==1,2], col="red")
points(y[y[,2]==2,1],y[y[,2]==2,2], col="cyan")
axis(1)
axis(2, at=1:2, labels=c("1","2"))
box()
par(op)

round(rhox,3)
round(rhoxd,3)

# Gumbel copula combined with Weibull margins
cop <- claytonCopula(dim=2, param=0.15)
model <- mvdc(cop,  c("weibull","weibull"), list(list(scale=2,shape=1.5), list(scale=2,shape=1.5)))
set.seed(01234)
x <- rMvdc(n=1000, model)
rhox <- cor(x)[1,2]
y <- x
y[y[,2]<qweibull(1/3,scale=2,shape=1.5),2] <- 1
y[y[,2]>qweibull(1/3,scale=2,shape=1.5) & y[,2]<qweibull(2/3,scale=2,shape=1.5),2] <- 2
y[y[,2]>qweibull(2/3,scale=2,shape=1.5),2] <- 3
cor(y)
cor(y)/cor(x)
rhoxd<-cor(y)[1,2]
op<-par()
par(mfrow=c(1,2), mai=c(0.8,0.8,0.4,0.1), mgp=c(2.5,1,0))
plot(x,xlab=expression(x),ylab=expression(y),main=bquote(rho[XY] == .(round(rhox, 3))))
points(x[y[,2]==1,1],x[y[,2]==1,2], col="orange")
points(x[y[,2]==2,1],x[y[,2]==2,2], col="green")
points(x[y[,2]==3,1],x[y[,2]==3,2], col="blue")
plot(y,xlab=expression(x),ylab=expression(y[d]),main=bquote(rho[XY[d]] == .(round(rhoxd, 3))),axes=FALSE)
points(y[y[,2]==1,1],y[y[,2]==1,2], col="orange")
points(y[y[,2]==2,1],y[y[,2]==2,2], col="green")
points(y[y[,2]==3,1],y[y[,2]==3,2], col="blue")
axis(1)
axis(2, at=1:3, labels=c("1","2","3"))
box()
par(op)

round(rhox,3)
round(rhoxd,3)
