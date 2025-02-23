# plot of the ratio between point-biserial and biserial correlation
# for a bivariate normal rv
# as a function of the threshold omega
biserial<-function(omega)
{
dnorm(omega)/sqrt(pnorm(omega)*(1-pnorm(omega)))
}
op<-par()
op<-par(mai=c(0.9,0.9,0.1,0.1),mgp=c(2.5,1,0),cex=1.25)
plot(function(x) biserial(x),xlim=c(-4,4),ylab=expression(rho[X[1]*Z[2]]/rho[Z[1]*Z[2]]),xlab=expression(omega))
abline(h=0)
par(op)
