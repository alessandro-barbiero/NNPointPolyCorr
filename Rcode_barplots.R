# Discretized distributions examined in the Section
# Point-polyserial correlation under non-normality
op <- par()
par(mfrow=c(3,3), mai=c(0.33,0.3,0.33,0.12), mgp=c(2,1,0),cex=.75)
ylim=c(0,1)
#sym2
f<-c(1/2,1/2)
barplot(f,names.arg = 1:2, ylim=ylim,main="2-sym")
box()
#asym2
f<-c(2/3,1/3)
plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
barplot(f,names.arg = 1:2, ylim=ylim,main="2-asym")
box()
#unif3
f<-c(1/3,1/3,1/3)
barplot(f,names.arg = 1:3, ylim=ylim,main="3-unif")
box()
#symm3
f<-c(1/6,2/3,1/6)
barplot(f,names.arg = 1:3, ylim=ylim,main="3-sym")
box()
#asym3
f<-c(1/2,1/3,1/6)
barplot(f,names.arg = 1:3, ylim=ylim,main="3-asym")
box()
#unif4
f<-rep(1/4,4)
barplot(f,names.arg = 1:4, ylim=ylim,main="4-unif")
box()
#symm3
f<-c(1/8,3/8,3/8,1/8)
barplot(f,names.arg = 1:4, ylim=ylim,main="4-sym")
box()
#asym3
f<-c(1/2,1/4,1/6,1/12)
barplot(f,names.arg = 1:4, ylim=ylim,main="4-asym")
box()
par(op)
# --> Figure 5