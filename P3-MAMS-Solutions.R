library("MAMS")

design.0<-mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=NULL, p0=NULL,
               delta=2, delta0=0.5, sd=4.4)
summary(design.0)

design.0.sims.null<-mams.sim(nsim=20000, nMat=matrix(c(117), nrow=1, ncol=5),pv=NULL,l=c(2.16),u=c(2.16), deltav=c(0,0,0,0), sd=4.4, ptest=1)
summary(design.0.sims.null)
design.0.sims.null$typeI

design.0.sims.lfc<-mams.sim(nsim=20000, nMat=matrix(c(117), nrow=1, ncol=5),pv=NULL,l=c(2.16),u=c(2.16), deltav=c(2,0.5,0.5,0.5), sd=4.4, ptest=1)
design.0.sims.lfc
design.0.sims.lfc$power


# Design 1
design.1<-mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=NULL, p0=NULL,
               delta=2, delta0=0.5, sd=4.4, lshape="fixed", lfix=0,ushape="fixed", ufix=50)
design.1


design.1.sims.null<-mams.sim(nsim=20000, nMat=matrix(c(59, 118), nrow=2, ncol=5),u=c(Inf, 2.157),
                             l=c(0.000, 2.157),pv=NULL, deltav=c(0,0,0,0), sd=4.4, ptest=NULL)
design.1.sims.null
design.1.sims.null$typeI

design.1.sims.lfc<-mams.sim(nsim=20000, nMat=matrix(c(59, 118), nrow=2, ncol=5),u=c(Inf, 2.157),
                            l=c(0.000, 2.157),pv=NULL, deltav=c(2,0.5,0.5,0.5), sd=4.4, ptest=NULL)
design.1.sims.lfc
design.1.sims.lfc$power


# Design 2


design.2<-mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=NULL, p0=NULL,
               delta=2, delta0=0.5, sd=4.4,lshape="fixed", lfix=0,ushape=function(x)(return(c(4*x/3,x))))
design.2


design.2.sims.null<-mams.sim(nsim=20000, nMat=matrix(c(61, 122), nrow=2, ncol=5),u=c(2.908, 2.181),
                             l=c(0, 2.181),pv=NULL, deltav=c(0,0,0,0), sd=4.4, ptest=NULL)
design.2.sims.null
design.2.sims.null$typeI


design.2.sims.lfc<-mams.sim(nsim=20000, nMat=matrix(c(61, 122), nrow=2, ncol=5),u=c(2.908, 2.181),
                            l=c(0, 2.181),pv=NULL, deltav=c(2,0.5,0.5,0.5), sd=4.4, ptest=NULL)
design.2.sims.lfc
design.2.sims.lfc$power
