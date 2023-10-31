# install.packages("MAMS")
library("MAMS")

design.0<-mams(K= , # number of experimental arms
               J= , # nuber of stages
               alpha= , # FWER
               power= ,  # Desired Power
               r=1, r0=1,  # Allocation ratio (1 to 1)
               delta= ,    # interesting value of effect 
               delta0= , # uninteresting value of effect
               p=NULL, p0=NULL, # interesting/uninterseting effects on probability scale (leave as is)
               sd=4.4)    # standard deviation
design.0



design.0.sims.null<-mams.sim(nsim=20000, # number of simulations
nMat=matrix(c(  ), nrow=1, ncol=5), # the matrix with sample sizes to be used at each stage
pv=NULL,
l=c(  ),   # lower critical values at each stage
u=c(  ),   # upper critical value at each stage
deltav=c(0,0,0,0), # the treatment effect at investigation arm under the null hypothesis
sd=  )      # standard deviation




