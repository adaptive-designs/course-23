library("PASWR")
library("rBeta2009")

### Equal Randomisation ###
two_arm_ER <- function(N=200, p = c(0,0),  burnin=3, Z_corrected=FALSE){
  # Generate Samples 
  n1 = min(max(rbinom(1, N, 0.5),burnin), N-burnin)
  n0 = N - n1
  x0 <- rbinom(n0,1, p[1]) 
  x1 <- rbinom(n1, 1, p[2])
  
  # Response
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  ##### Z - Test 
  # One-Point-Sample-Correction
  sigx <- sd(x0); sigy <- sd(x1)
  if(sigx==0 && sigy==0){
    if(mean(x0)==mean(x1)){Z = 0}
    if(mean(x0)<mean(x1)){Z = -Inf}
    if(mean(x0)>mean(x1)){Z = Inf}
    Z_P <-2*pnorm(-1*abs(Z))
  } else {
    Z_P <- PASWR::z.test(x0, x1, sigma.x =sigx, sigma.y = sigy )$p.value
  }
  return(c(response, Z_P, n1/N))
}

### PW ###
PW = function(N=200, p = c(0,0),  burnin=3, Z_corrected=FALSE){
  
  allocation = response = rep(0,N)
  if(burnin>0){
    allocation[1:burnin] = 0
    response[1:burnin] = rbinom(burnin, 1, p[1])
    allocation[(burnin+1):(2*burnin)] = 1
    response[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }
  
  for (i in (2*burnin+1):N){
    if(response[i-1]==1){
      allocation[i]=allocation[i-1]
    }
    else{ 
      allocation[i]=abs(allocation[i-1]-1)
    }
    response[i] = rbinom(1,1,p[allocation[i]+1])
  }
  
  # Response
  x0 <- response[allocation==0]; x1 <- response[allocation==1]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  ##### Z - Test 
  # One-Point-Sample-Correction
  sigx <- sd(x0); sigy <- sd(x1)
  if(sigx==0 && sigy==0){
    if(mean(x0)==mean(x1)){Z = 0}
    if(mean(x0)<mean(x1)){Z = -Inf}
    if(mean(x0)>mean(x1)){Z = Inf}
    Z_P <-2*pnorm(-1*abs(Z))
  } else {
    Z_P <- PASWR::z.test(x0, x1, sigma.x =sigx, sigma.y = sigy )$p.value
  }
  return(c(response, Z_P, n1/N))
}

### RPW ###
randomiseRPW = function(N=200, p = c(0,0),  burnin=3, Z_corrected=FALSE){
  
  allocation = response = rep(0,N)
  if(burnin>0){
    allocation[1:burnin] = 0
    response[1:burnin] = rbinom(burnin, 1, p[1])
    allocation[(burnin+1):(2*burnin)] = 1
    response[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }
  #urn = c(0,1, rep(NA,N)) without burnin
  urn = c(rep(0, burnin),rep(1, burnin), rep(NA,(N-2*burnin)))
  
  for (i in (2*burnin+1):N){
    ball = sample(urn[1:(i-1)], 1) # chose one ball
    allocation[i] = ball
    response[i] = rbinom(1,1,p[allocation[i]+1])
    urn[i] = ifelse(response[i], ball, 1-ball)
  }
  
  # Response
  x0 <- response[allocation==0]; x1 <- response[allocation==1]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  ##### Z - Test 
  # One-Point-Sample-Correction
  sigx <- sd(x0); sigy <- sd(x1)
  if(sigx==0 && sigy==0){
    if(mean(x0)==mean(x1)){Z = 0}
    if(mean(x0)<mean(x1)){Z = -Inf}
    if(mean(x0)>mean(x1)){Z = Inf}
    Z_P <-2*pnorm(-1*abs(Z))
  } else {
    Z_P <- PASWR::z.test(x0, x1, sigma.x =sigx, sigma.y = sigy )$p.value
  }
  return(c(response, Z_P, n1/N))
}

### Optimal Proportions ###
two_arm_OP <- function(N=200, p = c(0,0),  burnin=3, Z_corrected=FALSE, proportion="Neyman"){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = rbinom(burnin, 1, p[1])
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }
  
  for (i in (2*burnin+1):N){
    n1 = sum(A==1, na.rm = TRUE)
    n2 = sum(A==2, na.rm = TRUE)
    if(proportion=="minF"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt(p1.hat) + sqrt(p2.hat)) == 0){
          prop = 0.5
        } else {
          prop = sqrt(p2.hat)/(sqrt(p1.hat) + sqrt(p2.hat))
        }
    }
    if(proportion=="maxMR"){
        s1 = sum(X[A==1], na.rm = TRUE)
        p1.hat = (s1)/(n1)
        s2 = sum(X[A==2], na.rm = TRUE)
        p2.hat = (s2)/(n2)
        if((sqrt(1-p1.hat) + sqrt(1-p2.hat)) == 0){
          prop = 0.5
        } else {
          prop = sqrt(1-p1.hat)/(sqrt(1-p1.hat) + sqrt(1-p2.hat))
        }
    }
    if(proportion=="Neyman"){
      sig1_N <- sd(X[A==1], na.rm = TRUE)
      sig2_N <- sd(X[A==2], na.rm = TRUE)
      if((sig1_N+sig2_N)==0){
        prop <- 0.5
      } else {
        prop <- 1 - sig1_N/(sig1_N+sig2_N)
      }
    }
    if(prop==0){prop <- 1/(n1+n2)}
    if(prop==1){prop <- 1-1/(n1+n2)}
    A[i] = rbinom(1, 1, prop) + 1
    X[i]  <- rbinom(1,1,p[A[i]])
  }
  # Response
  x0 <- X[A==1]; x1 <- X[A==2]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  ##### Z - Test 
  # One-Point-Sample-Correction
  sigx <- sd(X[A==1], na.rm = TRUE); sigy <- sd(X[A==2], na.rm = TRUE)

  if(sigx==0 && sigy==0){
    if(mean(x0)==mean(x1)){Z = 0}
    if(mean(x0)<mean(x1)){Z = -Inf}
    if(mean(x0)>mean(x1)){Z = Inf}
    Z_P <-2*pnorm(-1*abs(Z))
  } else {
    Z_P <- PASWR::z.test(x0, x1, sigma.x =sigx, sigma.y = sigy )$p.value
  }
  return(c(response, Z_P, n1/N))
}


### Bayesian RAR ###

# define a function to select the arm with the highest expected value
select_arm <- function(num_arms, success_counts, failure_counts) {
  # generate a beta distribution for each arm
  values <- sapply(1:num_arms, function(i) generate_beta(success_counts[i] + 1, failure_counts[i] + 1))
  # select the arm with the highest expected value
  return(which.max(values))
}

generate_beta <- function(alpha, beta) {
  rbeta(1, alpha, beta)
}

BRAR = function(N=200, p = c(0,0),burnin=3, Z_corrected=FALSE){
  
  A = rep(NA, N)   # Vector of allocations
  X = rep(NA, N)   # Vector of responses
  
  #Burnin
  if(burnin>0){
    A[1:burnin] = 1
    X[1:burnin] = rbinom(burnin, 1, p[1])
    A[(burnin+1):(2*burnin)] = 2
    X[(burnin+1):(2*burnin)] = rbinom(burnin, 1, p[2])
  }
  
  for (i in (2*burnin+1):N){
    #Count Successes and Failures
    success_counts <- sapply(1:2, function(i) sum(X[A==i], na.rm=TRUE))
    failure_counts <- sapply(1:2, function(i) length(X[A==i])) - success_counts
    A[i] = select_arm(num_arms=2, success_counts, failure_counts)
    X[i] = rbinom(1,1,p[A[i]])
  }
  
  # Response
  x0 <- X[A==1]; x1 <- X[A==2]
  n0 <- length(x0); n1 <- length(x1)
  s0 <- sum(x0); s1 <- sum(x1)
  response <- (s0+s1)/N
  
  ##### Z - Test 
  # One-Point-Sample-Correction
  sigx <- sd(X[A==1], na.rm = TRUE); sigy <- sd(X[A==2], na.rm = TRUE)
  
  if(sigx==0 && sigy==0){
    if(mean(x0)==mean(x1)){Z = 0}
    if(mean(x0)<mean(x1)){Z = -Inf}
    if(mean(x0)>mean(x1)){Z = Inf}
    Z_P <-2*pnorm(-1*abs(Z))
  } else {
    Z_P <- PASWR::z.test(x0, x1, sigma.x =sigx, sigma.y = sigy )$p.value
  }
  return(c(response, Z_P, n1/N))
}


### Simulation ###
simulation <- function(N=3000, p = c(0.941,0.991), burnin=3, nsim= 10^4, method="ER", proportion=NULL){
  
  ########## Theoretical Variance Calculation ############################
  fun <- function(x){x*(1-x)}
  variance <-  data.frame(lapply(p,fun))
  
  sim = matrix(nrow = nsim, ncol = 3)
  ########## Select Method ############################
  for(i in 1:nsim){
    if(method=="ER"){
      sim[i,] <- two_arm_ER(N=N,  p =p,  burnin = burnin, Z_corrected=FALSE)
    } 
    else if(method=="PW"){
      sim[i,] <- PW(N=N,  p =p,  burnin = burnin, Z_corrected=FALSE)
    }
    else if(method=="RPW"){
      sim[i,] <- randomiseRPW(N=N,  p =p,  burnin = burnin, Z_corrected=FALSE)
    }
    else if(method=="OP"){
      sim[i,] <- two_arm_OP(N=N,  p =p,  burnin = burnin, Z_corrected=FALSE, proportion=proportion)
    }
    else if(method=="BRAR"){
      sim[i,] <- BRAR(N=N,  p =p,  burnin = burnin, Z_corrected=FALSE)
    }
  }
  
  ### Output
  reject_Z <- sum(sim[,2] < 0.05)
  power <- round(reject_Z/nsim*100,3)
  per.sup <- round(mean(sim[,3]),4)
  var.per.sup <- round(var(sim[,1]),6)
  emr <- round(mean(sim[,1]),4)
  var.emr <- round(var(sim[,1]),6)
  cat("Power:", power, "%", "\n")
  cat("Percentage of Patients assiegned to treatment arm 1:", per.sup, "\n")
  cat("Variance of percentage of Patients assiegned to treatment arm 1:", var.per.sup, "\n")
  cat("Expected Mean Response:", emr, "\n")
  cat("Variance of Expected Mean Response:", var.emr, "\n")
}

### Results ###
# Enter Burn-In per ARM!
simulation(N=366, p = c(0.941,0.991), burnin=2, nsim= 10^4, method="ER", proportion=NULL) 