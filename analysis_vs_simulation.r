library(diversitree)

##### THIS IS TO COMPARE SIMULATION AND ANALYTIC RESULTS FOR PURE BIRTH-DEATH PROCESS #####


#### Section 1 ####
lambda <- 0.1
mu <- 0.05

# simulate trees
NSIM <- 10000
NTIP <- rep(0,NSIM)
t <- 100

for( i in 1:NSIM){
  phy <- tree.bd(c(lambda,mu), max.t = t)
  nnull <- is.null(phy)
  if(nnull == FALSE) NTIP[i] <- Ntip(phy)
  else NTIP[i] <- 0
}

Msimulation <- mean(NTIP)
Vsimulation <- var(NTIP)

#this is analytical result that includes extinction
Manalytic <- exp((lambda - mu)*t )
Vanalytic <-  exp(2*(lambda - mu)*t)*(1 - exp((mu- lambda)*t))*(lambda + mu)/(lambda - mu)

#### Section 2 ####

# MAXTIPS <- 500
TIME <- 50
NSIM <- 200
NTIP <- rep(0, NSIM)

# the default parameter
ts <- 0.2
tc <- 0.2
tm <- 0.03

s0 <- ts*5/6
s1 <- ts*1/6
m0 <- tm/2
m1 <- tm/2
c0 <- tc/2
c1 <- tc/2

# for asymmetric speciation rate
param <- c(s0, s1 , m0, m1, c0, c1)

xhat <- 0.3486121811340026 # this is only true for the 5/6 and 1/6 coeff before, use mathematica to get such number.

for( i in 1:NSIM){
  phy <- tree.bisse(param, max.t = t, x0 = 0)
  nnull <- is.null(phy)
  if(nnull == FALSE) NTIP[i] <- Ntip(phy)
  else NTIP[i] <- 0
}

Msimulation <- mean(NTIP)

preN <- exp((xhat*(s0 - m0) + (1 - xhat)*(s1 - m1))*t)
g <- s0 - m0 - s1 + m1
f1 <- g*(1- xhat)^2/(g*(1 - xhat)^2 - c0)

Manalytic <- preN*(1 - f1) + exp((s0 - m0 - c0/(1- xhat))*t)*f1
