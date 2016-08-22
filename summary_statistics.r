# clear all variables
rm(list=ls(all=TRUE))

#### generating and saving trees with different parameters ####
setwd("C:/Users/tramiada/OneDrive/Projects/BiSSE_adequacy")


library(diversitree)
library(apTreeshape)
library(phytools)
library(moments)
#  other simulation parameters

initstate <- 0

# generate tree

# MAXTIPS <- 500
TIME <- 50
NSIM <- 200
stats <-matrix(rep(0, NSIM * 10), nrow = NSIM)

# the default parameter
s0 <- 0.1
s1 <- 0.1
m0 <- 0.03
m1 <- 0.03
c0 <- 0.1
c1 <- 0.1
ts <- s1 + s0
tc <- c1 + c0
tm <- m1 + m0

pars <- c(ts/2, ts/2, tm/2, tm/2, tc/2, tc/2)

# for asymmetric speciation rate
pars1 <- c(ts/3, ts*2/3 , m0, m1, c0, c1)
pars2 <- c(ts/6, ts*5/6 , m0, m1, c0, c1)

# for asymmetric colonization rate
pars3 <- c(ts/2, ts/2 , m0, m1, tc*2/3, tc/3)
pars4 <- c(ts/2, ts/2 , m0, m1, tc*5/6, tc/6)

# for asymmetric extinction rate
pars5 <- c(ts/2, ts/2 , tm*2/3, tm*1/3, tc/2, tc/2)
pars6 <- c(ts/2, ts/2 , tm*5/6, tm*1/6, tc/2, tc/2)

# for asymmetric extinction rate
pars7 <- c(ts/2, ts/2 , 3*tm*2/3, 3*tm*1/3, tc/2, tc/2)
pars8 <- c(ts/2, ts/2 , 3*tm*5/6, 3*tm*1/6, tc/2, tc/2)
pars9 <- c(ts/2, ts/2 , 3*tm*1/2, 3*tm*1/2, tc/2, tc/2)


########### generating dataset ################

# assign a seed for replicability
set.seed(1)

# this is for default
totaltime <- system.time(
  for(ii in 1:NSIM){
  #initialize the condition for each simulation
  nnull <- TRUE
  ntip <- 1
    while(nnull == TRUE || ntip < 20){
      phy <- tree.bisse(pars, max.t = TIME, x0 = initstate)
      nnull <- is.null(phy)
      if(nnull == FALSE) ntip <- Ntip(phy)
    }
    stats[ii,1] <- ntip
    stats[ii,2] <- sum(phy$tip.state)/ntip

    #branch length analysis
    el <- phy$edge.length
    stats[ii,3] <- ltt(phy, plot = FALSE)$gamma
    stats[ii,4] <- mean(el)
    stats[ii,5] <- var(el)
    stats[ii,6] <- skewness(el)
    stats[ii,7] <- kurtosis(el)

    #tree shape analysis
    mod.phy <-as.treeshape(phy)
    stats[ii,8] <- colless(mod.phy, norm = "pda")
    stats[ii,9] <- sackin(mod.phy, norm = "pda")

    # age of tips as function of state
    n <-length(phy$tip.label)
    tipage<-setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y==x),y=phy$edge[,2])],phy$tip.label)

    state1 <- which(phy$tip.state == 1)
    state0 <- which(phy$tip.state == 0)
    if(length(state1) == 0){
      stats[ii,10] <- Inf
    }
    else{
      tipage1 <- tipage[state1]
      tipage0 <- tipage[state0]
      stats[ii,10] <- mean(tipage0)/mean(tipage1)
    }
  }
)
totaltime

finalresult <- data.frame(stats)
colnames(finalresult)<- c("ntip","freq1","gamma", "mean","variance", "skweness",
 "kurtosis", "colless","sackin", "relative.tip.age")

filename <- paste("summary_pars_default.csv",sep = "")
write.csv(finalresult,filename)

########### PERTURB SPECIATION ##############
# assign a seed for replicability
set.seed(2)

totaltime <- system.time(
  for(ii in 1:NSIM){
    #initialize the condition for each simulation
    nnull <- TRUE
    ntip <- 1
    while(nnull == TRUE || ntip < 20){
      phy <- tree.bisse(pars1, max.t = TIME, x0 = initstate)
      # phy <- tree.bisse(pars2, max.t = TIME, x0 = initstate)
      head(phy)
      nnull <- is.null(phy)
      if(nnull == FALSE) ntip <- Ntip(phy)
    }
    stats[ii,1] <- ntip
    stats[ii,2] <- sum(phy$tip.state)/ntip

    #branch length analysis
    el <- phy$edge.length
    stats[ii,3] <- ltt(phy, plot = FALSE)$gamma
    stats[ii,4] <- mean(el)
    stats[ii,5] <- var(el)
    stats[ii,6] <- skewness(el)
    stats[ii,7] <- kurtosis(el)

    #tree shape analysis
    mod.phy <-as.treeshape(phy)
    stats[ii,8] <- colless(mod.phy, norm = "pda")
    stats[ii,9] <- sackin(mod.phy, norm = "pda")

    # age of tips as function of state
    n <-length(phy$tip.label)
    tipage<-setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y==x),y=phy$edge[,2])],phy$tip.label)

    state1 <- which(phy$tip.state == 1)
    state0 <- which(phy$tip.state == 0)
    if(length(state1) == 0){
      stats[ii,10] <- Inf
    }
    else{
      tipage1 <- tipage[state1]
      tipage0 <- tipage[state0]
      stats[ii,10] <- mean(tipage0)/mean(tipage1)
    }
  }
)

  totaltime

  finalresult <- data.frame(stats)
  colnames(finalresult)<- c("ntip","freq1","gamma", "mean","variance", "skweness",
  "kurtosis", "colless","sackin", "relative.tip.age")

  filename <- paste("summary_pars_2s1_t_",TIME,".csv",sep = "")
  write.csv(finalresult, filename)

  filename <- paste("summary_pars_5s1_t_",TIME,".csv",sep = "")
  write.csv(finalresult, filename)

########### PERTURB COLONIZATION ##############
# assign a seed for replicability
set.seed(2)

totaltime <- system.time(
  for(ii in 1:NSIM){
  #initialize the condition for each simulation
  nnull <- TRUE
  ntip <- 1
    while(nnull == TRUE || ntip < 20){
      phy <- tree.bisse(pars3, max.t = TIME, x0 = initstate)
      # phy <- tree.bisse(pars4, max.t = TIME, x0 = initstate)
      nnull <- is.null(phy)
      if(nnull == FALSE) ntip <- Ntip(phy)
    }
    stats[ii,1] <- ntip
    stats[ii,2] <- sum(phy$tip.state)/ntip

    #branch length analysis
    el <- phy$edge.length
    stats[ii,3] <- ltt(phy, plot = FALSE)$gamma
    stats[ii,4] <- mean(el)
    stats[ii,5] <- var(el)
    stats[ii,6] <- skewness(el)
    stats[ii,7] <- kurtosis(el)

    #tree shape analysis
    mod.phy <-as.treeshape(phy)
    stats[ii,8] <- colless(mod.phy, norm = "pda")
    stats[ii,9] <- sackin(mod.phy, norm = "pda")

    # age of tips as function of state
    n <-length(phy$tip.label)
    tipage<-setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y==x),y=phy$edge[,2])],phy$tip.label)

    state1 <- which(phy$tip.state == 1)
    state0 <- which(phy$tip.state == 0)
    if(length(state1) == 0){
      stats[ii,10] <- Inf
    }
    else{
      tipage1 <- tipage[state1]
      tipage0 <- tipage[state0]
      stats[ii,10] <- mean(tipage0)/mean(tipage1)
    }
  }
)
totaltime

finalresult <- data.frame(stats)
colnames(finalresult)<- c("ntip","freq1","gamma", "mean","variance", "skweness",
 "kurtosis", "colless","sackin", "relative.tip.age")

filename <- paste("summary_pars_2c0_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)

filename <- paste("summary_pars_5c0_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)


########### PERTURB EXTINCTION I ##############
# assign a seed for replicability
set.seed(2)

totaltime <- system.time(
  for(ii in 1:NSIM){
  #initialize the condition for each simulation
  nnull <- TRUE
  ntip <- 1
    while(nnull == TRUE || ntip < 20){
      phy <- tree.bisse(pars5, max.t = TIME, x0 = initstate)
      # phy <- tree.bisse(pars6, max.t = TIME, x0 = initstate)
      nnull <- is.null(phy)
      if(nnull == FALSE) ntip <- Ntip(phy)
    }
    stats[ii,1] <- ntip
    stats[ii,2] <- sum(phy$tip.state)/ntip

    #branch length analysis
    el <- phy$edge.length
    stats[ii,3] <- ltt(phy, plot = FALSE)$gamma
    stats[ii,4] <- mean(el)
    stats[ii,5] <- var(el)
    stats[ii,6] <- skewness(el)
    stats[ii,7] <- kurtosis(el)

    #tree shape analysis
    mod.phy <-as.treeshape(phy)
    stats[ii,8] <- colless(mod.phy, norm = "pda")
    stats[ii,9] <- sackin(mod.phy, norm = "pda")

    # age of tips as function of state
    n <-length(phy$tip.label)
    tipage<-setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y==x),y=phy$edge[,2])],phy$tip.label)

    state1 <- which(phy$tip.state == 1)
    state0 <- which(phy$tip.state == 0)
    if(length(state1) == 0){
      stats[ii,10] <- Inf
    }
    else{
      tipage1 <- tipage[state1]
      tipage0 <- tipage[state0]
      stats[ii,10] <- mean(tipage0)/mean(tipage1)
    }
  }
)
totaltime

finalresult <- data.frame(stats)
colnames(finalresult)<- c("ntip","freq1","gamma", "mean","variance", "skweness",
 "kurtosis", "colless","sackin", "relative.tip.age")

filename <- paste("summary_pars_2m0_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)

filename <- paste("summary_pars_5m0_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)


########### PERTURB EXTINCTION II ##############
# assign a seed for replicability
set.seed(2)

totaltime <- system.time(
  for(ii in 1:NSIM){
  #initialize the condition for each simulation
  nnull <- TRUE
  ntip <- 1
    while(nnull == TRUE || ntip < 20){
      phy <- tree.bisse(pars7, max.t = TIME, x0 = initstate)
      # phy <- tree.bisse(pars8, max.t = TIME, x0 = initstate)
      # phy <- tree.bisse(pars9, max.t = TIME, x0 = initstate)
      nnull <- is.null(phy)
      if(nnull == FALSE) ntip <- Ntip(phy)
    }
    stats[ii,1] <- ntip
    stats[ii,2] <- sum(phy$tip.state)/ntip

    #branch length analysis
    el <- phy$edge.length
    stats[ii,3] <- ltt(phy, plot = FALSE)$gamma
    stats[ii,4] <- mean(el)
    stats[ii,5] <- var(el)
    stats[ii,6] <- skewness(el)
    stats[ii,7] <- kurtosis(el)

    #tree shape analysis
    mod.phy <-as.treeshape(phy)
    stats[ii,8] <- colless(mod.phy, norm = "pda")
    stats[ii,9] <- sackin(mod.phy, norm = "pda")

    # age of tips as function of state
    n <-length(phy$tip.label)
    tipage<-setNames(phy$edge.length[sapply(1:n,function(x,y)   which(y==x),y=phy$edge[,2])],phy$tip.label)

    state1 <- which(phy$tip.state == 1)
    state0 <- which(phy$tip.state == 0)
    if(length(state1) == 0){
      stats[ii,10] <- Inf
    }
    else{
      tipage1 <- tipage[state1]
      tipage0 <- tipage[state0]
      stats[ii,10] <- mean(tipage0)/mean(tipage1)
    }
  }
)
totaltime

finalresult <- data.frame(stats)
colnames(finalresult)<- c("ntip","freq1","gamma", "mean","variance", "skweness",
 "kurtosis", "colless","sackin", "relative.tip.age")

filename <- paste("summary_pars_2m03_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)

filename <- paste("summary_pars_5m03_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)

filename <- paste("summary_pars_default3_t_",TIME,".csv",sep = "")
write.csv(finalresult, filename)


################# USE THIS FOR PLOTTING IN R
#
# ########### plotting dataset ################
#
# # importing dataset for different speciation
# d1 <- read.csv( paste("summary_pars_default.csv",sep = ""))
# d2 <- read.csv( paste("summary_pars_2s1.csv",sep = ""))
# d3 <- read.csv( paste("summary_pars_5s1.csv",sep = ""))
#
#
# alld <- data.frame(d1$freq1, d2$freq1, d3$freq1)
#
# col <- c("#eaab00", "#004165", "#618e02")
# profiles.plot(alld, col.line=col, las=1)
# abline(v=c(0.5, 2/3, 5/6), col=col, lty=2, lwd = 2)
#
#
# # importing dataset for different colonization
# d1 <- read.csv( paste("summary_pars_default.csv",sep = ""))
# d2 <- read.csv( paste("summary_pars_2c0.csv",sep = ""))
# d3 <- read.csv( paste("summary_pars_5c0.csv",sep = ""))
#
# alld <- data.frame(d1$freq1, d2$freq1, d3$freq1)
#
# col <- c("#eaab00", "#004165", "#618e02")
# profiles.plot(alld, col.line=col, las=1)
# abline(v=c(0.5, 2/3, 5/6), col=col, lty=2, lwd = 2)
#
# ## TO DO look if the results are robust?? ###
