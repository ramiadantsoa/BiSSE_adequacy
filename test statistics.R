library(diversitree)
set.seed(1)

### Generate a "real" tree

pars.real <- c(1,0.2)
max.t.real <- 6

phy.real <- tree.bd(pars.real, max.t=max.t.real)
plot(phy.real)
Ntip(phy.real) # 258

### Get its parameter estimates

lnL <- make.bd(phy.real)
ans <- find.mle(lnL, pars.real)

pars <- coef(ans)
max.t <- max(branching.times(phy.real))

### Generate simulated trees

phy.sim <- trees(pars, type="bd", n=100, max.t=max.t)
ntip.sim <- sapply(phy.sim, Ntip)
range(ntip.sim) # 3 323
hist(ntip.sim, breaks=20)

# Note that make.bd() conditions on non-extinction of the tree, as does
# trees().

### Computing phylogenetic diversity

# Faith's PD (total sum of branch length)
PD <- sum(phy.real$edge.length)
all.PD <- sapply(phy.sim, function(x) sum(x$edge.length))
hist(all.PD, breaks=20)


### Computing gamma-statistic (this should be a standard normal distribution)
library(phytools)

ltt(phy.real)$gamma
sim.gamma <- sapply(phy.sim, function(x) ltt(x, plot= FALSE)$gamma)
hist(sim.gamma,breaks = 20)

### Tree imbalance
library(apTreeshape)
t.phy.real <- as.treeshape(phy.real)
colless(t.phy.real)

#t.phy.sim <- sapply(phy.sim, function(x) as.treeshape(x))
# sim.colless <- sapply(t.phy.sim, function(x) colless(x)) this does not work

# prune tree to have the same number of tips
keep <- sample(phy.real$tip.label, 10, replace = FALSE)
pruned.tree <- drop.tip(phy.real, setdiff(phy.real$tip.label, keep))
plot(pruned.tree)

#play with colless for the same tree
set.seed(1)
real.colless <- rep(0, 100) # stores colless for each pruned tree 
for(i in 1:100){
	keep <- sample(phy.real$tip.label, 100, replace = FALSE)
	pruned.tree <- drop.tip(phy.real, setdiff(phy.real$tip.label, keep))
	temp <- as.treeshape(pruned.tree) 
	real.colless[i] <- colless(temp)
}
hist(real.colless,breaks = 20)

#this is for the simulated data
sim.colless <- rep(0, 100) # stores colless for each pruned tree 
sim.colless2 <- numeric()
NTIP <- 100
for(i in 1:100){
  temp.tree <- phy.sim[[i]]
	#keep <- sample(temp.tree$tip.label, min(NTIP, Ntip(temp.tree)), replace = FALSE)
	if(Ntip(temp.tree) >= NTIP){
	temp20 <- numeric()
		for(j in 1:20){
			keep <- sample(temp.tree$tip.label, NTIP, replace = FALSE)
			pruned.tree <- drop.tip(temp.tree, setdiff(temp.tree$tip.label, keep))
			temp <- as.treeshape(pruned.tree) 
			#sim.colless[i] <- colless(temp)
			temp20 <- c(temp20, colless(temp))			
		}
		mtemp <- mean(temp20)
		sim.colless2 <- c(sim.colless2, mtemp)
	}
}
#hist(sim.colless,breaks=20)
dev.off()
par(mfcol = c(1,2))
hist(real.colless,breaks = 20)
hist(sim.colless2,breaks=20)