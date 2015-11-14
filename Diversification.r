## Diversification exercise
rm(list = ls())
library(geiger)
library(phytools)
library(diversitree)
source('rabosky_functions.R')

# Using function simulateTree to simulate a birth death tree:
# Example, lambda = 10, mu = 5
pars <- c(10, 0); #in order: lambda, mu
tt <- simulateTree(pars, max.taxa=100)

##simulate under BD and fit gamma
## is confidence in speciation or extinction better?
## regression
## show loop construct
## pull lambda and mu

plot(tt,  show.tip.label = F)

# use function diversitree::make.bd  
# fit model with fitDiversitree
# Example:
#   my_likelihood_fxn <- make.bd(....)

tt.func <- make.bd(tt) #this creates a likelihood function based on the tree tt
fitDiversitree(tt.func) 
tt.func


##Simulate 100 trees

REPS <- 100
lambdaVal <- numeric(REPS)
muVal <- numeric(REPS)
for (i in 1:REPS){

  #  Here we would do the simulation: 
  pars <- c(10, 5); #in order: lambda, mu
  tt <- simulateTree(pars, max.taxa=100)
  tt.func <- make.bd(tt)
  tmp <- fitDiversitree(tt.func)
  
  lambdaVal[i] <- tmp $pars["lambda"]
  muVal[i] <- tmp $pars["mu"]
}
par(mfrow = c(2,1))
hist(lambdaVal)
hist(muVal)


#################################################################################
#Variation in bd process
REPS <- 100
taxon_count <- numeric(REPS)
for (i in 1:REPS){
  pars <- c(10, 5)
  tree <- simulateTree(pars, max.t = 1)
  taxon_count[i] <- length(tree$tip.label)
  
}

#    taxon_count is now a vector 
#     of richness: you can plot a histogram
#     with function hist(...)
#     compute mean, quantile, etc.

hist(taxon_count, main = "Stochasticity in birth-death process in taxon count")

####################################################################

#Lineage Through Time
## 
pars <- c(10, 5); #in order: lambda, mu
tt <- simulateTree(pars, max.taxa=100)
ttEx <- sim.bdtree(b = 10, d = 1, stop = "taxa", t = 100, n = 100 )
plot(ttEx, show.tip.label = F)
livingOnly <- drop.fossil(ttEx)  #prune out the extinct taxa
plot(livingOnly, show.tip.label = F)
par(mfrow = c(1,2))
plot(ttEx, show.tip.label = F, main = "extant and extinct")
plot(livingOnly, show.tip.label = F, main = "extant only")

#What does extinction do to the shape of the tree? Simulate 5 trees with and without extinction

pars <- c(10, 5)
treeEx <- simulateTree(pars, max.taxa = 100)
treeNoEx <- simulateTree(c(5, 0), max.taxa = 100)


par(mfrow = c(1,2))
plot(treeEx, show.tip.label = F, main = "extinction")
plot(treeNoEx, show.tip.label = F, main = "no extinction")

#Diversity-Dependence Model
library(DDD)   #package for diversity dependent diversification
ddTree <- dd_sim(c(0.12,0.02,100),100) # pars: c(lambda, mu, K)
#note that ddTree is a list. The tree with extinct members pruned is element 2
ddLivingOnly <- ddTree[[1]]
par(mfrow = c(1,2))
plot(ddLivingOnly, show.tip.label = F, main = "diversity dependence")
ltt(ddLivingOnly)


##Diversification Analysis
require(geiger) # loads geiger
require(laser) # loads laser
par(mfrow = c(1, 1))
snake_tree<-read.tree("homalops.phy") 
# plot the tree
plot(snake_tree)
pdf('snaketree.pdf')
plot(ladderize(snake_tree)) 
dev.off()
ltt.plot(snake_tree, log="y")

layout(matrix(1:2, 2, 1))
data(bird.families)
data(bird.orders)
#plot them with titles
ltt.plot(bird.families)
title("Lineages Through Time Plot of the Bird Families")
ltt.plot(bird.families, log = "y")
title(main = "Lineages Through Time Plot of the Bird Families",
      sub = "(with logarithmic transformation of the y-axis)")

layout(matrix(1:4, 2, 2))
plot(bird.families, show.tip.label = FALSE)
ltt.plot(bird.families, main = "Bird families")
plot(snake_tree, show.tip.label = FALSE)
ltt.plot(snake_tree, main = "Homalopsid species")

layout(matrix(1))
mltt.plot(bird.families, bird.orders)

#using laser to plot ltt
library(laser)
snake.times<-getBtimes('homalops.phy')
snake_tree<-read.tree("homalops.phy") 
snake.times2<-getBtimes(string = write.tree(snake_tree))
snake.times2
plotLtt(snake.times2)
#compare
layout(matrix(1:2, 2,1))
ltt.plot(snake_tree, log="y", main = "APE LTT")
plotLtt(snake.times) 

##gamma stat
snake_gamma <- gammaStat(snake_tree)
snake_gamma 
