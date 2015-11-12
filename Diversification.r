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




REPS <- 5
taxon_count <- numeric(REPS)
for (i in 1:REPS){
  pars <- c(10, 5)
  tree <- simulateTree(pars, max.taxa = 100)
  taxon_count[i] <- length(tree$tip.label)
  
}

#    taxon_count is now a vector 
#     of richness: you can plot a histogram
#     with function hist(...)
#     compute mean, quantile, etc.

hist(taxon_count)
