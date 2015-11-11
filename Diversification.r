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

tt.func <- make.bd(tt)
fitDiversitree(tt.func)
tt.func


##Simulate 1000 trees
# General template for simulations:
# In this example, we will just
#  generate a single random number 
#  for each replicate simulation, then
#  store the results in a vector.

REPS <- 5
myResults <- numeric(REPS)
for (i in 1:REPS){
  
  #  Here we would do the simulation: 
  tmp <- rnorm(1)
  
  #   tmp is now the result of our simulation.
  #   We will now store it in the results vector: 
  
  myResults[i] <- tmp
  
}

