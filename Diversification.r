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
par(mfrow = c(1, 1)) 
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


treeExlist <- list()  #make empty list for trees to fill in
treeNoExlist <- list()
for (i in 1:5) {
  parsbd <- c(10, 3)
  treeEx <- simulateTree(parsbd, max.taxa = 100, include.extinct = T)
  treeExPruned <- drop.fossil(treeEx)
  treeExlist[[i]] <- treeExPruned  #fill in list with trees
  parspb <- c(7, 0)
  treeNoEx <- simulateTree(parspb, max.taxa = 100)
  treeNoExlist[[i]] <- treeNoEx
}

#par(mfrow = c(5,2), mar=c(0.8, 0.8, 0.8, 0.8)) #beware of margins, (margins too big error)
layout(matrix(1:10, 5, 2))
for (i in 1:5) {
  plot(treeExlist[[i]], show.tip.label = F, main = "extinction")
  plot(treeNoExlist[[i]], show.tip.label = F, main = "no extinction")
}

for (i in 1:5) {
  ltt(treeExlist[[i]])
  ltt(treeNoExlist[[i]])
}


#Diversity-Dependence Model
library(DDD)   #package for diversity dependent diversification
ddTree <- dd_sim(c(0.12,0.02,100),100) # pars: c(lambda, mu, K)
#note that ddTree is a list. The tree with extinct members pruned is element 1
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

#############################################################
##gamma stat
snake_gamma <- gammaStat(snake_tree)
snake_gamma 

##mccrTest(CladeSize, NumberMissing, NumberOfReps, ObservedGamma = NULL, fname=NULL)
mccrResults <- mccrTest(34, 13, 100, ObservedGamma=snake_gamma)  #mccrTest requires laser package
mccrResults
hist(mccrResults$null.gamma)
abline(v = mccrResults$critical.value, col = 'red')


age <- 22
richness <- 34
snakebirth =  log(richness)/age
snakebirth

#this simulates gamma values when trees are undersampled.
#we will grow trees with n=34 and prune them down to 13 taxa

num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim_tree <- sim.bdtree(snakebirth, d=0, stop = "taxa", n=34)
  prune <- drop.random(sim_tree, 13) # prune down to the # of taxa in the phylogeny (34-13=21)
  g1_null[i] <- gammaStat(prune)
}
# create a histogram of the null distribution
hist(g1_null)

#arrow indicates where the observed gamma falls in the null you just generated
arrows(snake_gamma, 40, snake_gamma, 0, col="red", lwd=2) 

# Which of the null values are smaller (more negative) than the data?
smallerNull <- g1_null<=snake_gamma
# How many TRUEs are there?
count<-sum(smallerNull)
count
# finally, what is the p-value?
mccr_pval <- (count+1)/(num_simulations+1)
mccr_pval

#############################################################
##Fitting diversification models to branching times
#laser works with branching times instead of tree
snake.times <- getBtimes('homalops.phy')

pb <- pureBirth(snake.times) #pb
aic.pb <- pb$aic#store aic score
bd <- bd(snake.times) # birth-death
aic.bd <- bd$aic #store aic
ddl <- DDL(snake.times) # logistic density dependent diversification
aic.ddl <- ddl$aic
ddx <- DDX(snake.times) # exponential density dependent diversification
aic.ddx <- ddx$aic
#the MCCR test above cannot discriminate between early rapid speciation and early slow extinction.
#these models allow those to scenarios to be compared 
spvar <- fitSPVAR(snake.times, init=c(2, .2, .1))#variable speciation: exponentially declining speciation, constant extinction rate
aic.spvar <- spvar$aic
exvar <- fitEXVAR(snake.times, init=c(.3, .01, 1))#variable extinction: constant speciation, exponentially increasing extinction rate
aic.exvar <- exvar$aic
both <- fitBOTHVAR(snake.times, init=c(.3, .5, .1, .5))#variable extinction and speciation
aic.both <- both$aic

#make a table and compare them all
aic.all <- cbind(aic.pb, aic.bd, aic.ddl, aic.ddx, aic.spvar, aic.exvar, aic.both)
foo <- function(x) x-x[which(x==min(x))] #quick function to compare all models to the model with lowest aic score
daic <- t(apply(aic.all, 1, foo))#apply that function to everything in aic.all
cat("Table of delta-aic values; zero - best model\n")
print(daic, digits=2)

###################################################################

balistoid.tree <- read.tree("balistoids.phy")
balistoid.times <- getBtimes("balistoids.phy")

balistoid.gamma <- gammaStat(balistoid.tree)
balistoid.gamma 

age <- 48
richness <- 86
balistoidbirth =  log(richness)/age
balistoidbirth
num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim_tree <- sim.bdtree(balistoidbirth, d=0, stop = "taxa", n=100)
  prune <- drop.random(sim_tree, 100-86) # prune down to the # of taxa in the phylogeny (34-13=21)
  g1_null[i] <- gammaStat(prune)
}
hist(g1_null)

arrows(balistoid.gamma, 40, balistoid.gamma, 0, col="red", lwd=2) 

# Which of the null values are smaller (more negative) than the data?
smallerNull <- g1_null<=balistoid.gamma
# How many TRUEs are there?
count<-sum(smallerNull)
count
# finally, what is the p-value?
mccr_pval <- (count+1)/(num_simulations+1)
mccr_pval

pb <- pureBirth(balistoid.times) #pb
aic.pb <- pb$aic#store aic score
bd <- bd(balistoid.times) # birth-death
aic.bd <- bd$aic #store aic
ddl <- DDL(balistoid.times) # logistic density dependent diversification
aic.ddl <- ddl$aic
ddx <- DDX(balistoid.times) # exponential density dependent diversification
aic.ddx <- ddx$aic
#the MCCR test above cannot discriminate between early rapid speciation and early slow extinction.
#these models allow those to scenarios to be compared 
spvar <- fitSPVAR(balistoid.times, init=c(2, .2, .1))#variable speciation: exponentially declining speciation, constant extinction rate
aic.spvar <- spvar$aic
exvar <- fitEXVAR(balistoid.times, init=c(.3, .01, 1))#variable extinction: constant speciation, exponentially increasing extinction rate
aic.exvar <- exvar$aic
both <- fitBOTHVAR(balistoid.times, init=c(.3, .5, .1, .5))#variable extinction and speciation
aic.both <- both$aic

#make a table and compare them all
aic.all <- cbind(aic.pb, aic.bd, aic.ddl, aic.ddx, aic.spvar, aic.exvar, aic.both)
foo <- function(x) x-x[which(x==min(x))] #quick function to compare all models to the model with lowest aic score
daic <- t(apply(aic.all, 1, foo))#apply that function to everything in aic.all
cat("Table of delta-aic values; zero - best model\n")
print(daic, digits=2)


######################################################################
#Additional simulations on gamma statistic and the Colless tree imbalance statistic
library(ape)
library(phytools)
numtrees <- 1000
trees <- pbtree(n = 50, nsim = numtrees, ape = F)

foo <- function(tree, metric = "colless") {
  if (metric == "colless") {
    #xx <- as.treeshape(x)  # convert to apTreeshape format
    colless(tree)
    #colless(xx, norm = "yule")  # calculate colless' metric
  } else if (metric == "gamma") {
    gammaStat(tree)
  } else stop("metric should be one of colless or gamma")
}

theme_myblank <- function() {
  stopifnot(require(ggplot2))
  theme_blank <- ggplot2::theme_blank
  ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.background = element_blank(), plot.background = element_blank(), 
                 axis.title.x = element_text(colour = NA), axis.title.y = element_blank(), 
                 axis.text.x = element_blank(), axis.text.y = element_blank(), axis.line = element_blank(), 
                 axis.ticks = element_blank())
}
library(plyr)
library(apTreeshape)
colless_df <- ldply(trees, foo, metric = "colless")  # calculate metric for each tree

#calculate cutoffs 
col_percent_low <- round(length(colless_df[colless_df$V1 < -0.7, "V1"])/numtrees, 2) * 100
col_percent_high <- round(length(colless_df[colless_df$V1 > 0.7, "V1"])/numtrees, 2) * 100

library(ggplot2)