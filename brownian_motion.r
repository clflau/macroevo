#simulate brownian motion 
##Exercise 1
par(mfrow = c(1,1))
displace <- rnorm(1000)  #series of draws from normal distribution
plot(displace)
library(geiger)
x <- cumsum(displace)    #cumulative sum through time
plot(x, type = "l", xlab = "Time", ylab = "Trait Value")

cols<-rainbow(50) # put 50 samples from the rainbow pallette into a new vector called cols
plot(cumsum(rnorm(1000)), type="l", ylim=c(-100,100))
for(i in 1:length(cols)) 
  lines(cumsum(rnorm(1000)), col=cols[i]) # simple for loop for each of the 50 colour samples we draw a line on the plot that is the cumulative summation of 100 draws from a normal distribution. 

##Exercise 2
layout(matrix(1:2, 1, 2))
plot(cumsum(rnorm(50)), type="l", ylim=c(-50,50), main="Effect of time under BM")
for(i in 1:50)
  lines(cumsum(rnorm(50)), col="blue")
plot(cumsum(rnorm(100)), type="l", ylim=c(-50,50))
for(i in 1:50)
  lines(cumsum(rnorm(100)), col="red")

##Exercise 3
par(mfcol = c(1,1))
plot(cumsum(rnorm(100)), type="l", ylim=c(-80,80), ylab="Trait", col="blue")
for(i in 1:20) lines(cumsum(rnorm(100)), col="blue")
for(i in 1:20) lines(cumsum(rnorm(100, sd=5)), col="green")  # sigma=5, sigma^2=25

##Exercise 4
tree<-read.tree(text="((grizzlybear:1, polarbear:1):4, spectacledbear:5);")
plot(tree)
nodelabels() #adds node labels
axisPhylo() #adds branchlength axis

rootMass <- 250 # size of ancestor
sigmasq = 2.5 # Brownian rate
time = 5 # 5 million years of independent evolution from the root
sd <- sqrt(time * sigmasq) # Brownian evolution is proportional to rate X time)
specbearDeltaMass <- rnorm(1, mean = 0, sd = sd)
specbearMass <- rootMass + specbearDeltaMass
specbearMass

#another way to do this using the tree structure itself to supply the time argument
specbearDeltaMass<-rnorm(1, mean=0, sd=sqrt(sigmasq*tree$edge.length[1])) # in R, tree tip labels are numbered from top to bottomw, so in our tree spectacled bear = 1, polarbear = 2, and grizzly = 3
specbearMass <- rootMass + specbearDeltaMass
specbearMass
tree$edge.length[1]

rootMass <- 250 # size of ancestor
sigmasq = 2.5 # Brownian rate
#spectacledbear
time = 5 # 5 million years of independent evolution from the root
sd <- sqrt(time * sigmasq) # Brownian evolution is proportional to rate X time)
specbearDeltaMass <- rnorm(1, mean = 0, sd = sd)
specbearMass <- rootMass + specbearDeltaMass
cat("spectacled bear mass\n")
specbearMass
#polarbear + grizzlybear
time = 4
sd <- sqrt(time * sigmasq)
polgrizbearMass <- rootMass + rnorm(1, mean = 0, sd = sd)
#polarbear
time = 1
sd <- sqrt(time * sigmasq)
polbearMass <- polgrizbearMass + rnorm(1, mean = 0, sd = sd)
cat("polar bear mass\n")
polbearMass
#grizzlybear
time = 1
sd <- sqrt(time * sigmasq)
grizbearMass <- polgrizbearMass + rnorm(1, mean = 0, sd = sd)
cat("grizzly bear mass\n")
grizbearMass

#################################################################
#Trends
bmtrend <- function(time, mu, stdev){ # a function that requires mu which is the trend parameter showing how the mean of the normal distribution changes over time, time the length of the simulation and stdev the standard deviation of the normal distribution, it assumes the starting ancestral trait value is 0
  displace < -vector(, length=time) # set up an empty vector of the same length as time
  for(i in 1:time){
    displace[i] <- rnorm(1, mean=mu*i, sd=stdev) # sample from rnorm with a mean with a mean of mu multiplied by the iteration
  }
  traj <- cumsum(displace) # calculates the trajectory over time        
  return(traj) 
}
#Now we are going to run this bmtrend function 10 times using the replicate function, once with a mu of 0.02 and the other with mu=-0.02
par(mfrow = c(1,1))
bmtrend0.02 <- replicate(10, bmtrend(100, 0.02, 1))
bmtrendneg0.02 <- replicate(10, bmtrend(100, -0.02, 1))

cols1 <- rainbow(10) # sample 10 colours from the rainbow palette
cols2 <- grey(0:10/10) # sample 10 colours from the grey/gray level specification - it has a different syntax from the other palettes like rainbow, topo, terrain, heat.colors etc.

plot(bmtrend0.02[, 1], xlim=c(0,100), ylim=c(-150, 150), type="l", col=cols1[1], xlab="Time", ylab="Trait")
for(i in 2:10) 
  lines(bmtrend0.02[, i], col=cols1[i])
for(i in 1:10) 
  lines(bmtrendneg0.02[, i], col=cols2[i])
# simulate 15 runs with mu of 0.5 and -0.5 for 100my
par(mfrow = c(1,1))
bmtrend0.5 <- replicate(15, bmtrend(100, 0.5, 1))
bmtrendneg0.5 <- replicate(15, bmtrend(100, -0.5, 1))

cols1 <- rainbow(10) # sample 10 colours from the rainbow palette
cols2 <- grey(0:10/10) # sample 10 colours from the grey/gray level specification - it has a different syntax from the other palettes like rainbow, topo, terrain, heat.colors etc.

plot(bmtrend0.5[, 1], xlim=c(0,100), ylim=c(-2500, 2500), type="l", col=cols1[1], xlab="Time", ylab="Trait")
for(i in 2:10) 
  lines(bmtrend0.5[, i], col=cols1[i])
for(i in 1:10) 
  lines(bmtrendneg0.5[, i], col=cols2[i])

################################################################
#OU model

alpha<-0.4 #optimal value
theta<-1 #strength of pull to optimal
sigma<-0.05 #sigma^2 the stochastic (Brownian) motion parameter
x0<-0

OUalpha0.1_sample1<- x0 + alpha*(theta-x0)+sigma*(rnorm(1,mean=0))

# the next sample would be

OUalpha0.1_sample2<- OUalpha0.1_sample1 + alpha*(theta-OUalpha0.1_sample1)+sigma*(rnorm(1,mean=0))
OUalpha0.1_sample3<- OUalpha0.1_sample2 + alpha*(theta-OUalpha0.1_sample2)+sigma*(rnorm(1,mean=0))
plot(x=0:3, y=c(1,OUalpha0.1_sample1, OUalpha0.1_sample2, OUalpha0.1_sample3), type="l", main="Alpha=0.4, Theta=1, Sigma=0.05, x0=1")

#OU function
ornstein_uhlenbecksim <- function(n,theta,alpha,sigma2,x0){# n is the number of time units to simulate for, theta is the optima, alpha is the pull towards the optima, sigma is the rate and X0 is the starting trait value
  dw  <- rnorm(n, 0)
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma2*dw[i-1]
  }
  return(x);
}

par(mfrow = c(1, 2))
alpha<-0.4
theta<-1
sigma<-0.05
x0<-0
time = 100

alpha0.4sigma0.05theta1<-replicate(10, ornstein_uhlenbecksim(time, theta, alpha, sigma, x0)) # replicate() just runs it 10 times for us to illustrate the variance between runs

colr<-topo.colors(9)

#10 replicates under OU
plot(alpha0.4sigma0.05theta1[,1], type="l", main="OU alpha=0.4 sigma=0.05", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.4sigma0.05theta1)) lines(alpha0.4sigma0.05theta1[,i], type="l",col=colr[i-1] )

#10 replicates under BM
cols<-rainbow(10)
BMsigma0.5 <-cumsum(rnorm(100, mean = 0, sd = sigma))
plot(BMsigma0.5, ylim = c(-2, 2), type = "l", main = "BM alpha=0.4 sigma=0.05")
for(i in 1:length(cols))
  lines(cumsum(rnorm(100, mean = 0, sd = sigma)), col=cols[i])
  
###Several examples of OU
alpha1sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 1, 0.01, 0)) # replicate just runs it 10 times for us to illustrate the variance between runs

colr<-topo.colors(9)

plot(alpha1sigma0.01theta1[,1], type="l", main="alpha1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.01theta1)) lines(alpha1sigma0.01theta1[,i], type="l",col=colr[i-1] )

# play with different parameters.

alpha0.5sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.5, 0.01, 0))
alpha0.1sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.1, 0.01, 0))
alpha0sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0, 0.01, 0))
alpha1sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 1, 0.1, 0))
alpha0.5sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.5, 0.1, 0))
alpha0.1sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.1, 0.1, 0))
alpha0sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0, 0.1, 0))

par(mfrow=c(2,2))

plot(alpha1sigma0.01theta1[,1], type="l", main="alpha1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.01theta1)) lines(alpha1sigma0.01theta1[,i], type="l",col=colr[i-1] )
plot(alpha1sigma0.1theta1[,1], type="l", main="alpha1 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.1theta1)) lines(alpha1sigma0.1theta1[,i], type="l",col=colr[i-1] )
plot(alpha0.5sigma0.01theta1[,1], type="l", main="alpha0.5 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.5sigma0.01theta1)) lines(alpha0.5sigma0.01theta1[,i], type="l",col=colr[i-1] )
plot(alpha0.5sigma0.1theta1[,1], type="l", main="alpha0.5 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.5sigma0.1theta1)) lines(alpha0.5sigma0.1theta1[,i], type="l", col=colr[i-1] )

plot(alpha0.1sigma0.01theta1[,1], type="l", main="alpha0.1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.1sigma0.01theta1)) lines(alpha0.1sigma0.01theta1[,i], type="l", col=colr[i-1] )
plot(alpha0.1sigma0.1theta1[,1], type="l", main="alpha0.1 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.1sigma0.1theta1)) lines(alpha0.1sigma0.1theta1[,i], type="l", col=colr[i-1] )
plot(alpha0sigma0.01theta1[,1], type="l", main="alpha0 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0sigma0.01theta1)) lines(alpha0sigma0.01theta1[,i], type="l", col=colr[i-1] )
plot(alpha0sigma0.1theta1[,1], type="l", main="alpha0 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0sigma0.1theta1)) lines(alpha0sigma0.1theta1[,i], type="l", col=colr[i-1] )

######################################################################
#Model fitting

library(ape)
library(geiger)
library(parallel)

tre <- read.tree("Haemulidae4models.tre")
dat<-read.csv("Haemulidae4models.csv", stringsAsFactors=FALSE) # Contains data on standard length, eye width and buccal width
head(dat)
trait <- dat[,2] # Selecting a trait 
names(trait) <- tre$tip.label # Assigning the appropriate tip labels to the trait object
names(trait) <- dat$taxon

#fitting BM model with geiger function fitContinuous
BM.model <- fitContinuous(tre, trait, model="BM")
BM.model
BM.model$bnd
BM.model$opt # The ouput includes info on sigma squared and the root state (z0)
BM.model$opt$lnL
BM.model$opt$aicc
#fitting OU model
OU.model <- fitContinuous(tre, trait, model="OU") #could run into problems with bounds too narrow
OU.model$bnd
OU.model <- fitContinuous(tre, trait, model="OU", bounds=list(alpha=c(0, 150))) #change bounds

# How does the likelihood compare to the BM model?
BM.model$opt$lnL
OU.model$opt$lnL

# It looks like the OU model has the higher likelihood...

# One way to assess the difference between these two models more formally is the likelihood ratio test.
# BM is a special case of OU (for alpha=0), so the models are nested and the LR test is fine.
# Note that for thorough tests you should account for phylogenetic uncertainty and perform parameteric bootstrapping.
delta.BM.OU <- 1-pchisq(2*(OU.model$opt$lnL - BM.model$opt$lnL),1)
delta.BM.OU

#AIC and AICc
BM.model$opt$aicc # Higher AICc score is worse!
OU.model$opt$aicc # Lower AICc score is better!
all.aicc <- c(BM.model$opt$aicc, OU.model$opt$aicc)
delta.aicc <- all.aicc - min(all.aicc) # This may seem a bit circumstantial for just two models, but will pay off for additional models.
delta.aicc # delta AIC or AICC scores >2 are ususally considered to provide positive support
#akaike weights: AIC and AICc scores can also be expressed as Akaike weights, representing relative likelihood of the model (=exp( -0.5 * deltaAIC score for that model). Akaike weight for a given model are the rel. likelihood of the model divided by sum of all relative likelihoods across all models:
rel.L <- exp(-delta.aicc*0.5)
rel.L #relative likelihood
AK.weights <- rel.L/sum(rel.L)
AK.weights

#eye width
trait <- dat[,3]
names(trait) <- dat$taxon
BM.model <- fitContinuous(tre, trait, model="BM")
OU.model <- fitContinuous(tre, trait, model="OU")
EB.model <- fitContinuous(tre, trait, model="EB", bounds=list(a=c(-10, 0)))
EB.model$bnd
trend.model <- fitContinuous(tre, trait, model="trend")

#Report the likelihoods, AICcs, and delta AICs in a table
model <- c("BM", "OU", "EB", "trend") 
likelihood <- c(BM.model$opt$lnL, OU.model$opt$lnL, EB.model$opt$lnL, trend.model$opt$lnL) 
AICc <- c(BM.model$opt$aicc, OU.model$opt$aicc, EB.model$opt$aicc, trend.model$opt$aicc)
all.aicc <- c(BM.model$opt$aicc, OU.model$opt$aicc, EB.model$opt$aicc, trend.model$opt$aicc)
delta.aicc <- all.aicc - min(all.aicc)
rel.L <- exp(-delta.aicc*0.5)
AK.weights <- rel.L/sum(rel.L)
df <- data.frame(model, likelihood, AICc, delta.aicc, AK.weights) 
cat("Summary Table")
print(df)

