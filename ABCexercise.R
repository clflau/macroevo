rm(list = ls())
#Approximate Bayesian Computation

#ABC rejection sampling algorithm
##a) Assume that N can be any value between 100 and 100,000. Draw 1 million values from the prior distribution of N. Let's use a uniform distribution for N. Assume that N~U[100,100000].
priorNN <- runif(1000000, min = 100, max = 100000)

##b) For each value of N, simulate a TMRCA. Note, you should use the same approach as on last week's assignment (i.e. Draw times from an exponential distribution). Second hint: The total number of Y chromosomes in the population is N, rather than 2N (because the Y chromosome is haploid). Adjust your rates of coalescence accordingly.
CoalesTime <- rep(NA, 1000000)
for (ii in 1:1000000){
  CoalesTime[ii] <- rexp(1, 1/(priorNN[ii]))
}
mean(CoalesTime)

##c) For each TMRCA, add a Poisson number of mutations with the appropriate mutation rate. Again, this will follow from what you did last week.
bp <- 100000
mu <- 1e-8*bp
MutaRate <- 2*mu*CoalesTime
NumbSNP <- rep(NA, 1000000)
for (ii in 1:1000000){
  NumbSNP[ii] <- rpois(1, MutaRate[ii])
}
mean(NumbSNP)

##d) We now need to decide which of the million draws from the prior give data that are "close" to the observed number of SNPs in the actual data and should be accepted. To do this, let's accept all values of N that give somewhere between 45-55 SNPs.
sims <- cbind(priorNN, NumbSNP)
acceptNN <- subset(sims, sims[, 2] >= 45 & sims[, 2] <= 55)
posteriorNN <- acceptNN[, 1]

##2) make density plots
plot(density(priorNN),col = 2, ylim = c(0, 2e-5), lwd = 1,xlab = "population size", main = "Prior and Posterior")  
lines(density(posteriorNN),col = 4,lwd = 1)
legend(0, 2e-5 ,c("prior","posterior"),lwd=1,col=c(2,4),cex=1)

##3) median posterior
median(posteriorNN)

##4) 95% credible interval
abline(v=quantile(posteriorNN,0.95),lty=2,lwd=1,col=1)

##6) Repeat your ABC analysis, but change the prior distribution of N to be U~[1000,1000000]. What is the mean, median, and 95% credible interval for the posterior distribution. How does this differ from what you computed in questions 3-4 for the original prior distribution? What does this tell you about the effect of the prior distribution for in Bayesian statistics?


#Exercise 3: ABC approach to estimate t_split

##a) Assume that tsplit can be any value between 50,000 and 1,000,000 generations. Draw 1 million values from the prior distribution of tsplit. Assume that tsplit ~U[50000,1000000].
PriorTT_split <- runif(1000000, min = 50000, max = 1000000)

##b) For each value of tsplit, simulate a TMRCA. Note, you should use the same approach as on the first part of this assignment. But, keep in mind that the coalescent time is a function of both the split time (tsplit, which is drawn from your prior) and the coalescent time in the ancestral population. 
NN <- 100000
CoalesTime <- rep(NA, 1000000)
for (ii in 1:1000000){
  CoalesTime[ii] <- rexp(1, 1/NN)
}
TMRCA <- CoalesTime + PriorTT_split

##c) For each TMRCA, add a Poisson number of mutations with the appropriate mutation rate. Again, this will follow from what you did last week.
bp <- 100000
mu <- 1e-8*bp
MutaRate <- 2*mu*TMRCA
NumbSNP <- rep(NA, 1000000)
for (ii in 1:1000000){
  NumbSNP[ii] <- rpois(1, MutaRate[ii])
}
mean(NumbSNP)

##d) We now need to decide which of the million draws from the prior generate data that are "close" to the observed number of SNPs in the actual data and should be accepted. To do this, let's accept all values of tsplit that give somewhere between 550 and 650 SNPs.
sims <- cbind(PriorTT_split, NumbSNP)
acceptTT_split <- subset(sims, sims[, 2] >= 550 & sims[, 2] <= 650)
posteriorTT_split <- acceptTT_split[, 1]

##Exercise 4: Density Plot
plot(density(PriorTT_split),col = 2, ylim = c(0, 8e-6), lwd = 1,xlab = "population size", main = "Prior and Posterior")  
lines(density(posteriorTT_split),col = 4,lwd = 1)
legend(8e5, 8e-6 ,c("prior","posterior"),lwd = 1, col=c(2,4),cex=1)

##Exercise 5: Median posterior
median(posteriorTT_split)

##Exercise 6: 95% credible interval
abline(v=quantile(posteriorTT_split,0.95),lty=2,lwd=1,col=1)





