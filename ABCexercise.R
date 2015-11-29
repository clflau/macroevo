#Approximate Bayesian Computation

#ABC rejection sampling algorithm
##a) Assume that N can be any value between 100 and 100,000. Draw 1 million values from the prior distribution of N. Let's use a uniform distribution for N. Assume that N~U[100,100000].
NN <- runif(1000000, min = 100, max = 100000)

##b) For each value of N, simulate a TMRCA. Note, you should use the same approach as on last week's assignment (i.e. Draw times from an exponential distribution). Second hint: The total number of Y chromosomes in the population is N, rather than 2N (because the Y chromosome is haploid). Adjust your rates of coalescence accordingly.
bp <- 100000
mu <- 1e-8*bp
for (ii in 1:10000){
  CoalesTime[ii] <- rexp(1, 1/(NN[ii]))
}
mean(CoalesTime)

##c) For each TMRCA, add a Poisson number of mutations with the appropriate mutation rate. Again, this will follow from what you did last week.
MutaRate <- mu*CoalesTime
NumbSNP <- rep(NA, 10000)
for (ii in 1:10000){
  NumbSNP[ii] <- rpois(1, MutaRate[ii])
}
mean(NumbSNP)

##d) We now need to decide which of the million draws from the prior give data that are "close" to the observed number of SNPs in the actual data and should be accepted. To do this, let's accept all values of N that give somewhere between 45-55 SNPs.
lim45 <- NumbSNP >= 45





