rm(list = ls())

##Genetic Drift exercises

#expected number of success
mean(rbinom(10000, 10, 0.1))



#genetic drift
pp <- 0.1
NN <- 10 # number of diploid indiv; 2*NN = number of chromosome
count <- rbinom(1, 2*NN, pp)
count
count/(2*NN)


#function to simulate genetic drift (without mutations)
GenDriftFunc <- function(LL, NN, p0, TT){
  #matrix of allele freq of each SNP (col) in each gen (row)
  freqs <- matrix(nrow = TT, ncol = LL)
  #fill first row with p0
  freqs[1, ] <- p0
  #for loop to iterate rbinom for each subsequent row
  for (tt in 2:TT){
    for(ll in 1:LL){
        freqs[tt, ll] <- rbinom(1, 2*NN, freqs[tt-1, ll])/(2*NN)
    }
  }
  return(freqs)
}

#use function
NN=100
LL=1000
TT=10000
p0=0.1
sim1 <- GenDriftFunc(LL, NN, p0, TT)
#extinct?
sum(sim1[10000, ] == 0)
#fixation?
sum(sim1[10000, ] == 1)
#agree?
sum(sim1[10000, ] == 1)/1000
#plot
xx <- sim1[ , 1]
plot(xx, type = "l", ylim = c(0, 1), xlim = c(0, 2000), main = "N=100", xlab = "Generation", ylab = "Frequency")
for(ii in 2:100){
    lines(sim1[ , ii])
}

#repeat
NN=100
LL=1000
TT=10000
p0=0.6
sim2 <- GenDriftFunc(LL, NN, p0, TT)
#extinct?
sum(sim2[10000, ] == 0)
#fixation?
sum(sim2[10000, ] == 1)


# N=10
NN=10
LL=1000
TT=10000
p0=0.1
sim3 <- GenDriftFunc(LL, NN, p0, TT)

# N=500
NN=500
LL=1000
TT=10000
p0=0.1
sim4 <- GenDriftFunc(LL, NN, p0, TT)

# N=1000
NN=1000
LL=1000
TT=10000
p0=0.1
sim5 <- GenDriftFunc(LL, NN, p0, TT)

#plot sim1, sim3, sim4, sim5
#100=red; 10=green; 500=blue; 1000=cyan
plot(sim1[ , 1], type = "l", ylim = c(0, 1), xlim = c(0, 6000), main = "N=100", xlab = "Generation", ylab = "Frequency", col = 2, lwd = 0.1)
for(ii in 2:100){
  lines(sim1[ , ii], col = 2, lwd = 0.01)
  adjustcolor(col = 2, alpha.f = 0.7)
}
for(ii in 1:100){
  lines(sim3[ , ii], col = 3, lwd = 0.01)
  adjustcolor(col = 3, alpha.f = 0.9)
}
for(ii in 1:100){
  lines(sim4[ , ii], col = 4, lwd = 0.01)
  adjustcolor(col = 4, alpha.f = 0.3)
}
for(ii in 1:100){
  lines(sim5[ , ii], col = 5, lwd = 0.01)
  adjustcolor(col = 5, alpha.f = 0.1)
}


#fixation?
fix1 <- sum(sim1[10000, ] == 1)/1000
fix3 <- sum(sim3[10000, ] == 1)/1000
fix4 <- sum(sim4[10000, ] == 1)/1000
fix5 <- sum(sim5[10000, ] == 1)/1000
popsize <- c(10, 100, 500, 1000)
fixedproportion <- c(fix3, fix1, fix4, fix5)
table <- rbind(popsize, fixedproportion)
table


################################################################################
#Coalescence

#Draw 10 values from an exp dist with rate 2
rexp(10, 2)
mean(rexp(10000, 2))

#10,000 sims
NN <- 10000
CoalesTime <- rexp(10000, 1/(2*NN))
mean(CoalesTime)

#density plot
hist(CoalesTime)

sd(CoalesTime)


#########################
#poisson distribution
rpois(10, 2)
mean(rpois(10000, 2))

#add mutations, 10,000 sims, 2N=20,000
mu <- 1e-4
NN <- 10000
CoalesTime <- rexp(10000, 1/(2*NN))
MutaRate <- 2*mu*CoalesTime
<<<<<<< HEAD
NumbSNP <- rep(NA, 10000)
=======
mean(MutaRate)
NumbSNP <- rep(NA, NN)
>>>>>>> origin/master
for (ii in 1:10000){
  NumbSNP[ii] <- rpois(1, MutaRate[ii])
}
mean(NumbSNP)

rm(list = ls())

#theta=4*N*mu
theta <- 4*NN*mu
theta

#density plot
hist(NumbSNP)


