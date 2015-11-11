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
par(mfrow = c(1,2))
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
nodelabels()
axisPhylo()

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

