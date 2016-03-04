set.seed(8063)

#aim 1
# background mutation rate, our null hypothesis
p0s = c(0.01,0.02,0.05)
# cnot3 mutation rate (% samples mutated) in ion torrent data, this is our alternative hypothesis
p1 = 0.184 

for (p0 in p0s) { 
# from sample size 1:100, I generated 1000 binomial random variable , each random variable represents number of mutated samples at a given sample size and p
k = matrix(NA, nrow=10000, ncol=100)
for(n in 1:100){
  k[,n] = rbinom(1000, n, p1)
}

# peform binomial test on simulated random variable and store resulting p-values in p.mat
p.mat  = matrix(NA, nrow=1000, ncol=100)
for(i in 1:1000){
  for(n in 1:100){
    p.mat[i,n] = binom.test(x=k[i,n], n=n, p0)$p.value
  }
}

# for each sample size, count proportion of random variables that have significant results
power = c()
power = apply(p.mat, 2, FUN=function(x){ mean(x<.05) })
png(paste("power-",p0,".png", width = 1024, height = 768, units = 'px'))
plot(power, xlab="samples", ylab="power")
dev.off()
print(paste0("p0=",p0))
print(min(which(power>.8)))
print(min(which(power>.9)))
}

#mutsig doesn't work
#p0 = 1 − (1 − μfg)(3L/4). We used L = 1,500, and fg = 3.9 (representing the 90th percentile of fg Lg/1,500 across the approximately 18,000 genes and fg = 1 for the 50th percentile gene)
#u=3 # background mutation frequency
# f=3.9
# L=1500
# 1-(1-u*f)^(3*L/4)


# % increase of tumor numbers due to 
#inc= c(1.05, 1.1, 1.15, 1.2, 1.25, 1.30)
#inc= c(1.15, 1.2, 1.30)
inc= c(1.1, 1.2, 1.3)
for ( i in inc) {
  muA=29
  muB=ceiling(muA*i)
  kappa=1
  sdA=10
  sdB=10
  alpha=0.05
  beta=0.1
  (nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
  cat(paste0(muB, "\t", ceiling(nA),"\n"))
  #z=(muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)
  #(Power=pnorm(z-qnorm(1-alpha)))
}

inc= c(0.9, 0.8, 0.7)
for ( i in inc) {
  muA=119
  muB=ceiling(muA*i)
  kappa=1
  sdA=31
  sdB=31
  alpha=0.05
  beta=0.1 
  (nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
  cat(paste0(muB, "\t", ceiling(nA),"\n"))
  #z=(muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)
  #(Power=pnorm(z-qnorm(1-alpha)))
}

# aim 3
delta1=.5
library(PROPER)
sim.opts.Bottomly = RNAseq.SimOptions.2grp(ngenes = 31600, p.DE=0.04, lOD="Bottomly", lBaselineExpr="Bottomly")
simres = runSims(Nreps = c(5,7,10,12,15), sim.opts=sim.opts.Bottomly, DEmethod="edgeR", nsims=20)
powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1, stratify.by="expr", delta=delta1)
summary.power(powers)
plot.Power(powers)

sim.opts.Bottomly2 = RNAseq.SimOptions.2grp(ngenes = 31600, p.DE=0.07, lOD="Bottomly", lBaselineExpr="Bottomly")
simres2 = runSims(Nreps = c(5,7,10,12,15), sim.opts=sim.opts.Bottomly2, DEmethod="edgeR", nsims=20)
powers2 = comparePower(simres2, alpha.type="fdr", alpha.nominal=0.1, stratify.by="expr", delta=delta1)
summary.power(powers2)
plot.Power(powers2)

sim.opts.Bottomly3 = RNAseq.SimOptions.2grp(ngenes = 31600, p.DE=0.1, lOD="Bottomly", lBaselineExpr="Bottomly")
simres3 = runSims(Nreps = c(5,7,10,12,15), sim.opts=sim.opts.Bottomly3, DEmethod="edgeR", nsims=20)
powers3 = comparePower(simres3, alpha.type="fdr", alpha.nominal=0.1, stratify.by="expr", delta=delta1)
summary.power(powers3)
plot.Power(powers3)

sim.opts.Bottomly4 = RNAseq.SimOptions.2grp(ngenes = 31600, p.DE=0.5, lOD="Bottomly", lBaselineExpr="Bottomly")
simres4 = runSims(Nreps = c(5,7,10), sim.opts=sim.opts.Bottomly4, DEmethod="edgeR", nsims=20)
powers4 = comparePower(simres4, alpha.type="fdr", alpha.nominal=0.1, stratify.by="expr", delta=delta1)
summary.power(powers4)
plot.Power(powers4)


sim.opts.Bottomly5 = RNAseq.SimOptions.2grp(ngenes = 31600, p.DE=0.07, seqDepth=100000000)
simres2 = runSims(Nreps = c(5,7,10,12,15), sim.opts=sim.opts.Bottomly2, DEmethod="edgeR", nsims=20)
powers2 = comparePower(simres2, alpha.type="fdr", alpha.nominal=0.1, stratify.by="expr", delta=delta1)

# cnot3 
controls = c(1984.2445,  1104.558,  605.266,	1122.448)
controls = c(1104.558,  605.266,  1122.448)
cases = c(803.225,  898.128,	981.642,	247.45)
outlier = function(data){
  iqr = IQR(data)
  lowerq = quantile(data)[2]
  upperq = quantile(data)[4]
  mild.threshold.upper = (iqr * 1.5) + upperq
  mild.threshold.lower = lowerq - (iqr * 1.5)
  print(mild.threshold.upper)
  print(mild.threshold.lower)
}
outlier(controls)
outlier(cases)


inc= c(.9, .8, .7, .6, .5)
for ( i in inc) {
  muA=mean(controls)
  muB=ceiling(muA*i)
  kappa=1
  sdA=sd(controls)
  #sdB=sd(cases)
  sdB = sdA
  alpha=0.05
  beta=0.20 
  (nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
  cat(paste0(muB, "\t", ceiling(nA),"\n"))
  #z=(muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)
  #(Power=pnorm(z-qnorm(1-alpha)))
}


###################### TEST ##########################
p_iot = 7/38
p_ill = 2/25

p0 = .01

n=100
sum(dbinom(0:n,n,p0))
plot(dbinom(0:n,n,p0))

# number of success required to reject null hypothesis for different sample size
1-pbinom(3,n,p0)  
1-pbinom(2,50,p0)
1-pbinom(1,20,p0)

# p alternative 
1-pbinom(24,100,.18)

for (n in 1:100) {
  k = rbinom(1000, n, p)
}

# attempt 2
set.seed(8063)

# background mutation rate, our null hypothesis
p0 = 0.01
# cnot3 mutation rate (% samples mutated) in ion torrent data, this is our alternative hypothesis
p1 = 0.184 

# from sample size 1:100, I generated 1000 binomial random variable , each random variable represents number of mutated samples at a given sample size and p
k = matrix(NA, nrow=1000, ncol=100)
for(n in 1:100){
  k[,n] = rbinom(1000, n, p1)
}

# peform binomial test on simulated random variable and store resulting p-values in p.mat
p.mat  = matrix(NA, nrow=1000, ncol=100)
for(i in 1:1000){
  for(n in 1:100){
    p.mat[i,n] = binom.test(x=k[i,n], n=n, p0)$p.value
  }
}

# for each sample size, count proportion of random variables that have significant results
power = c()
power = apply(p.mat, 2, FUN=function(x){ mean(x<.05) })
png("power.png")
plot(power, xlab="samples", ylab="power")
dev.off()
min(which(power>.8))
min(which(power>.9))

# example of bionomial test, h0 = fair dice of 1/6 chance to roll a 6 
# significance of 51 "6" in 235 roll dice
# same as p(rolling 51 or greater "6") , the probabliy is less than .05, very unlikely,  so we reject null of fair dice and accept that it is biased.
binom.test(51,235,(1/6))
binom.test(51,235,(1/6), alternative="greater")
sum(dbinom(51:235,235,(1/6)))

# how many samples do I need to see mutated sample rate(p) at 80% power
# y axis= power, x axis = number of sample, at a mutated sample rate (p)
p = .3 # probability of a sample mutated at gene x
n=20  # number trials
size=100 # number of samples
k = rbinom(n, size, p) # returns number of mutated samples given size , p
# random variable x = number of mutated samples
pbinom(k, n, p) # p(x<=k) given n and p

# this will give a bionmial distribution, which half the time it exceeds p 
n=100000; hist(rbinom(n, 100, p)/100)

# number of success in trials (size) with probability = p
x = rbinom(n, size, p)
plot(pbinom(x, size, p_iot))

nn   <- 1:100      # sample sizes
pow  <- 1-pbinom(1, nn, .)  # 1 - p(x=0) at different sample size, which is p(x>0) at different sample size (1 or more success)
tStr <- expression(paste("Power for ", X>0, " given ", p[1]==0.3))
plot(nn, pow, type="l", xaxs="i", xlab="sample size", ylab="power",
     lwd=2, col="blue", main=tStr, cex.lab=1.4, cex.main=1.4)

# using library
library(pwr)
pwr.p.test(sig.level=0.05, power=.8, h = ES.h(.3, .01), alt="greater", n = NULL)


# example
n  <- 500                 # sample size
p1 <- 0.001               # success probability under alternative hypothesis
cc <- 1                   # threshold
dbinom(cc:n, n, p1)
sum(dbinom(cc:n, n, p1))  # power: probability for cc or more successes given p1
# same thing as saying culumlative distribution of P(1<=k<=n)
# which is p(x<=n) - p(x=0)
pbinom(n,n,p1) - pbinom(0,n,p1)

nn   <- 1:50                 # sample sizes
pow  <- 1-pbinom(cc-1, nn, p1)  # 1 - p(x=0) at different sample size, which is p(x>0) at different sample size (1 or more success)
tStr <- expression(paste("Power for ", X>0, " given ", p[1]==0.3))
plot(nn, pow, type="l", xaxs="i", xlab="sample size", ylab="power",
     lwd=2, col="blue", main=tStr, cex.lab=1.4, cex.main=1.4)


# coin example
k=9
#sum(dbinom(0:k, 10, .5)) = pbinom(k,10,.5)
sum(dbinom(0:k, 10, .5)) 
pbinom(k,10,.5)

# p(X<=k)
pbinom(9,10,.5)

# get it right 12 times or more
1-pbinom(11,16,.5)
1-pbinom(11,16,.75)
sum(dbinom(12:16,16,.75))
plot(dbinom(1:16,16,.75))


