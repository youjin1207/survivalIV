library(MASS)
library(parallel)
library(doParallel)
library(ranger)


source("Code/core/nuisance_binary_survival_para.R")
source("Code/core/binary_survival_para_sim.R")
source("Code/core/binary_survival_ranger_sim.R")

nsplits = 2

### generate simulation data
N = 1000
X = mvrnorm(N, rep(0,5), Sigma = diag(5))
deltas = c(-0.5, -0.5, 1, -1, -0.7)
alphas = c(0.1, 0.1, -0.2, -0.2, 0.3)*5
Z = rbinom(N, 1, plogis(as.numeric(X %*% deltas))) # define as a function of X
U = rnorm(N, 0, 1)
A = rbinom(N, 1, plogis(as.numeric(X %*% alphas) + 2*Z + 1.0*U - 0.1)) ## U : unmeasured confounder

betas = c(0.5, 0.5, -0.5, -0.5, -0.5)
AtoT = 1.5
UtoT = 2.5
gammas = -c(-0.3, -0.3, 0.3, 0.3, 0.3)*10
ZtoC = +0.0
AtoC = +0.5
t = c(1:180); t2 = t^2
tau = 180
coeff.Lambda = c(0.00, 0.0005, 0.0003)
Lambda.conti = coeff.Lambda[1] + coeff.Lambda[2]*t + coeff.Lambda[3]*t2 
Lambda.conti.censor = (coeff.Lambda[1]/2) + (coeff.Lambda[2]/4)*t + 2*coeff.Lambda[3]*t2 


## generate survival outcome under Cox model + Censoring
t = c(1:180); t2 = t^2
tau = 180
coeff.Lambda = c(0.00, 0.0005, 0.0003)
coeff.Lambda2 = c(0.00, 0.0025, 0.0002)
Lambda.conti = coeff.Lambda[1] + coeff.Lambda[2]*t + coeff.Lambda[3]*t2 
## cumulative hazards for T and C ##
print.Lambda = function(time){
  coeff.Lambda[1] + coeff.Lambda[2]*time + coeff.Lambda[3]*(time)^2
}

inverse.Lambda = function(input){
  (1/(coeff.Lambda[3]))*(-coeff.Lambda[2]/2 +  sqrt(coeff.Lambda[3]*input + coeff.Lambda[2]^2/4 -coeff.Lambda[1]*coeff.Lambda[3]))
}

inverse.Lambda2 = function(input){
  (1/(coeff.Lambda2[3]))*(-coeff.Lambda2[2]/2 +  sqrt(coeff.Lambda2[3]*input + coeff.Lambda2[2]^2/4 -coeff.Lambda2[1]*coeff.Lambda2[3]))
}


U.cf = runif(N, 0, 1)
obs.T = inverse.Lambda((-log(U.cf)*exp(-as.numeric(as.numeric(X %*% betas) + A*AtoT + U*UtoT))))

obs.T1 = inverse.Lambda((-log(U.cf)*exp(-as.numeric(as.numeric(X %*% betas) + 1*AtoT + U*UtoT))))
obs.T0 = inverse.Lambda((-log(U.cf)*exp(-as.numeric(as.numeric(X %*% betas) + 0*AtoT + U*UtoT))))

obs.C = obs.T + runif(N, -10, 50) # T-C is randomly distributed

time.list = seq(0, 30, 1)
obs.Y = pmin(obs.T, obs.C)
obs.T = as.integer(obs.T)
obs.C = as.integer(obs.C)
obs.Y = as.integer(obs.Y)
obs.T1= as.integer(obs.T1)
obs.T0 = as.integer(obs.T0)

surv1.obs = surv0.obs = c()
for(r in 1:length(time.list)){
  surv1.obs[r] = mean(obs.T1 > time.list[r])
  surv0.obs[r] = mean(obs.T0 > time.list[r])
}

## transform the first four baseline covariates
W = X
W[,1] = exp(X[,1]/2)
W[,2] = X[,2] / (1+exp(X[,1])) + 10
W[,3] = (X[,1]*X[,3]/25 + 0.6 )^3
W[,4] = (X[,2] + X[,4] + 20)^2

time.list = seq(0, 30, 1)
dat = data.frame(obs.T = obs.T, obs.C = obs.C, obs.Y = obs.Y, A = A, X = X, U = U, Z = Z, W = W)
dat$R = as.integer(dat$obs.T < dat$obs.C)

### parametric estimator
## (1) IF estimator
result.if = binary_survival_para_if(dat = dat, time.list = time.list, survival = "Cox")

result.if_omega_delta = binary_survival_para_if(dat = dat, 
	time.list = time.list, misspecify = c("censoring", "instrument"), survival = "Cox")

result.if_pi_mu = binary_survival_para_if(dat = dat, time.list = time.list,
                             misspecify = c("outcome", "treatment"), survival = "Cox")

result.if_pi_omega = binary_survival_para_if(dat = dat, time.list = time.list,
                            misspecify = c("treatment", "censoring"), survival = "Cox")

## (2) IF-hazard estimator
result.if_hazard = binary_survival_para_if_hazard(dat = dat, 
                            time.list = time.list, survival = "Cox")

result.if_hazard_omega_delta = binary_survival_para_if_hazard(dat = dat, 
                            time.list = time.list, misspecify = c("censoring", "instrument"),
                            survival = "Cox")

result.if_hazard_pi_mu = binary_survival_para_if_hazard(dat = dat, 
                            time.list = time.list, misspecify = c("hazard", "treatment"),
                            survival = "Cox")

result.if_hazard_pi_omega = binary_survival_para_if_hazard(dat = dat,
                            time.list = time.list, misspecify = c("treatment", "censoring"),
                            survival = "Cox")

## (3) naive-hazard estimator
naive.if =  binary_survival_para_naive(dat = dat, 
                            time.list = time.list, survival = "Cox")

naive.if_omega_delta = binary_survival_para_naive(dat = dat, 
                            time.list = time.list, misspecify = c("censoring", "instrument"),
                            survival = "Cox")

naive.if_pi_mu = binary_survival_para_naive(dat = dat, 
                            time.list = time.list, misspecify = c("hazard", "treatment"),  survival = "Cox")

naive.if_pi_omega = binary_survival_para_naive(dat = dat, 
                            time.list = time.list, misspecify = c("treatment", "censoring"),  survival = "Cox")



### nonparametric estimator
x.trt = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5, Z = dat$Z)
x.out = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5, Z = dat$Z, A = dat$A)
x.miss = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5, Z = dat$Z, A = dat$A)
x.instrument = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                          X4 = dat$X.4, X5 = dat$X.5)
## naive estimator does not use the instrument Z
x.trt.naive = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5)
x.out.naive = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5, A = dat$A)
x.miss.naive = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5, A = dat$A)
## misspecified covariates set
w.trt = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5, Z = dat$Z)
w.out = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5, Z = dat$Z, A = dat$A)
w.miss = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5, Z = dat$Z, A = dat$A)
w.instrument = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                          W4 = dat$W.4, W5 = dat$W.5)

w.trt.naive = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5)
w.out.naive = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5, A = dat$A)
w.miss.naive = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5, A = dat$A)

## (1) IF estimator
ranger.result.if = binary_survival_ranger_if(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, time.list = time.list, nsplits = nsplits)

ranger.result.if_omega_delta = binary_survival_ranger_if(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, 
	time.list = time.list, nsplits = nsplits, misspecify = c("censoring", "instrument"))

ranger.result.if_pi_mu = binary_survival_ranger_if(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, time.list = time.list, nsplits = nsplits,
                             misspecify = c("outcome", "treatment"))

ranger.result.if_pi_omega = binary_survival_ranger_if(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, time.list = time.list, nsplits = nsplits,
                            misspecify = c("treatment", "censoring"))

## (2) IF-hazard estimator
ranger.result.if_hazard = binary_survival_ranger_if_hazard(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list, nsplits = nsplits)

ranger.result.if_hazard_omega_delta = binary_survival_ranger_if_hazard(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, nsplits = nsplits,
                            time.list = time.list, misspecify = c("censoring", "instrument"))

ranger.result.if_hazard_pi_mu = binary_survival_ranger_if_hazard(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, nsplits = nsplits,
                            time.list = time.list, misspecify = c("hazard", "treatment"))

ranger.result.if_hazard_pi_omega = binary_survival_ranger_if_hazard(dat = dat, x.trt = x.trt, x.out = x.out, 
	x.miss = x.miss, x.instrument = x.instrument, 
	w.trt = w.trt, w.out = w.out, 
	w.miss = w.miss, w.instrument = w.instrument, nsplits = nsplits,
                            time.list = time.list, misspecify = c("treatment", "censoring"))

## (3) naive-hazard estimator
ranger.naive.if =  binary_survival_ranger_naive(dat = dat, x.trt = x.trt.naive, x.out = x.out.naive, 
	x.miss = x.miss.naive, x.instrument = NULL, 
	w.trt = w.trt.naive, w.out = w.out.naive, 
	w.miss = w.miss.naive, w.instrument = NULL,
                            time.list = time.list, nsplits = nsplits)

ranger.naive.if_omega_delta = binary_survival_ranger_naive(dat = dat,x.trt = x.trt.naive, x.out = x.out.naive, 
	x.miss = x.miss.naive, x.instrument = NULL, 
	w.trt = w.trt.naive, w.out = w.out.naive, 
	w.miss = w.miss.naive, w.instrument = NULL,
                            time.list = time.list, nsplits = nsplits, misspecify = c("censoring", "instrument"))

ranger.naive.if_pi_mu = binary_survival_ranger_naive(dat = dat, x.trt = x.trt.naive, x.out = x.out.naive, 
	x.miss = x.miss.naive, x.instrument = NULL, 
	w.trt = w.trt.naive, w.out = w.out.naive, 
	w.miss = w.miss.naive, w.instrument = NULL,
                            time.list = time.list, nsplits = nsplits, misspecify = c("hazard", "treatment"))

ranger.naive.if_pi_omega = binary_survival_ranger_naive(dat = dat, x.trt = x.trt.naive, x.out = x.out.naive, 
	x.miss = x.miss.naive, x.instrument = NULL, 
	w.trt = w.trt.naive, w.out = w.out.naive, 
	w.miss = w.miss.naive, w.instrument = NULL, 
                            time.list = time.list, nsplits = nsplits, misspecify = c("treatment", "censoring"))


