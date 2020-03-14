library(MASS)
library(ranger)
library(survival)
library(ahaz)

source("Code/core/nuisance_continuous_survival_para.R")
source("Code/core/continuous_survival_para_sim.R")
source("Code/core/continuous_survival_ranger_sim.R")

### generate the simulation data
N=1000
kappa = 2
X = mvrnorm(N, rep(0,5), Sigma = diag(5))
deltas = c(-0.5, -0.5, 1, -1, -0.7)
alphas = c(0.1, 0.1, -0.2, -0.2, 0.3)*5
Z = rnorm(N, as.numeric(X %*% deltas), 2) 
U = rnorm(N, 0, 1)
A = rbinom(N, 1, plogis(as.numeric(X %*% alphas) + 2*Z + 0.3*U + 0.1)) ## U : unmeasured confounder


betas = c(0.5, 0.5, -0.5, -0.5, -0.5)/50
AtoT = 1.5/50
UtoT = 0.5/50
gammas = c(-0.3, -0.3, 0.3, 0.3, 0.3)
ZtoC = 0.5
AtoC = -0.5
t = c(1:180); t2 = t^2
tau = 180
coeff.Lambda = c(0.00, 0.0005, 0.0003)
Lambda.conti = coeff.Lambda[1] + coeff.Lambda[2]*t + coeff.Lambda[3]*t2 
Lambda.conti.censor = (coeff.Lambda[1]/2) + (coeff.Lambda[2]/4)*t + 2*coeff.Lambda[3]*t2 

## cumulative hazards for T and C ##
print.Lambda = function(time){
coeff.Lambda[1] + coeff.Lambda[2]*time + coeff.Lambda[3]*(time)^2
}
print.Censor = function(time){
coeff.Lambda[1]/2 + (coeff.Lambda[2]/4)*time + (coeff.Lambda[3]/2)*(time)^2
}

inverse.Lambda = function(input, A){
(1/(coeff.Lambda[3]))*(-(coeff.Lambda[2] + as.numeric(X %*% betas) + A*AtoT + U*UtoT)/2 +  
                       sqrt(coeff.Lambda[3]*input + (coeff.Lambda[2] + as.numeric(X %*% betas) + A*AtoT + U*UtoT)^2/4 -coeff.Lambda[1]*coeff.Lambda[3]))
}


U.cf = runif(N, 0, 1)
obs.T1 = inverse.Lambda((-log(U.cf)), 1)
obs.T0 = inverse.Lambda((-log(U.cf)), 0)
obs.T = inverse.Lambda((-log(U.cf)), A)
time.list = seq(0, 30, 1)

## transform the first four baseline covariates
W = X
W[,1] = exp(X[,1]/2)
W[,2] = X[,2] / (1+exp(X[,1])) + 10
W[,3] = (X[,1]*X[,3]/25 + 0.6 )^3
W[,4] = (X[,2] + X[,4] + 20)^2


dat = data.frame(obs.Y = obs.T, A = A, X = X, U = U, Z = Z, W = W)
dat$R = rbinom(N, 1, plogis(as.numeric(X %*% gammas) + Z*ZtoC + A*AtoC + 2))


### parametric estimation
result.if = continuous_survival_para_if(dat = dat, kappa = kappa, 
                        time.list = time.list, misspecify = NULL,
                        survival = "additive")
result.ipw = continuous_survival_para_ipw(dat = dat, kappa = kappa, 
                        time.list = time.list, misspecify = NULL)
result.plugin = continuous_survival_para_plugin(dat = dat, kappa = kappa, 
                        time.list = time.list, misspecify = NULL,
                        survival = "additive")

### nonparametric estimation
x.trt = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                   X4 = dat$X.4, X5 = dat$X.5, Z = dat$Z)
x.out = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                   X4 = dat$X.4, X5 = dat$X.5, Z = dat$Z, A = dat$A)
x.miss = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                    X4 = dat$X.4, X5 = dat$X.5, Z = dat$Z, A = dat$A)
x.instrument = data.frame(X1 = dat$X.1, X2 = dat$X.2, X3 = dat$X.3, 
                          X4 = dat$X.4, X5 = dat$X.5)

w.trt = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                   W4 = dat$W.4, W5 = dat$W.5, Z = dat$Z)
w.out = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                   W4 = dat$W.4, W5 = dat$W.5, Z = dat$Z, A = dat$A)
w.miss = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                    W4 = dat$W.4, W5 = dat$W.5, Z = dat$Z, A = dat$A)
w.instrument = data.frame(W1 = dat$W.1, W2 = dat$W.2, W3 = dat$W.3, 
                          W4 = dat$W.4, W5 = dat$W.5)

ranger.result.if = continuous_survival_if(dat = dat, kappa = kappa, x.trt = x.trt, x.out = x.out,
                            x.miss = x.miss, x.instrument = x.instrument, 
							w.trt = w.trt, w.out = w.out,
                            w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list,
                            nsplits = 2)

ranger.result.ipw = continuous_survival_ipw(dat = dat, kappa = kappa, x.trt = x.trt, x.out = x.out,
                            x.miss = x.miss, x.instrument = x.instrument, 
							w.trt = w.trt, w.out = w.out,
                            w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list,
                            nsplits = 2)


ranger.result.plugin = continuous_survival_plugin(dat = dat, kappa = kappa, x.trt = x.trt, x.out = x.out,
                            x.miss = x.miss, x.instrument = x.instrument, 
							w.trt = w.trt, w.out = w.out,
                            w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list,
                            nsplits = 2)

