library(MASS)
library(ranger)
library(survival)
library(ahaz)

source("Code/core/nuisance_binary_survival_para.R")
source("Code/core/binary_survival_para_sim.R")
source("Code/core/binary_survival_ranger_sim.R")


### generate the simulation data
N = 1000
X = mvrnorm(N, rep(0,5), Sigma = diag(5))
deltas = c(-0.5, -0.5, 1, -1, -0.7)
alphas = c(0.1, 0.1, -0.2, -0.2, 0.3)*5
Z = rbinom(N, 1, plogis(as.numeric(X %*% deltas))) 
U = rnorm(N, 0, 1)
A = rbinom(N, 1, plogis(as.numeric(X %*% alphas) + 2*Z + 0.3*U + 0.1)) ## U : unmeasured confounder

betas = c(0.5, -0.5, 0.5, -0.5, 0.5)/50
AtoT = 1.0/50
UtoT = 0.5/50
gammas = c(0.3, 0.3, -0.3, -0.3, -0.3)/50
ZtoC = 0.5
AtoC = -0.5/50
t = c(1:180); t2 = t^2
tau = 180
coeff.Lambda = c(0.00, 0.0005, 0.0003)
coeff.Lambda2 = c(0.00, 0.0003, 0.0005)
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


inverse.Lambda2 = function(input, A){
  (1/(coeff.Lambda2[3]))*(-(coeff.Lambda2[2] + as.numeric(X %*% gammas) + A*AtoC)/2 +  
                            sqrt(coeff.Lambda2[3]*input + (coeff.Lambda2[2] + as.numeric(X %*% gammas) + A*AtoC)^2/4 -coeff.Lambda2[1]*coeff.Lambda2[3]))
}



U.cf = runif(N, 0, 1)
obs.T1 = inverse.Lambda((-log(U.cf)), 1)
obs.T0 = inverse.Lambda((-log(U.cf)), 0)
obs.T = inverse.Lambda((-log(U.cf)), A)

U.censor = runif(N, 0, 1)
obs.C = inverse.Lambda2(-log(U.censor), A)
#obs.C = runif(N, 50, 100)
obs.Y = pmin(obs.T, obs.C)


## transform the first four baseline covariates
W = X
W[,1] = exp(X[,1]/2)
W[,2] = X[,2] / (1+exp(X[,1])) + 10
W[,3] = (X[,1]*X[,3]/25 + 0.6 )^3
W[,4] = (X[,2] + X[,4] + 20)^2


time.list = seq(0, 30, 1)
R = as.integer(obs.T < obs.C)
obs.Y = pmin(obs.T, obs.C)
obs.T = as.integer(obs.T)
obs.C = as.integer(obs.C)
obs.Y = as.integer(obs.Y)
obs.T1= as.integer(obs.T1)
obs.T0 = as.integer(obs.T0)

#dat$R = rbinom(N, 1, plogis(as.numeric(X %*% gammas) + Z*ZtoC + A*AtoC))

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

dat = data.frame(obs.Y = obs.Y, obs.T = obs.T, obs.C = obs.C, R = R, R.censor = 1-R, A = A, X = X, U = U, Z = Z, W = W)


### parametric estimation
result.if = binary_survival_para_if(dat = dat, time.list = time.list, misspecify = NULL, survival = "additive", censor = "additive")

result.ipw = binary_survival_para_ipw(dat = dat, time.list = time.list, misspecify = NULL, censor = "additive")

result.plugin = binary_survival_para_plugin(dat = dat, time.list = time.list, misspecify = NULL, survival = "additive")


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

ranger.result.if = binary_survival_ranger_if_hazard(dat = dat, x.trt = x.trt, x.out = x.out,
                            x.miss = x.miss, x.instrument = x.instrument, 
              w.trt = w.trt, w.out = w.out,
                            w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list,
                            nsplits = 2, misspecify = NULL)

ranger.result.ipw = binary_survival_ranger_hazard_ipw(dat = dat, x.trt = x.trt, x.out = x.out,
                            x.miss = x.miss, x.instrument = x.instrument, 
              w.trt = w.trt, w.out = w.out,
                            w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list,
                            nsplits = 2, misspecify = NULL)

ranger.result.plugin = binary_survival_ranger_hazard_plugin(dat = dat, x.trt = x.trt, x.out = x.out,
                            x.miss = x.miss, x.instrument = x.instrument, 
              w.trt = w.trt, w.out = w.out,
                            w.miss = w.miss, w.instrument = w.instrument, 
                            time.list = time.list,
                            nsplits = 2, misspecify = NULL)


