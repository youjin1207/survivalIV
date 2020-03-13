## X : dimension of five
mu_fun = function(R = 1, Z, A, subdat, dat){ 
  fit = glm(Y ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + R + A, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), X.1 = dat$X.1, 
                                  X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                                  R = rep(1, nrow(dat)), A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}


w_fun = function(Z, A, subdat, dat){ 
  fit = glm(R ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + A, data = subdat, family = binomial())
  R.Z.A1 = predict(fit, data.frame(Z = rep(Z, nrow(dat)), 
                                   X.1 = dat$X.1, 
                                   X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                                   A = rep(1, nrow(dat))), type = "response") 
  R.Z.A0 = predict(fit, data.frame(Z = rep(Z, nrow(dat)), 
                                   X.1 = dat$X.1, 
                                   X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5, 
                                   A = rep(0, nrow(dat))), type = "response")  
  return(A*R.Z.A1 + (1-A)*R.Z.A0)
}


pi_fun = function(Z, subdat, dat){
  fit = glm(A ~ Z + X.1 + X.2 + X.3 + X.4 + X.5, data = subdat, family = binomial())
  A.Z= predict(fit, data.frame(Z = rep(Z, nrow(dat)), X.1 = dat$X.1, 
                               X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, 
                               X.5 = dat$X.5), type = "response")
  return(A.Z)
}


delta_fun = function(Z, subdat, dat){
  fit = glm(Z ~ X.1 + X.2 + X.3 + X.4 + X.5, data = subdat, family = binomial())
  Z.Z1 = predict(fit, data.frame(X.1 = dat$X.1, 
                                 X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, 
                                 X.5 = dat$X.5), type = "response")
  return(Z*Z.Z1 + (1-Z)*(1-Z.Z1))
}


hazard_funA = function(I = 1, Z, A, subdat, dat){ 
  fit = glm(Y ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + I + A, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), X.1 = dat$X.1, 
                                  X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                                  I = rep(1, nrow(dat)),
                                  A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}


censor_funA = function(I = 1, Z, A, subdat, dat){ 
  fit = glm(C ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + I + A, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), X.1 = dat$X.1, 
                                  X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                                  I = rep(1, nrow(dat)),
                                  A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}

pi_fun_naive = function(subdat, dat){
  fit = glm(A ~ X.1 + X.2 + X.3 + X.4 + X.5, data = subdat, family = binomial())
  A.Z = predict(fit, data.frame(X.1 = dat$X.1, 
                                  X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5), type = "response")
  return(A.Z)
}


hazard_fun_naive = function(I = 1, A, subdat, dat){ 
  fit = glm(Y ~ A + X.1 + X.2 + X.3 + X.4 + X.5 + I, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(A = rep(A, nrow(dat)), X.1 = dat$X.1, 
                                  X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                                  I = rep(1, nrow(dat))), type = "response")
  return(Y.fit)
}


censor_fun_naive = function(I = 1, A, subdat, dat){ 
  fit = glm(C ~ A + X.1 + X.2 + X.3 + X.4 + X.5 + I, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(A = rep(A, nrow(dat)), X.1 = dat$X.1, 
                                  X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                                  I = rep(1, nrow(dat))), type = "response")
  return(Y.fit)
}


##### misspecified models #####
mu_fun_wrong = function(R = 1, Z, A, subdat, dat){ 
  fit = glm(Y ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + R + A, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), W.1 = dat$W.1, 
                                  W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                                  R = rep(1, nrow(dat)), A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}


w_fun_wrong = function(Z, A, subdat, dat){ 
  fit = glm(R ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + A, data = subdat, family = binomial())
  R.Z.A1 = predict(fit, data.frame(Z = rep(Z, nrow(dat)), 
                                   W.1 = dat$W.1, 
                                   W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                                   A = rep(1, nrow(dat))), type = "response") 
  R.Z.A0 = predict(fit, data.frame(Z = rep(Z, nrow(dat)), 
                                   W.1 = dat$W.1, 
                                   W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5, 
                                   A = rep(0, nrow(dat))), type = "response")  
  return(A*R.Z.A1 + (1-A)*R.Z.A0)
}


pi_fun_wrong = function(Z, subdat, dat){
  fit = glm(A ~ Z + W.1 + W.2 + W.3 + W.4 + W.5, data = subdat, family = binomial())
  A.Z= predict(fit, data.frame(Z = rep(Z, nrow(dat)), W.1 = dat$W.1, 
                               W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, 
                               W.5 = dat$W.5), type = "response")
  return(A.Z)
}


delta_fun_wrong = function(Z, subdat, dat){
  fit = glm(Z ~ W.1 + W.2 + W.3 + W.4 + W.5, data = subdat, family = binomial())
  Z.Z1 = predict(fit, data.frame(W.1 = dat$W.1, 
                                 W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, 
                                 W.5 = dat$W.5), type = "response")
  return(Z*Z.Z1 + (1-Z)*(1-Z.Z1))
}


hazard_fun_wrongA = function(I = 1, Z, A, subdat, dat){ 
  fit = glm(Y ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + I + A, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), W.1 = dat$W.1, 
                                  W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                                  I = rep(1, nrow(dat)), 
                                  A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}


censor_fun_wrongA = function(I = 1, Z, A, subdat, dat){ 
  fit = glm(C ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + I + A, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), W.1 = dat$W.1, 
                                  W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                                  I = rep(1, nrow(dat)),
                                  A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}

pi_fun_naive_wrong = function(subdat, dat){
  fit = glm(A ~ W.1 + W.2 + W.3 + W.4 + W.5, data = subdat, family = binomial())
  A.Z = predict(fit, data.frame(W.1 = dat$W.1, 
                                  W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5), type = "response")
  return(A.Z)
}


hazard_fun_naive_wrong = function(I = 1, A, subdat, dat){ 
  fit = glm(Y ~ A + W.1 + W.2 + W.3 + W.4 + W.5 + I, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(A = rep(A, nrow(dat)),  W.1 = dat$W.1, 
                                  W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                                  I = rep(1, nrow(dat))), type = "response")
  return(Y.fit)
}


censor_fun_naive_wrong = function(I = 1, A, subdat, dat){ 
  fit = glm(C ~ A + W.1 + W.2 + W.3 + W.4 + W.5 + I, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(A = rep(A, nrow(dat)),  W.1 = dat$W.1, 
                                  W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                                  I = rep(1, nrow(dat))), type = "response")
  return(Y.fit)
}




