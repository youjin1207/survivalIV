mu_fun = function(R = 1, Z, A, dat){ 
  fit = glm(Y ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + A + R, data = dat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                    colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                    diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  R = rep(1, nrow(dat)), A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}


w_fun = function(Z, A, dat){ 
  fit = glm(R ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + A, data = dat, family = binomial())
  R.Z.A1 = predict(fit, data.frame(Z = rep(Z, nrow(dat)), 
                                   age = dat$age, sex = dat$sex, 
                                   fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                   colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                   polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                   diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                   A = rep(1, nrow(dat))), type = "response") 
  R.Z.A0 = predict(fit, data.frame(Z = rep(Z, nrow(dat)), 
                                   age = dat$age, sex = dat$sex, 
                                   fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                   colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                   polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                   diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                   A = rep(0, nrow(dat))), type = "response")  
  return(A*R.Z.A1 + (1-A)*R.Z.A0)
}


pi_fun = function(Z, dat){
  fit = glm(A ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no, data = dat, family = binomial())
  A.Z= predict(fit, data.frame(Z = rep(Z, nrow(dat)), age = dat$age, sex = dat$sex, 
                               fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                               colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                               polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                               diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no), type = "response")
  return(A.Z)
}


delta_fun = function(Z, dat){
  fit = glm(Z ~ age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no, data = dat, family = binomial())
  Z.Z1 = predict(fit, data.frame(age = dat$age, sex = dat$sex, 
                                 fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                 colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                 polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                 diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no), type = "response")
  return(Z*Z.Z1 + (1-Z)*(1-Z.Z1))
}


pi_fun_naive = function(dat){
  fit = glm(A ~ age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no, data = dat, family = binomial())
  A.Z= predict(fit, data.frame(age = dat$age, sex = dat$sex, 
                               fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                               colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                               polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                               diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no), type = "response")
  return(A.Z)
}




binary_survival_para_if_colo = function(dat, time.list){

  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  subdat = dat
  
  # fit a treatment model
  pi_1 = pi_fun(Z = 1, subdat = subdat, dat = dat)
  pi_0 = rep(0, nrow(dat))
  
  # fit an instrumental variable model
  delta_1 = delta_fun(1, subdat = subdat, dat = dat)
  delta_0 = 1-delta_1
  
  # fit a censoring model
  omega_Z1_A1 = w_fun(Z = 1, A = 1, subdat = subdat, dat = dat)
  omega_Z1_A0 = w_fun(Z = 1, A = 0, subdat = subdat, dat = dat)
  omega_Z0_A0 = w_fun(Z = 0, A = 0, subdat = subdat, dat = dat)
  
  ## fit outcome models (time-varying)
  for(j in 1:k){
    subdat$Y = as.integer(time.list[j] < subdat$obs.Y)
    dat$Y = as.integer(time.list[j] < dat$obs.Y)
    
    mu_Z1_A1 = mu_fun(R = 1, Z = 1, A = 1, subdat = subdat, dat = dat)
    mu_Z1_A0 = mu_fun(R = 1, Z = 1, A = 0, subdat = subdat, dat = dat)
    mu_Z0_A0 = mu_fun(R = 1, Z = 0, A = 0, subdat = subdat, dat = dat)
    
    
    phi1_1 = mu_Z1_A1*pi_1  + mu_Z1_A0*(1-pi_1) + 
      dat$R*dat$Z*(dat$A*dat$Y - dat$A*mu_Z1_A1)/(delta_1*omega_Z1_A1) + 
      dat$Z*mu_Z1_A1*(dat$A - pi_1)/delta_1 +
      dat$R*dat$Z*((1-dat$A)*dat$Y - (1-dat$A)*mu_Z1_A0)/(delta_1*omega_Z1_A0) + 
      dat$Z*mu_Z1_A0*((1-dat$A) - (1-pi_1))/delta_1
    
    phi1_0 = 0 + mu_Z0_A0*(1-pi_0) + 
      0 +
      dat$R*(1-dat$Z)*((1-dat$A)*dat$Y - (1-dat$A)*mu_Z0_A0)/(delta_0*omega_Z0_A0) + 
      (1-dat$Z)*mu_Z0_A0*((1-dat$A) - (1-pi_0))/delta_0 
    
    
    phi2_1 = pi_1 + 
      dat$Z*(dat$A - pi_1)/delta_1
    phi2_0 = 0
    
    ## influence-function based estimator:
    ifvals[, j] = mean(phi1_1 - phi1_0, na.rm = TRUE)/mean(phi2_1 - phi2_0, na.rm = TRUE) 
    
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
    est.var[j] = var( (phi1_1 - phi1_0 - est.eff[j]*(phi2_1))/ mean(phi2_1) ) / length(phi1_1)   
  }
  
  
  return(list(est.eff = est.eff, est.var = est.var))
}


binary_survival_para_if_hazard = function(dat, time.list){
	
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  
  subdat = dat
  
  # fit a treatment model
  pi_1 = pi_fun(Z = 1, subdat = subdat, dat = dat)	
  pi_0 = rep(0, nrow(dat))
  
  # fit an instrumental variable model
  delta_1 = delta_fun(1, subdat = subdat, dat = dat)
  delta_0 = 1-delta_1
  
  # control the minimum/maximum of the nuisance =s
  delta_1 = ifelse(delta_1 > 0.95, 0.95, ifelse(delta_1 < 0.05, 0.05, delta_1))
  delta_0 = ifelse(delta_0 > 0.95, 0.95, ifelse(delta_0 < 0.05, 0.05, delta_0))
  pi_1 = ifelse(pi_1 > 0.95, 0.95, ifelse(pi_1 < 0.05, 0.05, pi_1))
  
  
  phi2_1 = pi_1 + 
    dat$Z*(dat$A - pi_1)/delta_1
  phi2_0 = rep(0, nrow(dat))
  
  
  ## at the very first time
  S_tau_Z1_A1 = S_tau_Z1_A0 = S_tau_Z0_A1 = S_tau_Z0_A0 = 
    G_tau_Z1_A1 = G_tau_Z1_A0 = G_tau_Z0_A1 = G_tau_Z0_A0 = matrix(1, length(time.list), nrow(dat))
  h_Z1_A1 =  h_Z1_A0 =  h_Z0_A1 = h_Z0_A0 = 
    omega_Z1_A1 = omega_Z1_A0 = omega_Z0_A1 = omega_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  D_Z1_A1 = D_Z1_A0 = D_Z0_A1 = D_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  
  
  coxfit = coxph(Surv(obs.Y, R) ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
                   colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
                   diabetes_f_yes + diabetes_f_no + A, data = dat)
  basefit = basehaz(coxfit, centered = FALSE)
  hazardfit = matrix(0, nrow(basefit), 2)
  hazardfit[,2] = basefit$time
  hazardfit[1,1] = basefit$hazard[1]*(basefit$time[1] - 0)
  for(t in 2:nrow(hazardfit)){
    hazardfit[t,1] = (basefit$hazard[t]-basefit$hazard[t-1])*(basefit$time[t] - basefit$time[t-1])
  }
  
  dat$R.censor = 1-dat$R 
  coxfit.censor = coxph(Surv(obs.Y, R.censor) ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
                          colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
                          diabetes_f_yes + diabetes_f_no + A, data = dat)
  basefit.censor = basehaz(coxfit.censor, centered = FALSE)
  hazardfit.censor = matrix(0, nrow(basefit.censor), 2)
  hazardfit.censor[,2] = basefit.censor$time
  hazardfit.censor[1,1] = basefit.censor$hazard[1]*(basefit.censor$time[1] - 0)
  for(t in 2:nrow(hazardfit.censor)){
    hazardfit.censor[t,1] = (basefit.censor$hazard[t]-basefit.censor$hazard[t-1])*(basefit.censor$time[t] - basefit.censor$time[t-1])
  }
  
  
  newdata_Z1_A1 = newdata_Z1_A0 = newdata_Z0_A0 = data.frame(Z = dat$Z, age = dat$age, sex = dat$sex,
                                                             fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                                             colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, polyps_f_yes = dat$polyps_f_yes,
                                                             polyps_f_no = dat$polyps_f_no,  diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no, A = dat$A)
  newdata_Z1_A1$Z = 1; newdata_Z1_A1$A = 1
  newdata_Z1_A0$Z = 1; newdata_Z1_A0$A = 0
  newdata_Z0_A0$Z = 0; newdata_Z0_A0$A = 0
  ## fit hazard models (time-varying)
  for(j in 2:k){
    subdat$Y = as.integer(subdat$obs.Y  == time.list[j] & subdat$R == 1)
    subdat$C = as.integer(subdat$obs.Y  == time.list[j] & subdat$R == 0)
    subdat$I = as.integer(time.list[j-1]  < subdat$obs.Y)
    
    dat$Y = as.integer(dat$obs.Y  == time.list[j] & dat$R == 1)
    dat$C = as.integer(dat$obs.Y  == time.list[j] & dat$R == 0)
    dat$I = as.integer(time.list[j-1] < dat$obs.Y)
    
    # fit a hazard function
    baselines.h = hazardfit[which.min(abs(basefit$time - time.list[j])),1] # baseline
    baselines.h1 = basefit$hazard[which.min(abs(basefit$time - time.list[j]))] # baseline
    baselines.h0 = basefit$hazard[which.min(abs(basefit$time - time.list[j-1]))] # baseline
    
    
    h_Z1_A1[j,] =  1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit)$coefficients[1:12,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit)$coefficients[1:12,1])))) 
    h_Z1_A0[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit)$coefficients[1:12,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit)$coefficients[1:12,1])))) 
    h_Z0_A0[j,] =	 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit)$coefficients[1:12,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit)$coefficients[1:12,1])))) 
    
    
    S_tau_Z1_A1[j,] = S_tau_Z1_A1[j-1,]*(1-h_Z1_A1[j,])
    S_tau_Z1_A0[j,] = S_tau_Z1_A0[j-1,]*(1-h_Z1_A0[j,])
    S_tau_Z0_A0[j,] = S_tau_Z0_A0[j-1,]*(1-h_Z0_A0[j,])
    
    
    S_tau_Z1_A1[j,] = ifelse(S_tau_Z1_A1[j,] < 0.05, 0.05, S_tau_Z1_A1[j,])
    S_tau_Z1_A0[j,] = ifelse(S_tau_Z1_A0[j,] < 0.05, 0.05, S_tau_Z1_A0[j,])
    #S_tau_Z0_A1[j,] = ifelse(S_tau_Z0_A1[j,] < 0.05, 0.05, S_tau_Z0_A1[j,])
    S_tau_Z0_A0[j,] = ifelse(S_tau_Z0_A0[j,] < 0.05, 0.05, S_tau_Z0_A0[j,])
    
    # fit a censoring hazard function
    baselines.h = hazardfit.censor[which.min(abs(basefit.censor$time - time.list[j])),1] # baseline
    baselines.h1 = basefit.censor$hazard[which.min(abs(basefit.censor$time - time.list[j]))] # baseline
    baselines.h0 = basefit.censor$hazard[which.min(abs(basefit.censor$time - time.list[j-1]))] # baseline
    omega_Z1_A1[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit.censor)$coefficients[1:12,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit.censor)$coefficients[1:12,1])))) 
    omega_Z1_A0[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit.censor)$coefficients[1:12,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit.censor)$coefficients[1:12,1])))) 
    omega_Z0_A0[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit.censor)$coefficients[1:12,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit.censor)$coefficients[1:12,1])))) 
    
    
    G_tau_Z1_A1[j,] = G_tau_Z1_A1[j-1,]*(1-omega_Z1_A1[j,])
    G_tau_Z1_A0[j,] = G_tau_Z1_A0[j-1,]*(1-omega_Z1_A0[j,])
    G_tau_Z0_A0[j,] = G_tau_Z0_A0[j-1,]*(1-omega_Z0_A0[j,])
    
    
    for(t in 2:j){
      D_Z1_A1[j,] = D_Z1_A1[j,] +
        - pi_1*(dat$Z*dat$A*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z1_A1[j,]*(as.integer(dat$obs.Y  == time.list[t] &
                                                                                                 dat$R == 1) - h_Z1_A1[t,])/(S_tau_Z1_A1[t, ]*G_tau_Z1_A1[t-1,]*pi_1*delta_1))
      D_Z1_A0[j,] = D_Z1_A0[j,] +
        - (1-pi_1)*(dat$Z*(1-dat$A)*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z1_A0[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                                                                         dat$R == 1) - h_Z1_A0[t,])/(S_tau_Z1_A0[t, ]*G_tau_Z1_A0[t-1,]*(1-pi_1)*delta_1))
      
      D_Z0_A0[j,] = D_Z0_A0[j,] +
        - (1-pi_0)*((1-dat$Z)*(1-dat$A)*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z0_A0[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                                                                             dat$R == 1) - h_Z0_A0[t,])/(S_tau_Z0_A0[t, ]*G_tau_Z0_A0[t-1,]*(1-pi_0)*delta_0))  
    }
    
    
    D_Z1_A1[j,] = D_Z1_A1[j,] + S_tau_Z1_A1[j,]*dat$Z*(dat$A - pi_1)/delta_1 + S_tau_Z1_A1[j,]*pi_1
    D_Z1_A0[j,] = D_Z1_A0[j,] + S_tau_Z1_A0[j,]*dat$Z*((1-dat$A) - (1-pi_1))/delta_1 + S_tau_Z1_A0[j,]*(1-pi_1)
    D_Z0_A1[j,] = rep(0, nrow(dat))
    D_Z0_A0[j,] = D_Z0_A0[j,] + S_tau_Z0_A0[j,]*(1-dat$Z)*((1-dat$A) -(1-pi_0))/delta_0 + S_tau_Z0_A0[j,]*(1-pi_0)
    
    
    ## influence-function based estimator:
    ifvals[, j] = mean( (D_Z1_A1[j,] + D_Z1_A0[j, ]) - (D_Z0_A1[j,] + D_Z0_A0[j,]), na.rm = TRUE) / 
      mean(phi2_1, na.rm = TRUE) 
  }
  
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
    est.var[j] = var( ((D_Z1_A1[j,] + D_Z1_A0[j, ]) - (D_Z0_A1[j,] + D_Z0_A0[j,]) - est.eff[j]*(phi2_1))/ mean(phi2_1 ) ) / length(phi2_1)
  }
  
  
  return(list(est.eff = est.eff, est.var = est.var))
}


binary_survival_para_naive = function(dat, time.list){
	
  n = nrow(dat)
  k = length(time.list) 
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  
  subdat = dat
  # fit treatment model
  pi_1 = pi_fun_naive(subdat = subdat, dat = dat)
  pi_1 = ifelse(pi_1 > 0.95, 0.95, ifelse(pi_1 < 0.05, 0.05, pi_1))
  
  
  S_tau_A1 = S_tau_A0 = G_tau_A1 = G_tau_A0 = matrix(0, length(time.list), nrow(dat))
  h_A1 = h_A0 = omega_A1 = omega_A0 = matrix(0, length(time.list), nrow(dat))
  D_A1 = D_A0 = matrix(0, length(time.list), nrow(dat))
  ## when j = 1
  S_tau_A1[1,] = S_tau_A0[1,] = rep(1, nrow(dat))
  G_tau_A1[1,] = G_tau_A0[1,] = rep(1, nrow(dat))
  D_A1[1,] = D_A0[1,] = rep(1, nrow(dat))
  ifvals[, 1] = mean(D_A1[1,] - D_A0[1,], na.rm = TRUE) 
  
  coxfit = coxph(Surv(obs.Y, R) ~ A + age + sex + fh_cancer_yes + fh_cancer_no +
                   colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
                   diabetes_f_yes + diabetes_f_no, data = dat)
  basefit = basehaz(coxfit, centered = FALSE)
  hazardfit = matrix(0, nrow(basefit), 2)
  hazardfit[,2] = basefit$time
  hazardfit[1,1] = basefit$hazard[1]*(basefit$time[1] - 0)
  for(t in 2:nrow(hazardfit)){
    hazardfit[t,1] = (basefit$hazard[t]-basefit$hazard[t-1])*(basefit$time[t] - basefit$time[t-1])
  }
  
  dat$R.censor = 1-dat$R 
  coxfit.censor = coxph(Surv(obs.Y, R.censor) ~ A + age + sex + fh_cancer_yes + fh_cancer_no +
                          colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
                          diabetes_f_yes + diabetes_f_no, data = dat)
  basefit.censor = basehaz(coxfit.censor, centered = FALSE)
  hazardfit.censor = matrix(0, nrow(basefit.censor), 2)
  hazardfit.censor[,2] = basefit.censor$time
  hazardfit.censor[1,1] = basefit.censor$hazard[1]*(basefit.censor$time[1] - 0)
  for(t in 2:nrow(hazardfit.censor)){
    hazardfit.censor[t,1] = (basefit.censor$hazard[t]-basefit.censor$hazard[t-1])*(basefit.censor$time[t] - basefit.censor$time[t-1])
  }
  
  
  newdata_A1 = newdata_A0 = data.frame(A = dat$A, age = dat$age, sex = dat$sex,
                                       fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                       colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, polyps_f_yes = dat$polyps_f_yes,
                                       polyps_f_no = dat$polyps_f_no,  diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no)
  newdata_A1$A = 1; newdata_A0$A = 0
  
  ## fit hazard models (time-varying)
  for(j in 2:k){
    
    dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
    dat$C = as.integer(dat$obs.Y == time.list[j]  & dat$R == 0)
    dat$I = as.integer(time.list[j-1]  < dat$obs.Y)
    
    
    baselines.h1 = basefit$hazard[which.min(abs(basefit$time - time.list[j]))] # baseline
    baselines.h0 = basefit$hazard[which.min(abs(basefit$time - time.list[j-1]))] # baseline
    
    h_A1[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(coxfit)$coefficients[1:11,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(coxfit)$coefficients[1:11,1]))))             
    h_A0[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(coxfit)$coefficients[1:11,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(coxfit)$coefficients[1:11,1])))) 
    
    
    
    S_tau_A1[j,] = S_tau_A1[j-1,]*(1-h_A1[j,])
    S_tau_A0[j,] = S_tau_A0[j-1,]*(1-h_A0[j,])
    S_tau_A1[j,] = ifelse(S_tau_A1[j,] < 0.05, 0.05, S_tau_A1[j,])
    S_tau_A0[j,] = ifelse(S_tau_A0[j,] < 0.05, 0.05, S_tau_A0[j,])
    
    baselines.h1 = basefit.censor$hazard[which.min(abs(basefit$time - time.list[j]))] # baseline
    baselines.h0 = basefit.censor$hazard[which.min(abs(basefit$time - time.list[j-1]))] # baseline
    omega_A1[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(coxfit.censor)$coefficients[1:11,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(coxfit.censor)$coefficients[1:11,1])))) 
    
    omega_A0[j,] =  1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(coxfit.censor)$coefficients[1:11,1])))) /
      exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(coxfit.censor)$coefficients[1:11,1])))) 
    
    
    G_tau_A1[j,] = G_tau_A1[j-1,]*(1-omega_A1[j,])
    G_tau_A0[j,] = G_tau_A0[j-1,]*(1-omega_A0[j,])
    G_tau_A1[j,] = ifelse(G_tau_A1[j,] < 0.05, 0.05, G_tau_A1[j,])
    G_tau_A0[j,] = ifelse(G_tau_A0[j,] < 0.05, 0.05, G_tau_A0[j,])
    
    D_A1[j, ] = D_A0[j,] = rep(0, nrow(dat))
    for(t in 2:j){
      D_A1[j,] = D_A1[j,] +
        -(dat$A*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_A1[j,]*(as.integer(time.list[t] == dat$obs.Y &
                                                                                  dat$R == 1) - h_A1[t,])/(S_tau_A1[t, ]*G_tau_A1[t-1,]*pi_1))
      D_A0[j,] = D_A0[j,] +
        -((1-dat$A)*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_A0[j,]*(as.integer(time.list[t] == dat$obs.Y &
                                                                                      dat$R == 1) - h_A0[t,])/(S_tau_A0[t, ]*G_tau_A0[t-1,]*(1-pi_1)))
    }
    D_A1[j,] = D_A1[j,]+ S_tau_A1[j,]
    D_A0[j,] = D_A0[j,]+ S_tau_A0[j,]
    
    ifvals[, j] = mean(D_A1[j,] - D_A0[j,], na.rm = TRUE) 
  }
  
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
    est.var[j] = var(D_A1[j,] - D_A0[j,]) / nrow(dat)
  }
  
  return(list(est.eff = est.eff, est.var = est.var))
}


