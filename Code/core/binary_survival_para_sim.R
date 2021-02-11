source("survIV/nuisance_binary_survival_para.R")
binary_survival_para_if = function(dat, time.list, misspecify = NULL, survival = "Cox"){
  ## dat: dataframe containing 
  ## time.list: time point to measure the survival difference
  ## misspecify: c("treatment", "instrument", "censoring", "outcome"): a vector of characters that specify the misspecified nuisance functions
  ## survival : c("Cox", "additive"): a choice of parametric models for survival outcomes

  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  subdat = dat

  # fit treatment model
  if("treatment" %in% misspecify){
    pi_1 = pi_fun_wrong(Z = 1, subdat = subdat, dat = dat)
    pi_0 = pi_fun_wrong(Z = 0, subdat = subdat, dat = dat)
  }else{
    pi_1 = pi_fun(Z = 1, subdat = subdat, dat = dat)
    pi_0 = pi_fun(Z = 0, subdat = subdat, dat = dat)
  }

  if("instrument" %in% misspecify){
    delta_1 = delta_fun_wrong(1, subdat = subdat, dat = dat)
    delta_0 = 1-delta_1
  }else{
    delta_1 = delta_fun(1, subdat = subdat, dat = dat)
    delta_0 = 1-delta_1
  }

  if("censoring" %in% misspecify){
      omega_Z1_A1 = w_fun_wrong(Z = 1, A = 1, subdat = subdat, dat = dat)
      omega_Z1_A0 = w_fun_wrong(Z = 1, A = 0, subdat = subdat, dat = dat)
      omega_Z0_A1 = w_fun_wrong(Z = 0, A = 1, subdat = subdat, dat = dat)
      omega_Z0_A0 = w_fun_wrong(Z = 0, A = 0, subdat = subdat, dat = dat) 
  }else{
      omega_Z1_A1 = w_fun(Z = 1, A = 1, subdat = subdat, dat = dat)
      omega_Z1_A0 = w_fun(Z = 1, A = 0, subdat = subdat, dat = dat)
      omega_Z0_A1 = w_fun(Z = 0, A = 1, subdat = subdat, dat = dat)
      omega_Z0_A0 = w_fun(Z = 0, A = 0, subdat = subdat, dat = dat)
  } 


  
  if("outcome" %in% misspecify){
      newdata_Z1_A1 = data.frame(Z = rep(1, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z1_A0 = data.frame(Z = rep(1, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z0_A1 = data.frame(Z = rep(0, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z0_A0 = data.frame(Z = rep(0, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y) ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + A + R, data = dat)
        basefit = basehaz(coxfit, centered = FALSE)
    
      }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$Z, dat$W.1, dat$W.2, dat$W.3, dat$W.4, dat$W.5, dat$A, dat$R)
        colnames(design) = c("Z", "W.1", "W.2", "W.3", "W.4", "W.5", "A", "R")
        addfit = ahaz(Surv(dat$obs.YY), design)
        addbase = predict(addfit, type = "cumhaz")
      } 
  }else{
    newdata_Z1_A1 = data.frame(Z = rep(1, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
    newdata_Z1_A0 = data.frame(Z = rep(1, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
    newdata_Z0_A1 = data.frame(Z = rep(0, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
    newdata_Z0_A0 = data.frame(Z = rep(0, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))

    if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y) ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + A + R, data = dat)
        basefit = basehaz(coxfit, centered = FALSE)
    
    }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$Z, dat$X.1, dat$X.2, dat$X.3, dat$X.4, dat$X.5, dat$A, dat$R)
        colnames(design) = c("Z", "X.1", "X.2", "X.3", "X.4", "X.5", "A", "R")
        addfit = ahaz(Surv(dat$obs.YY), design)
        addbase = predict(addfit, type = "cumhaz")

    } 
  }  
  for(j in 1:k){
    subdat$Y = as.integer(time.list[j] < subdat$obs.Y)
    dat$Y = as.integer(time.list[j] < dat$obs.Y)
  
  
    if(survival == "Cox"){
        baseline = basefit[which.min(abs(basefit$time - time.list[j])),1] # baseline
        mu_Z1_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z1_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z0_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z0_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))

    }else if(survival == "additive"){
       addbaseline = addbase$cumhaz[which.min(abs(addbase$times - time.list[j]))] # baseline
       mu_Z1_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       mu_Z1_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       mu_Z0_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       mu_Z0_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
    } 
      
       
    phi1_1 = mu_Z1_A1*pi_1  + mu_Z1_A0*(1-pi_1) + 
         dat$R*dat$Z*(dat$A*dat$Y - dat$A*mu_Z1_A1)/(delta_1*omega_Z1_A1) + 
         dat$Z*mu_Z1_A1*(dat$A - pi_1)/delta_1 +
         dat$R*dat$Z*((1-dat$A)*dat$Y - (1-dat$A)*mu_Z1_A0)/(delta_1*omega_Z1_A0) + 
         dat$Z*mu_Z1_A0*((1-dat$A) - (1-pi_1))/delta_1



    phi1_0 = mu_Z0_A1*pi_0  + mu_Z0_A0*(1-pi_0) + 
         dat$R*(1-dat$Z)*(dat$A*dat$Y - dat$A*mu_Z0_A1)/(delta_0*omega_Z0_A1) + 
         (1-dat$Z)*mu_Z0_A1*(dat$A - pi_0)/delta_0 +
         dat$R*(1-dat$Z)*((1-dat$A)*dat$Y - (1-dat$A)*mu_Z0_A0)/(delta_0*omega_Z0_A0) + 
         (1-dat$Z)*mu_Z0_A0*((1-dat$A) - (1-pi_0))/delta_0 


      phi2_1 = pi_1 + 
          dat$Z*(dat$A - pi_1)/delta_1
      phi2_0 = pi_0  + 
          (1-dat$Z)*(dat$A - pi_0)/delta_0

      ## influence-function based estimator:
      ifvals[, j] = mean(phi1_1 - phi1_0, na.rm = TRUE) / 
                  mean(phi2_1 - phi2_0, na.rm = TRUE) 
      est.eff[j] = mean(ifvals[,j]) # for each delta vaule            
      est.var[j] = var( (phi1_1 - phi1_0 - est.eff[j]*(phi2_1-phi2_0))/ mean(phi2_1 - phi2_0) ) / length(phi1_1)  
    }

  
    return(list(est.eff = est.eff, est.var = est.var))
}


## add ipw and plug-in estimator
binary_survival_para_ipw = function(dat, time.list, misspecify = NULL){

  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = rep(NA, k)
  subdat = dat

  # fit treatment model
  if("treatment" %in% misspecify){
    pi_1 = pi_fun_wrong(Z = 1, subdat = subdat, dat = dat)
    pi_0 = pi_fun_wrong(Z = 0, subdat = subdat, dat = dat)
  }else{
    pi_1 = pi_fun(Z = 1, subdat = subdat, dat = dat)
    pi_0 = pi_fun(Z = 0, subdat = subdat, dat = dat)
  }

  if("instrument" %in% misspecify){
    delta_1 = delta_fun_wrong(1, subdat = subdat, dat = dat)
    delta_0 = 1-delta_1
  }else{
    delta_1 = delta_fun(1, subdat = subdat, dat = dat)
    delta_0 = 1-delta_1
  }

  if("censoring" %in% misspecify){
      omega_Z1_A1 = w_fun_wrong(Z = 1, A = 1, subdat = subdat, dat = dat)
      omega_Z1_A0 = w_fun_wrong(Z = 1, A = 0, subdat = subdat, dat = dat)
      omega_Z0_A1 = w_fun_wrong(Z = 0, A = 1, subdat = subdat, dat = dat)
      omega_Z0_A0 = w_fun_wrong(Z = 0, A = 0, subdat = subdat, dat = dat) 
  }else{
      omega_Z1_A1 = w_fun(Z = 1, A = 1, subdat = subdat, dat = dat)
      omega_Z1_A0 = w_fun(Z = 1, A = 0, subdat = subdat, dat = dat)
      omega_Z0_A1 = w_fun(Z = 0, A = 1, subdat = subdat, dat = dat)
      omega_Z0_A0 = w_fun(Z = 0, A = 0, subdat = subdat, dat = dat)
  } 

  ## fit outcome models (time-varying)
  for(j in 1:k){
      subdat$Y = as.integer(time.list[j] < subdat$obs.Y)
      dat$Y = as.integer(time.list[j] < dat$obs.Y)
            
           
      phi1_1 = mean((dat$Y*dat$R*dat$Z*dat$A/(omega_Z1_A1*pi_1*delta_1)))*
                mean((dat$A*dat$Z / (delta_1))) + 
                mean((dat$Y*dat$R*dat$Z*(1-dat$A)/(omega_Z1_A0*(1-pi_1)*delta_1)))*
                mean(((1-dat$A)*dat$Z / (delta_1)))
        
      phi1_0 =  mean((dat$Y*dat$R*(1-dat$Z)*dat$A/(omega_Z0_A1*pi_0*delta_0)))*
                mean((dat$A*(1-dat$Z)/ (delta_0)) ) + 
                mean((dat$Y*dat$R*(1-dat$Z)*(1-dat$A)/(omega_Z0_A0*(1-pi_0)*delta_0)))*
                mean(((1-dat$A)*(1-dat$Z) / (delta_0)))
        
      
      phi2_1 = dat$Z*dat$A/delta_1
      phi2_0 = (1-dat$Z)*dat$A/delta_0
       
      ## influence-function based estimator:
      ifvals[, j] = mean(phi1_1 - phi1_0) / 
                      mean(phi2_1- phi2_0) 
    }
    

    # compute estimator (using)
    for(j in 1:k){
      est.eff[j] = mean(ifvals[,j]) # for each delta vaule
    }


    return(list(est = est.eff))
}



binary_survival_para_plugin = function(dat, time.list, misspecify = NULL, survival = "Cox"){

  ## dat: data.frame(obs.T = obs.T, A = A, X = X, U = U, Z = Z)
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = rep(NA, k)
  subdat = dat

  # fit treatment model
  if("treatment" %in% misspecify){
    pi_1 = pi_fun_wrong(Z = 1, subdat = subdat, dat = dat)
    pi_0 = pi_fun_wrong(Z = 0, subdat = subdat, dat = dat)
  }else{
    pi_1 = pi_fun(Z = 1, subdat = subdat, dat = dat)
    pi_0 = pi_fun(Z = 0, subdat = subdat, dat = dat)
  }

   
  if("outcome" %in% misspecify){
      newdata_Z1_A1 = data.frame(Z = rep(1, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z1_A0 = data.frame(Z = rep(1, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z0_A1 = data.frame(Z = rep(0, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z0_A0 = data.frame(Z = rep(0, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y) ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + A + R, data = dat)
        basefit = basehaz(coxfit, centered = FALSE)
    
      }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$Z, dat$W.1, dat$W.2, dat$W.3, dat$W.4, dat$W.5, dat$A, dat$R)
        colnames(design) = c("Z", "W.1", "W.2", "W.3", "W.4", "W.5", "A", "R")
        addfit = ahaz(Surv(dat$obs.YY), design)
        addbase = predict(addfit, type = "cumhaz")
      } 
  }else{
    newdata_Z1_A1 = data.frame(Z = rep(1, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
    newdata_Z1_A0 = data.frame(Z = rep(1, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
    newdata_Z0_A1 = data.frame(Z = rep(0, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
    newdata_Z0_A0 = data.frame(Z = rep(0, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))

    if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y) ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + A + R, data = dat)
        basefit = basehaz(coxfit, centered = FALSE)
    
    }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$Z, dat$X.1, dat$X.2, dat$X.3, dat$X.4, dat$X.5, dat$A, dat$R)
        colnames(design) = c("Z", "X.1", "X.2", "X.3", "X.4", "X.5", "A", "R")
        addfit = ahaz(Surv(dat$obs.YY), design)
        addbase = predict(addfit, type = "cumhaz")

    } 
  }  


  ## fit outcome models (time-varying)
  for(j in 1:k){
      subdat$Y = as.integer(time.list[j] < subdat$obs.Y)
      dat$Y = as.integer(time.list[j] < dat$obs.Y)
        
     if(survival == "Cox"){
        baseline = basefit[which.min(abs(basefit$time - time.list[j])),1] # baseline
        mu_Z1_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z1_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z0_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z0_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))

    }else if(survival == "additive"){
       addbaseline = addbase$cumhaz[which.min(abs(addbase$times - time.list[j]))] # baseline
       mu_Z1_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       mu_Z1_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       mu_Z0_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       mu_Z0_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
    } 
      

       
      phi1_1 = mu_Z1_A1*pi_1  + mu_Z1_A0*(1-pi_1) 
  
      phi1_0 = mu_Z0_A1*pi_0  + mu_Z0_A0*(1-pi_0) 
    
      phi2_1 = pi_1 
      phi2_0 = pi_0  

      ifvals[, j] = mean(phi1_1 - phi1_0) / 
          mean(phi2_1 - phi2_0) 
  }
    

  for(j in 1:k){
      est.eff[j] = mean(ifvals[,j]) # for each delta vaule
  }

  return(list(est = est.eff))
}



binary_survival_para_if_hazard = function(dat, time.list,  misspecify = NULL, survival = "Cox"){
  
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  subdat = dat

  # fit a treatment model
  if("treatment" %in% misspecify){
      pi_1 = pi_fun_wrong(Z = 1, subdat = subdat, dat = dat)
      pi_0 = pi_fun_wrong(Z = 0, subdat = subdat, dat = dat)
  }else{
      pi_1 = pi_fun(Z = 1, subdat = subdat, dat = dat)
      pi_0 = pi_fun(Z = 0, subdat = subdat, dat = dat)
  }
  
  # fit an instrument model
  if("instrument" %in% misspecify){
    delta_1 = delta_fun_wrong(1, subdat = subdat, dat = dat)
    delta_0 = 1-delta_1
  }else{
    delta_1 = delta_fun(1, subdat = subdat, dat = dat)
    delta_0 = 1-delta_1
  }

  #delta_1 = ifelse(delta_1 > 0.95, 0.95, ifelse(delta_1 < 0.05, 0.05, delta_1))
  #delta_0 = ifelse(delta_0 > 0.95, 0.95, ifelse(delta_0 < 0.05, 0.05, delta_0))
  #pi_1 = ifelse(pi_1 > 0.95, 0.95, ifelse(pi_1 < 0.05, 0.05, pi_1))
  #pi_0 = ifelse(pi_0 > 0.95, 0.95, ifelse(pi_0 < 0.05, 0.05, pi_0))
  
  phi2_1 = pi_1 + 
          dat$Z*(dat$A - pi_1)/delta_1
  phi2_0 = pi_0  + 
          (1-dat$Z)*(dat$A - pi_0)/delta_0
  

  S_tau_Z1_A1 = S_tau_Z1_A0 = S_tau_Z0_A1 = S_tau_Z0_A0 = 
  G_tau_Z1_A1 = G_tau_Z1_A0 = G_tau_Z0_A1 = G_tau_Z0_A0 = matrix(1, length(time.list), nrow(dat))
  h_Z1_A1 =  h_Z1_A0 =  h_Z0_A1 = h_Z0_A0 = 
  omega_Z1_A1 = omega_Z1_A0 = omega_Z0_A1 = omega_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  D_Z1_A1 = D_Z1_A0 = D_Z0_A1 = D_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  ## when j = 1
  D_Z1_A1[1,] = D_Z1_A0[1,] = D_Z0_A1[1,] = D_Z0_A0[1,] = rep(1, nrow(dat))
  ifvals[, 1] = mean(D_Z1_A1[1,] + D_Z1_A0[1, ] - (D_Z0_A1[1,] + D_Z0_A0[1,] )) /
      mean(phi2_1 - phi2_0)

 if("hazard" %in% misspecify){
      newdata_Z1_A1 = data.frame(Z = rep(1, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)))
      newdata_Z1_A0 = data.frame(Z = rep(1, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)))
      newdata_Z0_A1 = data.frame(Z = rep(0, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)))
      newdata_Z0_A0 = data.frame(Z = rep(0, nrow(dat)), W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)))
      if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y, R) ~ Z + W.1 + W.2 + W.3 + W.4 + W.5 + A, data = dat)
        basefit = basehaz(coxfit, centered = FALSE)
        hazardfit = matrix(0, nrow(basefit), 2)
        hazardfit[,2] = basefit$time
        hazardfit[1,1] = basefit$hazard[1]*(basefit$time[1] - 0)
        for(t in 2:nrow(hazardfit)){
            hazardfit[t,1] = (basefit$hazard[t]-basefit$hazard[t-1])*(basefit$time[t] - basefit$time[t-1])
        }
    
      }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$Z, dat$W.1, dat$W.2, dat$W.3, dat$W.4, dat$W.5, dat$A)
        colnames(design) = c("Z", "W.1", "W.2", "W.3", "W.4", "W.5", "A")
        addfit = ahaz(Surv(dat$obs.YY, dat$R), design)
        addbase = predict(addfit, type = "cumhaz")
        addbasehazard = matrix(0, length(addbase$times), 2)
        addbasehazard [,2] = addbase$times
        addbasehazard [1,1] = addbase$cumhaz[1]*(addbase$times[1] - 0)
         for(t in 2:nrow(addbasehazard )){
           addbasehazard[t,1] = (addbase$cumhaz[t] - addbase$cumhaz[t-1])*(addbase$times[t] - addbase$times[t-1])
        }
      } 
  }else{
    newdata_Z1_A1 = data.frame(Z = rep(1, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)))
    newdata_Z1_A0 = data.frame(Z = rep(1, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)))
    newdata_Z0_A1 = data.frame(Z = rep(0, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)))
    newdata_Z0_A0 = data.frame(Z = rep(0, nrow(dat)), X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)))

    if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y, R) ~ Z + X.1 + X.2 + X.3 + X.4 + X.5 + A, data = dat)
         basefit = basehaz(coxfit, centered = FALSE)
         hazardfit = matrix(0, nrow(basefit), 2)
         hazardfit[,2] = basefit$time
         hazardfit[1,1] = basefit$hazard[1]*(basefit$time[1] - 0)
         for(t in 2:nrow(hazardfit)){
            hazardfit[t,1] = (basefit$hazard[t]-basefit$hazard[t-1])*(basefit$time[t] - basefit$time[t-1])
         }
    }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$Z, dat$X.1, dat$X.2, dat$X.3, dat$X.4, dat$X.5, dat$A)
        colnames(design) = c("Z", "X.1", "X.2", "X.3", "X.4", "X.5", "A")
        addfit = ahaz(Surv(dat$obs.YY, dat$R), design)
        addbase = predict(addfit, type = "cumhaz")
        
        tau = as.integer(max(dat$obs.Y)); tmp.cum = matrix(0, tau, 2)
        for(t in 1:tau){
         dummy = which(as.integer(addbase$times) == t)[1]
         tmp.cum[t,] = c(addbase$cumhaz[dummy], addbase$times[dummy])
        }
        tmp.cum = na.omit(tmp.cum)
        addbasehazard = matrix(0, nrow(tmp.cum), 2)
        addbasehazard[,2] = as.integer(tmp.cum[,2])
        addbasehazard[1,1] = tmp.cum[1,1]*(tmp.cum[1,2] - 0)
         for(t in 2:nrow(tmp.cum)){
           addbasehazard[t,1] = (tmp.cum[t,1]-tmp.cum[t-1,1])*(tmp.cum[t,2] - tmp.cum[t-1,2])
        }
    } 
  }  


  for(j in 2:k){

        subdat$Y = as.integer(subdat$obs.Y == time.list[j] & subdat$R == 1)
        subdat$C = as.integer(subdat$obs.Y == time.list[j] & subdat$R == 0)
        subdat$I = as.integer(time.list[j-1]  < subdat$obs.Y)
        

        dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
        dat$C = as.integer(dat$obs.Y == time.list[j] & dat$R == 0)
        dat$I = as.integer(time.list[j-1]  < dat$obs.Y)

        if(survival == "Cox"){
          baselines.h = hazardfit[which.min(abs(basefit$time - time.list[j])),1] # baseline

          baselines.h1 = basefit$hazard[which.min(abs(basefit$time - time.list[j]))] # baseline
          baselines.h0 = basefit$hazard[which.min(abs(basefit$time - time.list[j-1]))] # baseline
          h_Z1_A1[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) /
          exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) 
          h_Z1_A0[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) /
          exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) 
                       
          h_Z0_A1[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) /
          exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) 
          h_Z0_A0[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) /
          exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(coxfit)$coefficients[1:7,1])))) 

        }else if(survival == "additive"){
         baselines.h = addbasehazard[which.min(abs( addbasehazard[,2] - time.list[j])),1] # baseline
          h_Z1_A1[j,] = baselines.h + as.numeric(as.matrix(newdata_Z1_A1) %*% as.matrix(summary(addfit)$coefficients[1:7,1]))
          h_Z1_A0[j,] = baselines.h + as.numeric(as.matrix(newdata_Z1_A0) %*% as.matrix(summary(addfit)$coefficients[1:7,1]))
          h_Z0_A1[j,] = baselines.h + as.numeric(as.matrix(newdata_Z0_A1) %*% as.matrix(summary(addfit)$coefficients[1:7,1]))
          h_Z0_A0[j,] = baselines.h + as.numeric(as.matrix(newdata_Z0_A0) %*% as.matrix(summary(addfit)$coefficients[1:7,1]))
        } 
      
          
      S_tau_Z1_A1[j,] = S_tau_Z1_A1[j-1,]*(1-h_Z1_A1[j,])
      S_tau_Z1_A0[j,] = S_tau_Z1_A0[j-1,]*(1-h_Z1_A0[j,])
      S_tau_Z0_A1[j,] = S_tau_Z0_A1[j-1,]*(1-h_Z0_A1[j,])
      S_tau_Z0_A0[j,] = S_tau_Z0_A0[j-1,]*(1-h_Z0_A0[j,])
      
      #S_tau_Z1_A1[j,] = ifelse(S_tau_Z1_A1[j,] < 0.05, 0.05, S_tau_Z1_A1[j,])
      #S_tau_Z1_A0[j,] = ifelse(S_tau_Z1_A0[j,] < 0.05, 0.05, S_tau_Z1_A0[j,])
      #S_tau_Z0_A1[j,] = ifelse(S_tau_Z0_A1[j,] < 0.05, 0.05, S_tau_Z0_A1[j,])
      #S_tau_Z0_A0[j,] = ifelse(S_tau_Z0_A0[j,] < 0.05, 0.05, S_tau_Z0_A0[j,])

      if("censoring" %in% misspecify){
          omega_Z1_A1[j,] = censor_fun_wrongA(I = 1, Z = 1, A = 1, subdat = subdat, dat = dat)
          omega_Z1_A0[j,] = censor_fun_wrongA(I = 1, Z = 1, A = 0, subdat = subdat, dat = dat)
          omega_Z0_A1[j,] = censor_fun_wrongA(I = 1, Z = 0, A = 1, subdat = subdat, dat = dat)
          omega_Z0_A0[j,] = censor_fun_wrongA(I = 1, Z = 0, A = 0, subdat = subdat, dat = dat)
      }else{
          omega_Z1_A1[j,] = censor_funA(I = 1, Z = 1, A = 1, subdat = subdat, dat = dat)
          omega_Z1_A0[j,] = censor_funA(I = 1, Z = 1, A = 0, subdat = subdat, dat = dat)
          omega_Z0_A1[j,] = censor_funA(I = 1, Z = 0, A = 1, subdat = subdat, dat = dat)
          omega_Z0_A0[j,] = censor_funA(I = 1, Z = 0, A = 0, subdat = subdat, dat = dat)
      }

      G_tau_Z1_A1[j,] = G_tau_Z1_A1[j-1,]*(1-omega_Z1_A1[j,])
      G_tau_Z1_A0[j,] = G_tau_Z1_A0[j-1,]*(1-omega_Z1_A0[j,])
      G_tau_Z0_A1[j,] = G_tau_Z0_A1[j-1,]*(1-omega_Z0_A1[j,])
      G_tau_Z0_A0[j,] = G_tau_Z0_A0[j-1,]*(1-omega_Z0_A0[j,])

      #G_tau_Z1_A1[j,] = ifelse(G_tau_Z1_A1[j,] < 0.05, 0.05, G_tau_Z1_A1[j,])
      #G_tau_Z1_A0[j,] = ifelse(G_tau_Z1_A0[j,] < 0.05, 0.05, G_tau_Z1_A0[j,])
      #G_tau_Z0_A1[j,] = ifelse(G_tau_Z0_A1[j,] < 0.05, 0.05, G_tau_Z0_A1[j,])
      #G_tau_Z0_A0[j,] = ifelse(G_tau_Z0_A0[j,] < 0.05, 0.05, G_tau_Z0_A0[j,])
    
    
      for(t in 2:j){
          D_Z1_A1[j,] = D_Z1_A1[j,] +
            - pi_1*(dat$Z*dat$A*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z1_A1[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                 dat$R == 1) - h_Z1_A1[t,])/(S_tau_Z1_A1[t, ]*G_tau_Z1_A1[t-1,]*pi_1*delta_1))
          D_Z1_A0[j,] = D_Z1_A0[j,] +
            - (1-pi_1)*(dat$Z*(1-dat$A)*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z1_A0[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                 dat$R == 1) - h_Z1_A0[t,])/(S_tau_Z1_A0[t, ]*G_tau_Z1_A0[t-1,]*(1-pi_1)*delta_1))
          D_Z0_A1[j,] = D_Z0_A1[j,] +
            - pi_0*((1-dat$Z)*dat$A*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z0_A1[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                     dat$R == 1) - h_Z0_A1[t,])/(S_tau_Z0_A1[t, ]*G_tau_Z0_A1[t-1,]*pi_0*delta_0))
          D_Z0_A0[j,] = D_Z0_A0[j,] +
            - (1-pi_0)*((1-dat$Z)*(1-dat$A)*as.integer(time.list[t-1] < dat$obs.Y)*S_tau_Z0_A0[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                     dat$R == 1) - h_Z0_A0[t,])/(S_tau_Z0_A0[t, ]*G_tau_Z0_A0[t-1,]*(1-pi_0)*delta_0))  
       }
   
     
       D_Z1_A1[j,] = D_Z1_A1[j,] + S_tau_Z1_A1[j,]*dat$Z*(dat$A - pi_1)/delta_1 + S_tau_Z1_A1[j,]*pi_1
       D_Z1_A0[j,] = D_Z1_A0[j,] + S_tau_Z1_A0[j,]*dat$Z*((1-dat$A) - (1-pi_1))/delta_1 + S_tau_Z1_A0[j,]*(1-pi_1)
       D_Z0_A1[j,] = D_Z0_A1[j,] + S_tau_Z0_A1[j,]*(1-dat$Z)*(dat$A -pi_0)/delta_0 + S_tau_Z0_A1[j,]*pi_0
       D_Z0_A0[j,] = D_Z0_A0[j,] + S_tau_Z0_A0[j,]*(1-dat$Z)*((1-dat$A) -(1-pi_0))/delta_0 + S_tau_Z0_A0[j,]*(1-pi_0)


      ifvals[, j] = mean( (D_Z1_A1[j,] + D_Z1_A0[j, ]) - (D_Z0_A1[j,] + D_Z0_A0[j,]), na.rm = TRUE) / 
                    mean(phi2_1 - phi2_0, na.rm = TRUE) 
      est.eff[j] = mean(ifvals[,j]) # for each delta vaule              
      est.var[j] = var( ((D_Z1_A1[j,] + D_Z1_A0[j, ]) - (D_Z0_A1[j,] + D_Z0_A0[j,]) - est.eff[j]*(phi2_1-phi2_0))/ mean(phi2_1 - phi2_0) ) / length(phi2_1)
    }



   return(list(est.eff = est.eff, est.var = est.var))
}


binary_survival_para_naive = function(dat, time.list, misspecify = NULL, survival = "Cox"){
  
  n = nrow(dat)
  k = length(time.list) 
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)

  subdat = dat
  # fit a treatment model
  if("treatment" %in% misspecify){
      pi_1 = pi_fun_naive_wrong(subdat = subdat, dat = dat)
  }else{
      pi_1 = pi_fun_naive(subdat = subdat, dat = dat)
  } 
    
   S_tau_A1 = S_tau_A0 = G_tau_A1 = G_tau_A0 = matrix(0, length(time.list), nrow(dat))
   h_A1 = h_A0 = omega_A1 = omega_A0 = matrix(0, length(time.list), nrow(dat))
   D_A1 = D_A0 = matrix(0, length(time.list), nrow(dat))
   ## when j = 1
   S_tau_A1[1,] = S_tau_A0[1,] = rep(1, nrow(dat))
   G_tau_A1[1,] = G_tau_A0[1,] = rep(1, nrow(dat))
   D_A1[1,] = D_A0[1,] = rep(1, nrow(dat))
   ifvals[, 1] = mean(D_A1[1,] - D_A0[1,], na.rm = TRUE) 

   if("hazard" %in% misspecify){
      newdata_A1 = data.frame(W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)))
      newdata_A0 = data.frame(W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)))

      if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y, R) ~ W.1 + W.2 + W.3 + W.4 + W.5 + A, data = dat)
        basefit = basehaz(coxfit, centered = FALSE)
        hazardfit = matrix(0, nrow(basefit), 2)
        hazardfit[,2] = basefit$time
        hazardfit[1,1] = basefit$hazard[1]*(basefit$time[1] - 0)
        for(t in 2:nrow(hazardfit)){
            hazardfit[t,1] = (basefit$hazard[t]-basefit$hazard[t-1])*(basefit$time[t] - basefit$time[t-1])
        }
    
      }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$W.1, dat$W.2, dat$W.3, dat$W.4, dat$W.5, dat$A)
        colnames(design) = c("W.1", "W.2", "W.3", "W.4", "W.5", "A")
        addfit = ahaz(Surv(dat$obs.YY, dat$R), design)
        addbase = predict(addfit, type = "cumhaz")
        addbasehazard = matrix(0, length(addbase$times), 2)
        addbasehazard [,2] = addbase$times
        addbasehazard [1,1] = addbase$cumhaz[1]*(addbase$times[1] - 0)
         for(t in 2:nrow(addbasehazard )){
           addbasehazard[t,1] = (addbase$cumhaz[t] - addbase$cumhaz[t-1])*(addbase$times[t] - addbase$times[t-1])
        }
      } 
  }else{
    newdata_A1 = data.frame(X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(1, nrow(dat)))
    newdata_A0 = data.frame(X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                         A = rep(0, nrow(dat)))
    
    if(survival == "Cox"){
        coxfit = coxph(Surv(obs.Y, R) ~ X.1 + X.2 + X.3 + X.4 + X.5 + A, data = dat)
         basefit = basehaz(coxfit, centered = FALSE)
         #hazardfit = matrix(0, nrow(basefit), 2)
         #hazardfit[,2] = basefit$time
         #hazardfit[1,1] = basefit$hazard[1]*(basefit$time[1] - 0)
         #for(t in 2:nrow(hazardfit)){
         #   hazardfit[t,1] = (basefit$hazard[t]-basefit$hazard[t-1])*(basefit$time[t] - basefit$time[t-1])
         #}#
    }else if(survival == "additive"){
        dat$obs.YY = dat$obs.Y + runif(dim(dat)[1], 0, 1)*10^(-3)
        design = cbind(dat$X.1, dat$X.2, dat$X.3, dat$X.4, dat$X.5, dat$A)
        colnames(design) = c("X.1", "X.2", "X.3", "X.4", "X.5", "A")
        addfit = ahaz(Surv(dat$obs.YY, dat$R), design)
        addbase = predict(addfit, type = "cumhaz")
        
        tau = as.integer(max(dat$obs.Y)); tmp.cum = matrix(0, tau, 2)
        for(t in 1:tau){
         dummy = which(as.integer(addbase$times) == t)[1]
         tmp.cum[t,] = c(addbase$cumhaz[dummy], addbase$times[dummy])
        }
        tmp.cum = na.omit(tmp.cum)
        addbasehazard = matrix(0, nrow(tmp.cum), 2)
        addbasehazard[,2] = as.integer(tmp.cum[,2])
        addbasehazard[1,1] = tmp.cum[1,1]*(tmp.cum[1,2] - 0)
        for(t in 2:nrow(tmp.cum)){
           addbasehazard[t,1] = (tmp.cum[t,1]-tmp.cum[t-1,1])*(tmp.cum[t,2] - tmp.cum[t-1,2])
        }
    } 
  }  

  
   for(j in 2:k){
      subdat$Y = as.integer(subdat$obs.Y == time.list[j] & subdat$R == 1)
      subdat$C = as.integer(subdat$obs.Y  == time.list[j] & subdat$R == 0)
      subdat$I = as.integer(time.list[j-1]  < subdat$obs.Y) 

      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(dat$obs.Y == time.list[j]  & dat$R == 0)
      dat$I = as.integer(time.list[j-1]  < dat$obs.Y)

      if(survival == "Cox"){
          baselines.h1 = basefit$hazard[which.min(abs(basefit$time - time.list[j]))] # baseline
          baselines.h0 = basefit$hazard[which.min(abs(basefit$time - time.list[j-1]))] # baseline
          h_A1[j,] = 1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(coxfit)$coefficients[1:6,1])))) /
          exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(coxfit)$coefficients[1:6,1])))) 
                       
          h_A0[j,] =  1 - exp(-(baselines.h1)*exp(as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(coxfit)$coefficients[1:6,1])))) /
          exp(-(baselines.h0)*exp(as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(coxfit)$coefficients[1:6,1])))) 
        }else if(survival == "additive"){
          baselines.h = addbasehazard[which.min(abs( addbasehazard[,2] - time.list[j])),1] # baseline
          h_A1[j,] = baselines.h + as.numeric(as.matrix(newdata_A1) %*% as.matrix(summary(addfit)$coefficients[1:6,1]))
          h_A0[j,] = baselines.h + as.numeric(as.matrix(newdata_A0) %*% as.matrix(summary(addfit)$coefficients[1:6,1]))
      } 

      S_tau_A1[j,] = S_tau_A1[j-1,]*(1-h_A1[j,])
      S_tau_A0[j,] = S_tau_A0[j-1,]*(1-h_A0[j,])

      S_tau_A1[j,] = ifelse(S_tau_A1[j,] < 0.05, 0.05, S_tau_A1[j,])
      S_tau_A0[j,] = ifelse(S_tau_A0[j,] < 0.05, 0.05, S_tau_A0[j,])


      if("censoring" %in% misspecify){
          omega_A1[j,] = censor_fun_naive_wrong(I = 1, A = 1, subdat = subdat, dat = dat)
          omega_A0[j,] = censor_fun_naive_wrong(I = 1, A = 0, subdat = subdat, dat = dat)
      }else{
          omega_A1[j,] = censor_fun_naive(I = 1, A = 1, subdat = subdat, dat = dat)
          omega_A0[j,] = censor_fun_naive(I = 1, A = 0, subdat = subdat, dat = dat)
      }

      G_tau_A1[j,] = G_tau_A1[j-1,]*(1-omega_A1[j,])
      G_tau_A0[j,] = G_tau_A0[j-1,]*(1-omega_A0[j,])

      #G_tau_A1[j,] = ifelse(G_tau_A1[j,] < 0.05, 0.05, G_tau_A1[j,])
      #G_tau_A0[j,] = ifelse(G_tau_A0[j,] < 0.05, 0.05, G_tau_A0[j,])

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
      est.eff[j] = mean(ifvals[,j]) # for each delta vaule              
      est.var[j] = var(D_A1[j,] - D_A0[j,]) / nrow(dat)
   }
    

  return(list(est.eff = est.eff, est.var = est.var))
}