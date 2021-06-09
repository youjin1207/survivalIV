#source("binary_survival_para_nuisance.R")
binary_survival_ranger_if = function(dat, x.trt, x.out, x.miss, x.instrument, 
                                         w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  
  x.trt_Z1 = x.trt_Z0 = x.trt
  x.trt_Z1$Z = rep(1, nrow(x.trt_Z1)); x.trt_Z0$Z = rep(0, nrow(x.trt_Z0))
  
  x.out_Z1_A1 = x.out_Z1_A0 = x.out_Z0_A1 = x.out_Z0_A0 = x.out
  x.out_Z1_A1$Z = rep(1, nrow(x.out)); x.out_Z1_A1$A = rep(1, nrow(x.out))
  x.out_Z1_A0$Z = rep(1, nrow(x.out)); x.out_Z1_A0$A = rep(0, nrow(x.out))
  x.out_Z0_A1$Z = rep(0, nrow(x.out)); x.out_Z0_A1$A = rep(1, nrow(x.out))
  x.out_Z0_A0$Z = rep(0, nrow(x.out)); x.out_Z0_A0$A = rep(0, nrow(x.out))
  
  x.miss_Z1_A1 = x.miss_Z1_A0 = x.miss_Z0_A1 = x.miss_Z0_A0 = x.miss
  x.miss_Z1_A1$Z = rep(1, nrow(x.miss)); x.miss_Z1_A1$A = rep(1, nrow(x.miss))
  x.miss_Z1_A0$Z = rep(1, nrow(x.miss)); x.miss_Z1_A0$A = rep(0, nrow(x.miss))
  x.miss_Z0_A1$Z = rep(0, nrow(x.miss)); x.miss_Z0_A1$A = rep(1, nrow(x.miss))
  x.miss_Z0_A0$Z = rep(0, nrow(x.miss)); x.miss_Z0_A0$A = rep(0, nrow(x.miss))
  
  ## dat: data.frame(obs.T = obs.T, A = A, X = X, U = U, Z = Z)
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s
  tmp.var = matrix(NA, nsplits, k)
  for(split in 1:nsplits){
    # within 9 groups
    print(paste("split", split));
    flush.console() 
    
    # fit treatment model
    if("treatment" %in% misspecify){
      w.trt_Z1 = w.trt_Z0 = w.trt
      w.trt_Z1$Z = rep(1, nrow(w.trt_Z1)); w.trt_Z0$Z = rep(0, nrow(w.trt_Z0)) 
      trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = w.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = w.trt_Z0)$predictions 
      
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }else{
      trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = x.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = x.trt_Z0)$predictions  
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }
    
    if("instrument" %in% misspecify){
      instmod = ranger(Z ~ . , dat = cbind(w.instrument, Z = dat$Z)[slong!=split,])
      delta_1 = predict(instmod, data = w.instrument)$predictions
      delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
      delta_0 = (1-delta_1)
    }else{
      instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,])
      delta_1 = predict(instmod, data = x.instrument)$predictions
      delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
      delta_0 = (1-delta_1)
    }
    
    
    if("censoring" %in% misspecify){
      w.miss_Z1_A1 = w.miss_Z1_A0 = w.miss_Z0_A1 = w.miss_Z0_A0 = w.miss
      w.miss_Z1_A1$Z = rep(1, nrow(w.miss)); w.miss_Z1_A1$A = rep(1, nrow(w.miss))
      w.miss_Z1_A0$Z = rep(1, nrow(w.miss)); w.miss_Z1_A0$A = rep(0, nrow(w.miss))
      w.miss_Z0_A1$Z = rep(0, nrow(w.miss)); w.miss_Z0_A1$A = rep(1, nrow(w.miss))
      w.miss_Z0_A0$Z = rep(0, nrow(w.miss)); w.miss_Z0_A0$A = rep(0, nrow(w.miss))
      
      missmod = ranger(R ~., data = cbind(w.miss, R = dat$R)[slong!=split,])
      omega_Z1_A1 = predict(missmod, data = w.miss_Z1_A1)$predictions
      omega_Z1_A0 = predict(missmod, data = w.miss_Z1_A0)$predictions
      omega_Z0_A1 = predict(missmod, data = w.miss_Z0_A1)$predictions
      omega_Z0_A0 = predict(missmod, data = w.miss_Z0_A0)$predictions
    }else{
      missmod = ranger(R ~., data = cbind(x.miss, R = dat$R)[slong!=split,])
      omega_Z1_A1 = predict(missmod, data = x.miss_Z1_A1)$predictions
      omega_Z1_A0 = predict(missmod, data = x.miss_Z1_A0)$predictions
      omega_Z0_A1 = predict(missmod, data = x.miss_Z0_A1)$predictions
      omega_Z0_A0 = predict(missmod, data = x.miss_Z0_A0)$predictions
    }
    
    
    ## fit outcome models (time-varying)
    for(j in 1:k){
      dat$Y = as.integer(time.list[j] < dat$obs.Y)
      
      if("outcome" %in% misspecify){
        w.out_Z1_A1 = w.out_Z1_A0 = w.out_Z0_A1 = w.out_Z0_A0 = w.out
        w.out_Z1_A1$Z = rep(1, nrow(w.out)); w.out_Z1_A1$A = rep(1, nrow(w.out))
        w.out_Z1_A0$Z = rep(1, nrow(w.out)); w.out_Z1_A0$A = rep(0, nrow(w.out))
        w.out_Z0_A1$Z = rep(0, nrow(w.out)); w.out_Z0_A1$A = rep(1, nrow(w.out))
        w.out_Z0_A0$Z = rep(0, nrow(w.out)); w.out_Z0_A0$A = rep(0, nrow(w.out))
        
        outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
        mu_Z1_A1 = predict(outmod, data = w.out_Z1_A1)$predictions
        mu_Z1_A0 = predict(outmod, data = w.out_Z1_A0)$predictions
        mu_Z0_A1 = predict(outmod, data = w.out_Z0_A1)$predictions
        mu_Z0_A0 = predict(outmod, data = w.out_Z0_A0)$predictions
      }else{
        outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
        mu_Z1_A1 = predict(outmod, data = x.out_Z1_A1)$predictions
        mu_Z1_A0 = predict(outmod, data = x.out_Z1_A0)$predictions
        mu_Z0_A1 = predict(outmod, data = x.out_Z0_A1)$predictions
        mu_Z0_A0 = predict(outmod, data = x.out_Z0_A0)$predictions
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
      ifvals[s == split, j] = mean(phi1_1[s==split] - phi1_0[s==split]) / 
        mean(phi2_1[s==split] - phi2_0[s==split]) 
      tmp.var[split,j] = var( (phi1_1[s==split] - phi1_0[s==split] - mean(ifvals[s==split, j])*(phi2_1[s==split]-phi2_0[s==split]))/ mean(phi2_1[s==split] - phi2_0[s==split]) ) / nrow(dat[s==split,])  
    }
  }
  
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
    est.var[j] = sum(tmp.var[,j])/(nsplits)^2
  }
  
  return(list(est.eff = est.eff, est.var = est.var))
}


binary_survival_ranger_if_hazard = function(dat, x.trt, x.out, x.miss, x.instrument, 
                                                w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  
  x.trt_Z1 = x.trt_Z0 = x.trt
  x.trt_Z1$Z = rep(1, nrow(x.trt_Z1)); x.trt_Z0$Z = rep(0, nrow(x.trt_Z0))
  
  x.out_Z1_A1 = x.out_Z1_A0 = x.out_Z0_A1 = x.out_Z0_A0 = x.out
  x.out_Z1_A1$Z = rep(1, nrow(x.out)); x.out_Z1_A1$A = rep(1, nrow(x.out))
  x.out_Z1_A0$Z = rep(1, nrow(x.out)); x.out_Z1_A0$A = rep(0, nrow(x.out))
  x.out_Z0_A1$Z = rep(0, nrow(x.out)); x.out_Z0_A1$A = rep(1, nrow(x.out))
  x.out_Z0_A0$Z = rep(0, nrow(x.out)); x.out_Z0_A0$A = rep(0, nrow(x.out))
  
  x.miss_Z1_A1 = x.miss_Z1_A0 = x.miss_Z0_A1 = x.miss_Z0_A0 = x.miss
  x.miss_Z1_A1$Z = rep(1, nrow(x.miss)); x.miss_Z1_A1$A = rep(1, nrow(x.miss))
  x.miss_Z1_A0$Z = rep(1, nrow(x.miss)); x.miss_Z1_A0$A = rep(0, nrow(x.miss))
  x.miss_Z0_A1$Z = rep(0, nrow(x.miss)); x.miss_Z0_A1$A = rep(1, nrow(x.miss))
  x.miss_Z0_A0$Z = rep(0, nrow(x.miss)); x.miss_Z0_A0$A = rep(0, nrow(x.miss))
  
  ## dat: data.frame(obs.T = obs.T, A = A, X = X, U = U, Z = Z)
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s
  tmp.var = matrix(NA, nsplits, k)
  
  S_tau_Z1_A1 = S_tau_Z1_A0 = S_tau_Z0_A1 = S_tau_Z0_A0 = 
    G_tau_Z1_A1 = G_tau_Z1_A0 = G_tau_Z0_A1 = G_tau_Z0_A0 = matrix(1, length(time.list), nrow(dat))
  h_Z1_A1 =  h_Z1_A0 =  h_Z0_A1 = h_Z0_A0 = 
    omega_Z1_A1 = omega_Z1_A0 = omega_Z0_A1 = omega_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  D_Z1_A1 = D_Z1_A0 = D_Z0_A1 = D_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  
  
  for(split in 1:nsplits){
    # within 9 groups
    print(paste("split", split));
    flush.console() 
    
    # fit a treatment model
    if("treatment" %in% misspecify){
      w.trt_Z1 = w.trt_Z0 = w.trt
      w.trt_Z1$Z = rep(1, nrow(w.trt_Z1)); w.trt_Z0$Z = rep(0, nrow(w.trt_Z0)) 
      trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = w.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = w.trt_Z0)$predictions 
      
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }else{
      trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = x.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = x.trt_Z0)$predictions  
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }
    
    
    # fit an instrumental variable model
    if("instrument" %in% misspecify){
      instmod = ranger(Z ~ . , dat = cbind(w.instrument, Z = dat$Z)[slong!=split,])
      delta_1 = predict(instmod, data = w.instrument)$predictions
      delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
      delta_0 = (1-delta_1)
    }else{
      instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,])
      delta_1 = predict(instmod, data = x.instrument)$predictions
      delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
      delta_0 = (1-delta_1)
    }
    
    phi2_1 = pi_1 + 
      dat$Z*(dat$A - pi_1)/delta_1
    phi2_0 = pi_0  + 
      (1-dat$Z)*(dat$A - pi_0)/delta_0
    
    
    
    
    ## fit outcome models (time-varying)
    for(j in 2:k){
      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(time.list[j] == dat$obs.Y & dat$R == 0)
      dat$I = as.integer(time.list[j-1] < dat$obs.Y)
      
      if("outcome" %in% misspecify){
        w.out_Z1_A1 = w.out_Z1_A0 = w.out_Z0_A1 = w.out_Z0_A0 = w.out
        w.out_Z1_A1$Z = rep(1, nrow(w.out)); w.out_Z1_A1$A = rep(1, nrow(w.out))
        w.out_Z1_A0$Z = rep(1, nrow(w.out)); w.out_Z1_A0$A = rep(0, nrow(w.out))
        w.out_Z0_A1$Z = rep(0, nrow(w.out)); w.out_Z0_A1$A = rep(1, nrow(w.out))
        w.out_Z0_A0$Z = rep(0, nrow(w.out)); w.out_Z0_A0$A = rep(0, nrow(w.out))
        w.out$I = dat$I;
        w.out_Z1_A1$I = w.out_Z1_A0$I = w.out_Z0_A1$I = w.out_Z0_A0$I = rep(1, nrow(dat))
        
        outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
        h_Z1_A1[j,slong==split] = predict(outmod, data = w.out_Z1_A1)$predictions[slong==split]
        h_Z1_A0[j,slong==split] = predict(outmod, data = w.out_Z1_A0)$predictions[slong==split]
        h_Z0_A1[j,slong==split] = predict(outmod, data = w.out_Z0_A1)$predictions[slong==split]
        h_Z0_A0[j,slong==split] = predict(outmod, data = w.out_Z0_A0)$predictions[slong==split]
      }else{
        x.out$I = dat$I;
        x.out_Z1_A1$I = x.out_Z1_A0$I = x.out_Z0_A1$I = x.out_Z0_A0$I = rep(1, nrow(dat))
        outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
        h_Z1_A1[j,slong==split] = predict(outmod, data = x.out_Z1_A1)$predictions[slong==split]
        h_Z1_A0[j,slong==split] = predict(outmod, data = x.out_Z1_A0)$predictions[slong==split]
        h_Z0_A1[j,slong==split] = predict(outmod, data = x.out_Z0_A1)$predictions[slong==split]
        h_Z0_A0[j,slong==split] = predict(outmod, data = x.out_Z0_A0)$predictions[slong==split]
      }
      
      if("censoring" %in% misspecify){
        w.miss_Z1_A1 = w.miss_Z1_A0 = w.miss_Z0_A1 = w.miss_Z0_A0 = w.miss
        w.miss_Z1_A1$Z = rep(1, nrow(w.miss)); w.miss_Z1_A1$A = rep(1, nrow(w.miss))
        w.miss_Z1_A0$Z = rep(1, nrow(w.miss)); w.miss_Z1_A0$A = rep(0, nrow(w.miss))
        w.miss_Z0_A1$Z = rep(0, nrow(w.miss)); w.miss_Z0_A1$A = rep(1, nrow(w.miss))
        w.miss_Z0_A0$Z = rep(0, nrow(w.miss)); w.miss_Z0_A0$A = rep(0, nrow(w.miss))
        w.miss$I = dat$I;
        w.miss_Z1_A1$I = w.miss_Z1_A0$I = w.miss_Z0_A1$I = w.miss_Z0_A0$I = rep(1, nrow(dat))
        
        missmod = ranger(C ~., data = cbind(w.miss, C = dat$C)[slong!=split,])
        omega_Z1_A1[j,slong==split] = predict(missmod, data = w.miss_Z1_A1)$predictions[slong==split]
        omega_Z1_A0[j,slong==split] = predict(missmod, data = w.miss_Z1_A0)$predictions[slong==split]
        omega_Z0_A1[j,slong==split] = predict(missmod, data = w.miss_Z0_A1)$predictions[slong==split]
        omega_Z0_A0[j,slong==split] = predict(missmod, data = w.miss_Z0_A0)$predictions[slong==split]
      }else{
        x.miss$I = dat$I;
        x.miss_Z1_A1$I = x.miss_Z1_A0$I = x.miss_Z0_A1$I = x.miss_Z0_A0$I = rep(1, nrow(dat))
        missmod = ranger(C ~., data = cbind(x.miss, C = dat$C)[slong!=split,])
        omega_Z1_A1[j,slong==split] = predict(missmod, data = x.miss_Z1_A1)$predictions[slong==split]
        omega_Z1_A0[j,slong==split] = predict(missmod, data = x.miss_Z1_A0)$predictions[slong==split]
        omega_Z0_A1[j,slong==split] = predict(missmod, data = x.miss_Z0_A1)$predictions[slong==split]
        omega_Z0_A0[j,slong==split] = predict(missmod, data = x.miss_Z0_A0)$predictions[slong==split]
      }
      
      
      S_tau_Z1_A1[j,slong==split] = S_tau_Z1_A1[j-1,slong==split]*(1-h_Z1_A1[j,slong==split])
      S_tau_Z1_A0[j,slong==split] = S_tau_Z1_A0[j-1,slong==split]*(1-h_Z1_A0[j,slong==split])
      S_tau_Z0_A1[j,slong==split] = S_tau_Z0_A1[j-1,slong==split]*(1-h_Z0_A1[j,slong==split])
      S_tau_Z0_A0[j,slong==split] = S_tau_Z0_A0[j-1,slong==split]*(1-h_Z0_A0[j,slong==split])
      
      
      G_tau_Z1_A1[j,slong==split] = G_tau_Z1_A1[j-1,slong==split]*(1-omega_Z1_A1[j,slong==split])
      G_tau_Z1_A0[j,slong==split] = G_tau_Z1_A0[j-1,slong==split]*(1-omega_Z1_A0[j,slong==split])
      G_tau_Z0_A1[j,slong==split] = G_tau_Z0_A1[j-1,slong==split]*(1-omega_Z0_A1[j,slong==split])
      G_tau_Z0_A0[j,slong==split] = G_tau_Z0_A0[j-1,slong==split]*(1-omega_Z0_A0[j,slong==split])
      
      
      S_tau_Z1_A1[j,slong==split] = ifelse(S_tau_Z1_A1[j,slong==split] < 0.05, 0.05, S_tau_Z1_A1[j,slong==split])
      S_tau_Z1_A0[j,slong==split] = ifelse(S_tau_Z1_A0[j,slong==split] < 0.05, 0.05, S_tau_Z1_A0[j,slong==split])
      S_tau_Z0_A1[j,slong==split] = ifelse(S_tau_Z0_A1[j,slong==split] < 0.05, 0.05, S_tau_Z0_A1[j,slong==split])
      S_tau_Z0_A0[j,slong==split] = ifelse(S_tau_Z0_A0[j,slong==split] < 0.05, 0.05, S_tau_Z0_A0[j,slong==split])
      
      
      G_tau_Z1_A1[j,slong==split] = ifelse(G_tau_Z1_A1[j,slong==split] < 0.05, 0.05, G_tau_Z1_A1[j,slong==split])
      G_tau_Z1_A0[j,slong==split] = ifelse(G_tau_Z1_A0[j,slong==split] < 0.05, 0.05, G_tau_Z1_A0[j,slong==split])
      G_tau_Z0_A1[j,slong==split] = ifelse(G_tau_Z0_A1[j,slong==split] < 0.05, 0.05, G_tau_Z0_A1[j,slong==split])
      G_tau_Z0_A0[j,slong==split] = ifelse(G_tau_Z0_A0[j,slong==split] < 0.05, 0.05, G_tau_Z0_A0[j,slong==split])
      
      
      for(t in 2:j){
        D_Z1_A1[j,slong==split] = D_Z1_A1[j,slong==split] +
          - pi_1[slong==split]*(dat$Z[slong==split]*dat$A[slong==split]*as.integer(time.list[t-1] < dat$obs.Y[slong==split])*S_tau_Z1_A1[j,slong==split]*(as.integer(dat$obs.Y[slong==split] == time.list[t] &
                                                                                                                                                                       dat$R[slong==split] == 1) - h_Z1_A1[t,slong==split])/(S_tau_Z1_A1[t, slong==split]*G_tau_Z1_A1[t-1,slong==split]*pi_1[slong==split]*delta_1[slong==split]))
        D_Z1_A0[j,slong==split] = D_Z1_A0[j,slong==split] +
          - (1-pi_1[slong==split])*(dat$Z[slong==split]*(1-dat$A[slong==split])*as.integer(time.list[t-1] < dat$obs.Y[slong==split])*S_tau_Z1_A0[j,slong==split]*(as.integer(dat$obs.Y[slong==split] == time.list[t] &
                                                                                                                                                                               dat$R[slong==split] == 1) - h_Z1_A0[t,slong==split])/(S_tau_Z1_A0[t, slong==split]*G_tau_Z1_A0[t-1,slong==split]*(1-pi_1[slong==split])*delta_1[slong==split]))
        D_Z0_A1[j,slong==split] = D_Z0_A1[j,slong==split] +
          - pi_0[slong==split]*((1-dat$Z[slong==split])*dat$A[slong==split]*as.integer(time.list[t-1] < dat$obs.Y[slong==split])*S_tau_Z0_A1[j,slong==split]*(as.integer(dat$obs.Y[slong==split] == time.list[t] &
                                                                                                                                                                           dat$R[slong==split] == 1) - h_Z0_A1[t,slong==split])/(S_tau_Z0_A1[t, slong==split]*G_tau_Z0_A1[t-1,slong==split]*pi_0[slong==split]*delta_0[slong==split]))
        D_Z0_A0[j,slong==split] = D_Z0_A0[j,slong==split] +
          - (1-pi_0[slong==split])*((1-dat$Z[slong==split])*(1-dat$A[slong==split])*as.integer(time.list[t-1] < dat$obs.Y[slong==split])*S_tau_Z0_A0[j,slong==split]*(as.integer(dat$obs.Y[slong==split] == time.list[t] &
                                                                                                                                                                                   dat$R[slong==split] == 1) - h_Z0_A0[t,slong==split])/(S_tau_Z0_A0[t, slong==split]*G_tau_Z0_A0[t-1,slong==split]*(1-pi_0[slong==split])*delta_0[slong==split]))  
      }
      
      
      D_Z1_A1[j,slong==split] = D_Z1_A1[j,slong==split] + S_tau_Z1_A1[j,slong==split]*dat$Z[slong==split]*(dat$A[slong==split] - pi_1[slong==split])/delta_1[slong==split] + S_tau_Z1_A1[j,slong==split]*pi_1[slong==split]
      D_Z1_A0[j,slong==split] = D_Z1_A0[j,slong==split] + S_tau_Z1_A0[j,slong==split]*dat$Z[slong==split]*((1-dat$A[slong==split]) - (1-pi_1[slong==split]))/delta_1[slong==split] + S_tau_Z1_A0[j,slong==split]*(1-pi_1[slong==split])
      D_Z0_A1[j,slong==split] = D_Z0_A1[j,slong==split] + S_tau_Z0_A1[j,slong==split]*(1-dat$Z[slong==split])*(dat$A[slong==split] -pi_0[slong==split])/delta_0[slong==split] + S_tau_Z0_A1[j,slong==split]*pi_0[slong==split]
      D_Z0_A0[j,slong==split] = D_Z0_A0[j,slong==split] + S_tau_Z0_A0[j,slong==split]*(1-dat$Z[slong==split])*((1-dat$A[slong==split]) -(1-pi_0[slong==split]))/delta_0[slong==split] + S_tau_Z0_A0[j,slong==split]*(1-pi_0[slong==split])
      
      
      ifvals[slong==split, j] = mean( (D_Z1_A1[j,slong==split] + D_Z1_A0[j, slong==split]) - (D_Z0_A1[j,slong==split] + D_Z0_A0[j,slong==split]), na.rm = TRUE) / 
        mean(phi2_1[slong==split] - phi2_0[slong==split], na.rm = TRUE) 
      tmp.var[split,j] =  var( ((D_Z1_A1[j,slong==split] + D_Z1_A0[j, slong==split]) - (D_Z0_A1[j,slong==split] + D_Z0_A0[j,slong==split]) - 
                                  mean(ifvals[slong==split, j])*(phi2_1[slong==split]-phi2_0[slong==split]))/ mean(phi2_1[slong==split] - phi2_0[slong==split]) ) / nrow(dat[slong==split,])
    }
  } 
  
  
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) 
    est.var[j] = sum(tmp.var[,j])/(nsplits)^2
  }
  
  return(list(est.eff = est.eff, est.var = est.var))
}


binary_survival_ranger_hazard_plugin = function(dat, x.trt, x.out, x.miss, x.instrument, 
                                                w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  
  x.trt_Z1 = x.trt_Z0 = x.trt
  x.trt_Z1$Z = rep(1, nrow(x.trt_Z1)); x.trt_Z0$Z = rep(0, nrow(x.trt_Z0))
  
  x.out_Z1_A1 = x.out_Z1_A0 = x.out_Z0_A1 = x.out_Z0_A0 = x.out
  x.out_Z1_A1$Z = rep(1, nrow(x.out)); x.out_Z1_A1$A = rep(1, nrow(x.out))
  x.out_Z1_A0$Z = rep(1, nrow(x.out)); x.out_Z1_A0$A = rep(0, nrow(x.out))
  x.out_Z0_A1$Z = rep(0, nrow(x.out)); x.out_Z0_A1$A = rep(1, nrow(x.out))
  x.out_Z0_A0$Z = rep(0, nrow(x.out)); x.out_Z0_A0$A = rep(0, nrow(x.out))
  
  x.miss_Z1_A1 = x.miss_Z1_A0 = x.miss_Z0_A1 = x.miss_Z0_A0 = x.miss
  x.miss_Z1_A1$Z = rep(1, nrow(x.miss)); x.miss_Z1_A1$A = rep(1, nrow(x.miss))
  x.miss_Z1_A0$Z = rep(1, nrow(x.miss)); x.miss_Z1_A0$A = rep(0, nrow(x.miss))
  x.miss_Z0_A1$Z = rep(0, nrow(x.miss)); x.miss_Z0_A1$A = rep(1, nrow(x.miss))
  x.miss_Z0_A0$Z = rep(0, nrow(x.miss)); x.miss_Z0_A0$A = rep(0, nrow(x.miss))
  
  ## dat: data.frame(obs.T = obs.T, A = A, X = X, U = U, Z = Z)
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(nrow = n, ncol = k)
  est.eff = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s
  
  for(split in 1:nsplits){
    # within 9 groups
    print(paste("split", split));
    flush.console() 
    
    # fit treatment model
    if("treatment" %in% misspecify){
      w.trt_Z1 = w.trt_Z0 = w.trt
      w.trt_Z1$Z = rep(1, nrow(w.trt_Z1)); w.trt_Z0$Z = rep(0, nrow(w.trt_Z0)) 
      trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = w.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = w.trt_Z0)$predictions 
      
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }else{
      trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = x.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = x.trt_Z0)$predictions  
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }
    
    S_tau_Z1_A1 = S_tau_Z1_A0 = S_tau_Z0_A1 = S_tau_Z0_A0 =  matrix(1, length(time.list), nrow(dat))
    h_Z1_A1 =  h_Z1_A0 =  h_Z0_A1 = h_Z0_A0 = matrix(0, length(time.list), nrow(dat))
    
    
    
    ## fit outcome models (time-varying)
    for(j in 2:k){
      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(time.list[j] == dat$obs.Y & dat$R == 0)
      dat$I = as.integer(time.list[j-1] < dat$obs.Y)
      
      if("outcome" %in% misspecify){
        w.out_Z1_A1 = w.out_Z1_A0 = w.out_Z0_A1 = w.out_Z0_A0 = w.out
        w.out_Z1_A1$Z = rep(1, nrow(w.out)); w.out_Z1_A1$A = rep(1, nrow(w.out))
        w.out_Z1_A0$Z = rep(1, nrow(w.out)); w.out_Z1_A0$A = rep(0, nrow(w.out))
        w.out_Z0_A1$Z = rep(0, nrow(w.out)); w.out_Z0_A1$A = rep(1, nrow(w.out))
        w.out_Z0_A0$Z = rep(0, nrow(w.out)); w.out_Z0_A0$A = rep(0, nrow(w.out))
        w.out$I = dat$I;
        w.out_Z1_A1$I = w.out_Z1_A0$I = w.out_Z0_A1$I = w.out_Z0_A0$I = rep(1, nrow(dat))
        
        outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
        h_Z1_A1[j,] = predict(outmod, data = w.out_Z1_A1)$predictions
        h_Z1_A0[j,] = predict(outmod, data = w.out_Z1_A0)$predictions
        h_Z0_A1[j,] = predict(outmod, data = w.out_Z0_A1)$predictions
        h_Z0_A0[j,] = predict(outmod, data = w.out_Z0_A0)$predictions
      }else{
        x.out$I = dat$I;
        x.out_Z1_A1$I = x.out_Z1_A0$I = x.out_Z0_A1$I = x.out_Z0_A0$I = rep(1, nrow(dat))
        outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
        h_Z1_A1[j,] = predict(outmod, data = x.out_Z1_A1)$predictions
        h_Z1_A0[j,] = predict(outmod, data = x.out_Z1_A0)$predictions
        h_Z0_A1[j,] = predict(outmod, data = x.out_Z0_A1)$predictions
        h_Z0_A0[j,] = predict(outmod, data = x.out_Z0_A0)$predictions
      }
      
      S_tau_Z1_A1[j,] = S_tau_Z1_A1[j-1,]*(1-h_Z1_A1[j,])
      S_tau_Z1_A0[j,] = S_tau_Z1_A0[j-1,]*(1-h_Z1_A0[j,])
      S_tau_Z0_A1[j,] = S_tau_Z0_A1[j-1,]*(1-h_Z0_A1[j,])
      S_tau_Z0_A0[j,] = S_tau_Z0_A0[j-1,]*(1-h_Z0_A0[j,])
      
      
      
      
      phi1_1 = S_tau_Z1_A1[j,]*pi_1  + S_tau_Z1_A0[j,]*(1-pi_1) 
      phi1_0 = S_tau_Z0_A1[j,]*pi_0  + S_tau_Z0_A0[j,]*(1-pi_0) 
      
      phi2_1 = pi_1 
      phi2_0 = pi_0  
      ## influence-function based estimator:
      ifvals[s==split, j] = mean(phi1_1[s==split] - phi1_0[s==split]) / 
        mean(phi2_1[s==split] - phi2_0[s==split]) 
    }
  }
  
  # compute estimator (using)
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
  }
  
  return(list(est = est.eff))
}



binary_survival_ranger_hazard_ipw = function(dat, x.trt, x.out, x.miss, x.instrument, 
                                             w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  
  x.trt_Z1 = x.trt_Z0 = x.trt
  x.trt_Z1$Z = rep(1, nrow(x.trt_Z1)); x.trt_Z0$Z = rep(0, nrow(x.trt_Z0))
  
  x.out_Z1_A1 = x.out_Z1_A0 = x.out_Z0_A1 = x.out_Z0_A0 = x.out
  x.out_Z1_A1$Z = rep(1, nrow(x.out)); x.out_Z1_A1$A = rep(1, nrow(x.out))
  x.out_Z1_A0$Z = rep(1, nrow(x.out)); x.out_Z1_A0$A = rep(0, nrow(x.out))
  x.out_Z0_A1$Z = rep(0, nrow(x.out)); x.out_Z0_A1$A = rep(1, nrow(x.out))
  x.out_Z0_A0$Z = rep(0, nrow(x.out)); x.out_Z0_A0$A = rep(0, nrow(x.out))
  
  x.miss_Z1_A1 = x.miss_Z1_A0 = x.miss_Z0_A1 = x.miss_Z0_A0 = x.miss
  x.miss_Z1_A1$Z = rep(1, nrow(x.miss)); x.miss_Z1_A1$A = rep(1, nrow(x.miss))
  x.miss_Z1_A0$Z = rep(1, nrow(x.miss)); x.miss_Z1_A0$A = rep(0, nrow(x.miss))
  x.miss_Z0_A1$Z = rep(0, nrow(x.miss)); x.miss_Z0_A1$A = rep(1, nrow(x.miss))
  x.miss_Z0_A0$Z = rep(0, nrow(x.miss)); x.miss_Z0_A0$A = rep(0, nrow(x.miss))
  
  ## dat: data.frame(obs.T = obs.T, A = A, X = X, U = U, Z = Z)
  n = nrow(dat)
  k = length(time.list) # for treatment density
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s
  
  G_tau_Z1_A1 = G_tau_Z1_A0 = G_tau_Z0_A1 = G_tau_Z0_A0 = matrix(1, length(time.list), nrow(dat))
  omega_Z1_A1 = omega_Z1_A0 = omega_Z0_A1 = omega_Z0_A0 = matrix(0, length(time.list), nrow(dat))
  
  for(split in 1:nsplits){
    # within 9 groups
    print(paste("split", split));
    flush.console() 
    
    # fit a treatment model
    if("treatment" %in% misspecify){
      w.trt_Z1 = w.trt_Z0 = w.trt
      w.trt_Z1$Z = rep(1, nrow(w.trt_Z1)); w.trt_Z0$Z = rep(0, nrow(w.trt_Z0)) 
      trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = w.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = w.trt_Z0)$predictions 
      
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }else{
      trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = x.trt_Z1)$predictions
      pi_0 = predict(trtmod, data = x.trt_Z0)$predictions  
      
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
      pi_0 = ifelse(pi_0 < 0.05, 0.05, ifelse(pi_0 > 0.95, 0.95, pi_0)) 
    }
    
    
    # fit an instrumental variable model
    if("instrument" %in% misspecify){
      instmod = ranger(Z ~ . , dat = cbind(w.instrument, Z = dat$Z)[slong!=split,])
      delta_1 = predict(instmod, data = w.instrument)$predictions
      delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
      delta_0 = (1-delta_1)
    }else{
      instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,])
      delta_1 = predict(instmod, data = x.instrument)$predictions
      delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
      delta_0 = (1-delta_1)
    }
    
    
    ## fit outcome models (time-varying)
    for(j in 2:k){
      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(time.list[j] == dat$obs.Y & dat$R == 0)
      dat$I = as.integer(time.list[j-1] < dat$obs.Y)
      
      
      if("censoring" %in% misspecify){
        w.miss_Z1_A1 = w.miss_Z1_A0 = w.miss_Z0_A1 = w.miss_Z0_A0 = w.miss
        w.miss_Z1_A1$Z = rep(1, nrow(w.miss)); w.miss_Z1_A1$A = rep(1, nrow(w.miss))
        w.miss_Z1_A0$Z = rep(1, nrow(w.miss)); w.miss_Z1_A0$A = rep(0, nrow(w.miss))
        w.miss_Z0_A1$Z = rep(0, nrow(w.miss)); w.miss_Z0_A1$A = rep(1, nrow(w.miss))
        w.miss_Z0_A0$Z = rep(0, nrow(w.miss)); w.miss_Z0_A0$A = rep(0, nrow(w.miss))
        w.miss$I = dat$I;
        w.miss_Z1_A1$I = w.miss_Z1_A0$I = w.miss_Z0_A1$I = w.miss_Z0_A0$I = rep(1, nrow(dat))
        
        missmod = ranger(C ~., data = cbind(w.miss, C = dat$C)[slong!=split,])
        omega_Z1_A1[j,slong==split] = predict(missmod, data = w.miss_Z1_A1)$predictions[slong==split]
        omega_Z1_A0[j,slong==split] = predict(missmod, data = w.miss_Z1_A0)$predictions[slong==split]
        omega_Z0_A1[j,slong==split] = predict(missmod, data = w.miss_Z0_A1)$predictions[slong==split]
        omega_Z0_A0[j,slong==split] = predict(missmod, data = w.miss_Z0_A0)$predictions[slong==split]
      }else{
        x.miss$I = dat$I;
        x.miss_Z1_A1$I = x.miss_Z1_A0$I = x.miss_Z0_A1$I = x.miss_Z0_A0$I = rep(1, nrow(dat))
        missmod = ranger(C ~., data = cbind(x.miss, C = dat$C)[slong!=split,])
        omega_Z1_A1[j,slong==split] = predict(missmod, data = x.miss_Z1_A1)$predictions[slong==split]
        omega_Z1_A0[j,slong==split] = predict(missmod, data = x.miss_Z1_A0)$predictions[slong==split]
        omega_Z0_A1[j,slong==split] = predict(missmod, data = x.miss_Z0_A1)$predictions[slong==split]
        omega_Z0_A0[j,slong==split] = predict(missmod, data = x.miss_Z0_A0)$predictions[slong==split]
      }
      
      
      G_tau_Z1_A1[j,slong==split] = G_tau_Z1_A1[j-1,slong==split]*(1-omega_Z1_A1[j,slong==split])
      G_tau_Z1_A0[j,slong==split] = G_tau_Z1_A0[j-1,slong==split]*(1-omega_Z1_A0[j,slong==split])
      G_tau_Z0_A1[j,slong==split] = G_tau_Z0_A1[j-1,slong==split]*(1-omega_Z0_A1[j,slong==split])
      G_tau_Z0_A0[j,slong==split] = G_tau_Z0_A0[j-1,slong==split]*(1-omega_Z0_A0[j,slong==split])
      
      G_tau_Z1_A1[j,slong==split] = ifelse(G_tau_Z1_A1[j,slong==split] < 0.05, 0.05, G_tau_Z1_A1[j,slong==split])
      G_tau_Z1_A0[j,slong==split] = ifelse(G_tau_Z1_A0[j,slong==split] < 0.05, 0.05, G_tau_Z1_A0[j,slong==split])
      G_tau_Z0_A1[j,slong==split] = ifelse(G_tau_Z0_A1[j,slong==split] < 0.05, 0.05, G_tau_Z0_A1[j,slong==split])
      G_tau_Z0_A0[j,slong==split] = ifelse(G_tau_Z0_A0[j,slong==split] < 0.05, 0.05, G_tau_Z0_A0[j,slong==split])
      
      
      
      
      phi1_1 = mean(((dat$obs.Y >= time.list[j])*dat$Z*dat$A/(G_tau_Z1_A1[j,]*pi_1*delta_1))[s == split])*
        mean((dat$A*dat$Z / (delta_1))[s == split]) + 
        mean(((dat$obs.Y >= time.list[j])*dat$Z*(1-dat$A)/(G_tau_Z1_A0[j,]*(1-pi_1)*delta_1))[s == split])*
        mean(((1-dat$A)*dat$Z / (delta_1))[s == split])
      
      phi1_0 =  mean(((dat$obs.Y >= time.list[j])*(1-dat$Z)*dat$A/(G_tau_Z0_A1[j,]*pi_0*delta_0))[s == split])*
        mean((dat$A*(1-dat$Z)/ (delta_0))[s == split] ) + 
        mean(((dat$obs.Y >= time.list[j])*(1-dat$Z)*(1-dat$A)/(G_tau_Z0_A0[j,]*(1-pi_0)*delta_0))[s == split])*
        mean(((1-dat$A)*(1-dat$Z) / (delta_0))[s == split])
      
      
      phi2_1 = dat$Z*dat$A/delta_1
      phi2_0 = (1-dat$Z)*dat$A/delta_0
      
      
      ifvals[s==split, j] = mean(phi1_1- phi1_0) / mean(phi2_1[s==split]- phi2_0[s==split]) 
    }
  } 
  
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
  }
  
  return(list(est = est.eff))
}


binary_survival_ranger_naive = function(dat, x.trt, x.out, x.miss, x.instrument = NULL, 
                                            w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  
  x.out_A1 = x.out_A0 = x.out
  x.out_A1$A = rep(1, nrow(x.out))
  x.out_A0$A = rep(0, nrow(x.out))
  
  x.miss_A1 = x.miss_A0 = x.miss
  x.miss_A1$A = rep(1, nrow(x.miss))
  x.miss_A0$A = rep(0, nrow(x.miss))
  
  n = nrow(dat)
  k = length(time.list)
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = est.var = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s
  tmp.var = matrix(NA, nsplits, k)
  
  S_tau_A1 = S_tau_A0 = 
    G_tau_A1 = G_tau_A0 = matrix(1, length(time.list), nrow(dat))
  h_A1 =  h_A0 = omega_A1 = omega_A0 = matrix(0, length(time.list), nrow(dat))
  D_A1 = D_A0 = matrix(0, length(time.list), nrow(dat))
  
  for(split in 1:nsplits){
    print(paste("split", split));
    flush.console() 
    
    # fit treatment model
    if("treatment" %in% misspecify){
      trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = w.trt)$predictions
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
    }else{
      trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
      pi_1 = predict(trtmod, data = x.trt)$predictions
      pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
    }
    
    
    for(j in 2:k){
      
      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(dat$obs.Y == time.list[j]  & dat$R == 0)
      dat$I = as.integer(time.list[j-1]  < dat$obs.Y)
      
      if("outcome" %in% misspecify){
        w.out_A1 = w.out_A0 =  w.out
        w.out_A1$A = rep(1, nrow(w.out))
        w.out_A0$A = rep(0, nrow(w.out))
        w.out$I = dat$I;
        w.out_A1$I = w.out_A0$I = rep(1, nrow(dat))
        
        outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
        h_A1[j,slong==split] = predict(outmod, data = w.out_A1)$predictions[slong==split]
        h_A0[j,slong==split] = predict(outmod, data = w.out_A0)$predictions[slong==split]
      }else{
        x.out$I = dat$I;
        x.out_A1$I = x.out_A0$I = rep(1, nrow(dat))
        outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
        h_A1[j,slong==split] = predict(outmod, data = x.out_A1)$predictions[slong==split]
        h_A0[j,slong==split] = predict(outmod, data = x.out_A0)$predictions[slong==split]
      }
      
      if("censoring" %in% misspecify){
        w.miss_A1 = w.miss_A0 = w.miss
        w.miss_A1$A = rep(1, nrow(w.miss))
        w.miss_A0$A = rep(0, nrow(w.miss))
        
        w.miss$I = dat$I;
        w.miss_A1$I = w.miss_A0$I = rep(1, nrow(dat))
        
        missmod = ranger(C ~., data = cbind(w.miss, C = dat$C)[slong!=split,])
        omega_A1[j,slong==split] = predict(missmod, data = w.miss_A1)$predictions[slong==split]
        omega_A0[j,slong==split] = predict(missmod, data = w.miss_A0)$predictions[slong==split]
        
      }else{
        x.miss$I = dat$I;
        x.miss_A1$I = x.miss_A0$I = rep(1, nrow(dat))
        missmod = ranger(C ~., data = cbind(x.miss, C = dat$C)[slong!=split,])
        omega_A1[j,slong==split] = predict(missmod, data = x.miss_A1)$predictions[slong==split]
        omega_A0[j,slong==split] = predict(missmod, data = x.miss_A0)$predictions[slong==split]
      }
      
      
      S_tau_A1[j,slong==split] = S_tau_A1[j-1,slong==split]*(1-h_A1[j,slong==split])
      S_tau_A0[j,slong==split] = S_tau_A0[j-1,slong==split]*(1-h_A0[j,slong==split])
      
      S_tau_A1[j,slong==split] = ifelse(S_tau_A1[j,slong==split] < 0.05, 0.05, S_tau_A1[j,slong==split])
      S_tau_A0[j,slong==split] = ifelse(S_tau_A0[j,slong==split] < 0.05, 0.05, S_tau_A0[j,slong==split])
      
      G_tau_A1[j,slong==split] = G_tau_A1[j-1,slong==split]*(1-omega_A1[j,slong==split])
      G_tau_A0[j,slong==split] = G_tau_A0[j-1,slong==split]*(1-omega_A0[j,slong==split])
      
      G_tau_A1[j,slong==split] = ifelse(G_tau_A1[j,slong==split] < 0.05, 0.05, G_tau_A1[j,slong==split])
      G_tau_A0[j,slong==split] = ifelse(G_tau_A0[j,slong==split] < 0.05, 0.05, G_tau_A0[j,slong==split])
      
      D_A1[j, ] = D_A0[j,] = rep(0, nrow(dat))
      for(t in 2:j){
        D_A1[j,slong==split] = D_A1[j,slong==split] +
          -(dat$A[slong==split]*as.integer(time.list[t-1] < dat$obs.Y[slong==split])*S_tau_A1[j,slong==split]*(as.integer(time.list[t] == dat$obs.Y[slong==split] &
                                                                                                                            dat$R[slong==split] == 1) - h_A1[t,slong==split])/(S_tau_A1[t,slong==split]*G_tau_A1[t-1,slong==split]*pi_1[slong==split]))
        D_A0[j,slong==split] = D_A0[j,slong==split] +
          -((1-dat$A[slong==split])*as.integer(time.list[t-1] < dat$obs.Y[slong==split])*S_tau_A0[j,slong==split]*(as.integer(time.list[t] == dat$obs.Y[slong==split] &
                                                                                                                                dat$R[slong==split] == 1) - h_A0[t,slong==split])/(S_tau_A0[t, slong==split]*G_tau_A0[t-1,slong==split]*(1-pi_1[slong==split])))
      }
      D_A1[j,slong==split] = D_A1[j,slong==split]+ S_tau_A1[j,slong==split]
      D_A0[j,slong==split] = D_A0[j,slong==split]+ S_tau_A0[j,slong==split]
      
      ifvals[slong==split, j] = mean(D_A1[j,slong==split] - D_A0[j,slong==split], na.rm = TRUE) 
      tmp.var[split,j] =  var(D_A1[j,slong==split] - D_A0[j,slong==split]) / nrow(dat[slong==split,])
    }
  }
  
  
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
    est.var[j] = sum(tmp.var[,j])/(nsplits)^2
  }
  
  return(list(est.eff = est.eff, est.var = est.var))
}