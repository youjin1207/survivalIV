binary_survival_ranger_if_colo = function(dat, x.trt, x.out, x.miss, x.instrument, 
   time.list, nsplits){
  ## x.trt : covariates for estimating the treatment assignment;
  ## x.out : covariates for estimating the outcome density; 
  ## x.miss : covariates for estimating the missing indicator;
  ## x. instrument : covariates for estimating the instrumental variable indicator

  x.trt_Z1 = x.trt_Z0 = x.trt
  x.trt_Z1$Z = rep(1, nrow(x.trt_Z1)); x.trt_Z0$Z = rep(0, nrow(x.trt_Z0))

  x.out_Z1_A1 = x.out_Z1_A0 = x.out_Z0_A1 = x.out_Z0_A0 = x.out
  x.out_Z1_A1$Z = rep(1, nrow(x.out)); x.out_Z1_A1$A = rep(1, nrow(x.out))
  x.out_Z1_A0$Z = rep(1, nrow(x.out)); x.out_Z1_A0$A = rep(0, nrow(x.out))
  #x.out_Z0_A1$Z = rep(0, nrow(x.out)); x.out_Z0_A1$A = rep(1, nrow(x.out))
  x.out_Z0_A0$Z = rep(0, nrow(x.out)); x.out_Z0_A0$A = rep(0, nrow(x.out))

  x.miss_Z1_A1 = x.miss_Z1_A0 = x.miss_Z0_A1 = x.miss_Z0_A0 = x.miss
  x.miss_Z1_A1$Z = rep(1, nrow(x.miss)); x.miss_Z1_A1$A = rep(1, nrow(x.miss))
  x.miss_Z1_A0$Z = rep(1, nrow(x.miss)); x.miss_Z1_A0$A = rep(0, nrow(x.miss))
  #x.miss_Z0_A1$Z = rep(0, nrow(x.miss)); x.miss_Z0_A1$A = rep(1, nrow(x.miss))
  x.miss_Z0_A0$Z = rep(0, nrow(x.miss)); x.miss_Z0_A0$A = rep(0, nrow(x.miss))

  n = nrow(dat)
  k = length(time.list)
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s

   for(split in 1:nsplits){
    print(paste("split", split));
    flush.console() 
    
    # fit a treatment model
    trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
    pi_1 = predict(trtmod, data = x.trt_Z1)$predictions
    pi_0 = rep(0, nrow(dat))
    pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
  
    # fit an instrumental variable model  
    instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,])
    delta_1 = predict(instmod, data = x.instrument)$predictions
    delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
    delta_0 = (1-delta_1)
    
    # fit a censoring model
    missmod = ranger(R ~., data = cbind(x.miss, R = dat$R)[slong!=split,])
    omega_Z1_A1 = predict(missmod, data = x.miss_Z1_A1)$predictions
    omega_Z1_A0 = predict(missmod, data = x.miss_Z1_A0)$predictions
    #omega_Z0_A1 = predict(missmod, data = x.miss_Z0_A1)$predictions
    omega_Z0_A0 = predict(missmod, data = x.miss_Z0_A0)$predictions
    
   
    ## fit outcome models (time-varying)
    for(j in 1:k){
      dat$Y = as.integer(time.list[j] < dat$obs.Y)
     
      # fit a survival probability model
      outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
      mu_Z1_A1 = predict(outmod, data = x.out_Z1_A1)$predictions
      mu_Z1_A0 = predict(outmod, data = x.out_Z1_A0)$predictions
      #mu_Z0_A1 = predict(outmod, data = x.out_Z0_A1)$predictions
      mu_Z0_A0 = predict(outmod, data = x.out_Z0_A0)$predictions
      
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
        phi2_0 = rep(0, nrow(dat))


        ## influence-function based estimator:
        ifvals[s==split, j] = mean(phi1_1[s==split] - phi1_0[s==split]) / 
                      mean(phi2_1[s==split] - phi2_0[s==split])            
      }
    }

    for(j in 1:k){
      est.eff[j] = mean(ifvals[,j])
    }


    return(list(est = est.eff))

  }



binary_survival_ranger_if_hazard_colo = function(dat, x.trt, x.out, x.miss, x.instrument, 
   time.list, nsplits){
  ## x.trt : covariates for estimating the treatment assignment;
  ## x.out : covariates for estimating the outcome density; 
  ## x.miss : covariates for estimating the missing indicator;
  ## x. instrument : covariates for estimating the instrumental variable indicator

  x.trt_Z1 = x.trt_Z0 = x.trt
  x.trt_Z1$Z = rep(1, nrow(x.trt_Z1)); x.trt_Z0$Z = rep(0, nrow(x.trt_Z0))

  x.out_Z1_A1 = x.out_Z1_A0 = x.out_Z0_A1 = x.out_Z0_A0 = x.out
  x.out_Z1_A1$Z = rep(1, nrow(x.out)); x.out_Z1_A1$A = rep(1, nrow(x.out))
  x.out_Z1_A0$Z = rep(1, nrow(x.out)); x.out_Z1_A0$A = rep(0, nrow(x.out))
  #x.out_Z0_A1$Z = rep(0, nrow(x.out)); x.out_Z0_A1$A = rep(1, nrow(x.out))
  x.out_Z0_A0$Z = rep(0, nrow(x.out)); x.out_Z0_A0$A = rep(0, nrow(x.out))

  x.miss_Z1_A1 = x.miss_Z1_A0 = x.miss_Z0_A1 = x.miss_Z0_A0 = x.miss
  x.miss_Z1_A1$Z = rep(1, nrow(x.miss)); x.miss_Z1_A1$A = rep(1, nrow(x.miss))
  x.miss_Z1_A0$Z = rep(1, nrow(x.miss)); x.miss_Z1_A0$A = rep(0, nrow(x.miss))
  #x.miss_Z0_A1$Z = rep(0, nrow(x.miss)); x.miss_Z0_A1$A = rep(1, nrow(x.miss))
  x.miss_Z0_A0$Z = rep(0, nrow(x.miss)); x.miss_Z0_A0$A = rep(0, nrow(x.miss))

  n = nrow(dat)
  k = length(time.list)
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s

   for(split in 1:nsplits){
    print(paste("split", split));
    flush.console() 
    
    # fit a treatment model
    trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
    pi_1 = predict(trtmod, data = x.trt_Z1)$predictions
    pi_0 = rep(0, nrow(dat))
    pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
  
    # fit an instrumental variable model  
    instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,])
    delta_1 = predict(instmod, data = x.instrument)$predictions
    delta_1 = ifelse(delta_1 < 0.05, 0.05, ifelse(delta_1 > 0.95, 0.95, delta_1))
    delta_0 = (1-delta_1)

    phi2_1 = pi_1 + 
          dat$Z*(dat$A - pi_1)/delta_1
  
  
    S_tau_Z1_A1 = S_tau_Z1_A0 = S_tau_Z0_A1 = S_tau_Z0_A0 = 
    G_tau_Z1_A1 = G_tau_Z1_A0 = G_tau_Z0_A1 = G_tau_Z0_A0 = matrix(1, length(time.list), nrow(dat))
    h_Z1_A1 =  h_Z1_A0 =  h_Z0_A1 = h_Z0_A0 = 
      omega_Z1_A1 = omega_Z1_A0 = omega_Z0_A1 = omega_Z0_A0 = matrix(0, length(time.list), nrow(dat))
    D_Z1_A1 = D_Z1_A0 = D_Z0_A1 = D_Z0_A0 = matrix(0, length(time.list), nrow(dat))
   
    ## fit outcome models (time-varying)
    for(j in 2:k){
    
      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(time.list[j] == dat$obs.Y & dat$R == 0)
      dat$I = as.integer(time.list[j-1] < dat$obs.Y)
     
      # fit a discrete hazard model
      x.out$I = dat$I;
      x.out_Z1_A1$I = x.out_Z1_A0$I = x.out_Z0_A0$I = rep(1, nrow(dat))
      outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
      h_Z1_A1[j,] = predict(outmod, data = x.out_Z1_A1)$predictions
      h_Z1_A0[j,] = predict(outmod, data = x.out_Z1_A0)$predictions
      h_Z0_A0[j,] = predict(outmod, data = x.out_Z0_A0)$predictions

      # fit a censoring model
      x.miss$I = dat$I;
      x.miss_Z1_A1$I = x.miss_Z1_A0$I = x.miss_Z0_A0$I = rep(1, nrow(dat))
      missmod = ranger(C ~., data = cbind(x.miss, C = dat$C)[slong!=split,])
      omega_Z1_A1[j,] = predict(missmod, data = x.miss_Z1_A1)$predictions
      omega_Z1_A0[j,] = predict(missmod, data = x.miss_Z1_A0)$predictions
      omega_Z0_A0[j,] = predict(missmod, data = x.miss_Z0_A0)$predictions

      
      S_tau_Z1_A1[j,] = S_tau_Z1_A1[j-1,]*(1-h_Z1_A1[j,])
      S_tau_Z1_A0[j,] = S_tau_Z1_A0[j-1,]*(1-h_Z1_A0[j,])
      S_tau_Z0_A0[j,] = S_tau_Z0_A0[j-1,]*(1-h_Z0_A0[j,])
     
    
      G_tau_Z1_A1[j,] = G_tau_Z1_A1[j-1,]*(1-omega_Z1_A1[j,])
      G_tau_Z1_A0[j,] = G_tau_Z1_A0[j-1,]*(1-omega_Z1_A0[j,])
      G_tau_Z0_A0[j,] = G_tau_Z0_A0[j-1,]*(1-omega_Z0_A0[j,])
      
      
      for(t in 2:j){
            D_Z1_A1[j,] = D_Z1_A1[j,] +
              - pi_1*(dat$Z*dat$A*as.integer(dat$obs.Y > time.list[t-1])*S_tau_Z1_A1[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                   dat$R == 1) - h_Z1_A1[t,])/(S_tau_Z1_A1[t, ]*G_tau_Z1_A1[t-1,]*pi_1*delta_1))
            D_Z1_A0[j,] = D_Z1_A0[j,] +
              - (1-pi_1)*(dat$Z*(1-dat$A)*as.integer(dat$obs.Y > time.list[t-1])*S_tau_Z1_A0[j,]*(as.integer(dat$obs.Y == time.list[t] &
                                                   dat$R == 1) - h_Z1_A0[t,])/(S_tau_Z1_A0[t, ]*G_tau_Z1_A0[t-1,]*(1-pi_1)*delta_1))
        
            D_Z0_A0[j,] = D_Z0_A0[j,] +
              - (1-pi_0)*((1-dat$Z)*(1-dat$A)*as.integer(dat$obs.Y > time.list[t-1])*S_tau_Z0_A0[j,]*(as.integer(dat$obs.Y == time.list[t] &
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
    }

    for(j in 1:k){
      est.eff[j] = mean(ifvals[,j])
    }


    return(list(est = est.eff))

  }

binary_survival_ranger_naive_colo = function(dat, x.trt, x.out, x.miss, x.instrument = NULL, 
   time.list, nsplits){
  ## x.trt : covariates for estimating the treatment assignment;
  ## x.out : covariates for estimating the outcome density; 
  ## x.miss : covariates for estimating the missing indicator;
  ## x. instrument : covariates for estimating the instrumental variable indicator

  x.out_A1 = x.out_A0 = x.out
  x.out_A1$A = rep(1, nrow(x.out))
  x.out_A0$A = rep(0, nrow(x.out))
 
  x.miss_A1 = x.miss_A0 = x.miss
  x.miss_A1$A = rep(1, nrow(x.miss))
  x.miss_A0$A = rep(0, nrow(x.miss))

  n = nrow(dat)
  k = length(time.list)
  ifvals = matrix(0, nrow = n, ncol = k)
  est.eff = rep(NA, k)
  s = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  slong = s

   for(split in 1:nsplits){
    print(paste("split", split));
    flush.console() 
    
    # fit a treatment model
    trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
    pi_1 = predict(trtmod, data = x.trt)$predictions
    pi_0 = rep(0, nrow(dat))
    pi_1 = ifelse(pi_1 < 0.05, 0.05, ifelse(pi_1 > 0.95, 0.95, pi_1)) 
  

    S_tau_A1 = S_tau_A0 = 
    G_tau_A1 = G_tau_A0 = matrix(1, length(time.list), nrow(dat))
    h_A1 =  h_A0 = omega_A1 = omega_A0 = matrix(0, length(time.list), nrow(dat))
    D_A1 = D_A0 = matrix(0, length(time.list), nrow(dat))
   
    ## fit outcome models (time-varying)
    for(j in 2:k){
    
      dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
      dat$C = as.integer(time.list[j] == dat$obs.Y & dat$R == 0)
      dat$I = as.integer(time.list[j-1] < dat$obs.Y)
     
      # fit a survival probability model
      x.out$I = dat$I;
      x.out_A1$I = x.out_A0$I = rep(1, nrow(dat))
      outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
      h_A1[j,] = predict(outmod, data = x.out_A1)$predictions
      h_A0[j,] = predict(outmod, data = x.out_A0)$predictions
    
      # fit a censoring model
      x.miss$I = dat$I;
      x.miss_A1$I = x.miss_A0$I =  rep(1, nrow(dat))
      missmod = ranger(C ~., data = cbind(x.miss, C = dat$C)[slong!=split,])
      omega_A1[j,] = predict(missmod, data = x.miss_A1)$predictions
      omega_A0[j,] = predict(missmod, data = x.miss_A0)$predictions
  
      if(j > 1){
        S_tau_A1[j,] = S_tau_A1[j-1,]*(1-h_A1[j,])
        S_tau_A0[j,] = S_tau_A0[j-1,]*(1-h_A0[j,])
      }else{
        S_tau_A1[j,] = 1*(1-h_A1[j,])
        S_tau_A0[j,] = 1*(1-h_A0[j,])
      }

      if(j>1){
        G_tau_A1[j,] = G_tau_A1[j-1,]*(1-omega_A1[j,])
        G_tau_A0[j,] = G_tau_A0[j-1,]*(1-omega_A0[j,])
      }else{
        G_tau_A1[j,] = 1*(1-omega_A1[j,])
        G_tau_A0[j,] = 1*(1-omega_A0[j,])
      }
      
     
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
  }

    for(j in 1:k){
      est.eff[j] = mean(ifvals[,j])
    }

    return(list(est = est.eff))

  }  