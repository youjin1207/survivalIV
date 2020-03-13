continuous_survival_ipw = function(dat, kappa, x.trt, x.out, x.miss, x.instrument, 
  w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  ## x.trt : covariates for estimating the treatment assignment;
  ## x.out : covariates for estimating the outcome density; 
  ## x.miss : covariates for estimating the missing indicator;
  ## x. instrument : covariates for estimating the instrumental variable indicator
  x.trt_Z = x.trt_Z_minus = x.trt_Z_plus = x.trt
  x.trt_Z_minus$Z = dat$Z - kappa; x.trt_Z_plus$Z = dat$Z + kappa 

  x.out_Z_A1 = x.out_Z_minus_A1 = x.out_Z_plus_A1 = x.out 
  x.out_Z_A0 = x.out_Z_minus_A0 = x.out_Z_plus_A0 = x.out 
  x.out_Z_A1$A = rep(1, nrow(x.out));  x.out_Z_A0$A = rep(0, nrow(x.out))
  x.out_Z_minus_A1$Z = dat$Z-kappa; x.out_Z_minus_A1$A = rep(1, nrow(x.out))
  x.out_Z_plus_A1$Z = dat$Z+kappa; x.out_Z_plus_A1$A = rep(1, nrow(x.out))
  x.out_Z_minus_A0$Z = dat$Z-kappa; x.out_Z_minus_A0$A = rep(0, nrow(x.out))
  x.out_Z_plus_A0$Z = dat$Z+kappa; x.out_Z_plus_A0$A = rep(0, nrow(x.out))

  x.miss_Z_A1 = x.miss_Z_minus_A1 = x.miss_Z_plus_A1 = x.out 
  x.miss_Z_A0 = x.miss_Z_minus_A0 = x.miss_Z_plus_A0 = x.out 
  x.miss_Z_A1$A = rep(1, nrow(x.out));  x.miss_Z_A0$A = rep(0, nrow(x.out))
  x.miss_Z_minus_A1$Z = dat$Z-kappa; x.miss_Z_minus_A1$A = rep(1, nrow(x.out))
  x.miss_Z_plus_A1$Z = dat$Z+kappa; x.miss_Z_plus_A1$A = rep(1, nrow(x.out))
  x.miss_Z_minus_A0$Z = dat$Z-kappa; x.miss_Z_minus_A0$A = rep(0, nrow(x.out))
  x.miss_Z_plus_A0$Z = dat$Z+kappa; x.miss_Z_plus_A0$A = rep(0, nrow(x.out))



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
        w.trt_Z = w.trt_Z_minus = w.trt_Z_plus = w.trt
        w.trt_Z_minus$Z = dat$Z - kappa; w.trt_Z_plus$Z = dat$Z + kappa 
        trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
        pi_Z = predict(trtmod, data = w.trt_Z)$predictions
        pi_Z_minus = predict(trtmod, data = w.trt_Z_minus)$predictions
        pi_Z_plus = predict(trtmod, data = w.trt_Z_plus)$predictions
      }else{
        trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
        pi_Z = predict(trtmod, data = x.trt_Z)$predictions
        pi_Z_minus = predict(trtmod, data = x.trt_Z_minus)$predictions
        pi_Z_plus = predict(trtmod, data = x.trt_Z_plus)$predictions
      }
      pi_Z = ifelse(pi_Z < 0.01, 0.01, ifelse(pi_Z > 0.99, 0.99, pi_Z))
      pi_Z_plus = ifelse(pi_Z_plus < 0.01, 0.01, ifelse(pi_Z_plus > 0.99, 0.99, pi_Z_plus))
      pi_Z_minus = ifelse(pi_Z_minus < 0.01, 0.01, ifelse(pi_Z_minus > 0.99, 0.99, pi_Z_minus))
    

      if("instrument" %in% misspecify){
        instmod = ranger(Z ~ . , dat = cbind(w.instrument, Z = dat$Z)[slong!=split,],
          probability = TRUE) # this takes long time
         Zvalues = dat$Z[slong!=split]
         Z_predict = predict(instmod, data = cbind(w.instrument))$predictions     
       }else{
           instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,],
          probability = TRUE) # this takes long time   
          Zvalues = dat$Z[slong!=split]
          Z_predict = predict(instmod, data = cbind(x.instrument))$predictions
       }
      

        delta_Z = delta_Z_plus = delta_Z_minus = c()
       for(i in 1:nrow(dat)){
        tmp.smooth = density(Zvalues, weights = Z_predict[i,])
        delta_Z[i] = tmp.smooth$y[which.min( abs(tmp.smooth$x - dat$Z[i]))]
        delta_Z_plus[i] = tmp.smooth$y[which.min( abs(tmp.smooth$x - (dat$Z[i] + kappa) ))]
        delta_Z_minus[i] = tmp.smooth$y[which.min( abs(tmp.smooth$x - (dat$Z[i] - kappa) ))]
       }
       
      delta_Z = ifelse(delta_Z < 0.01, 0.01, ifelse(delta_Z > 0.99, 0.99, delta_Z))
      delta_Z_plus = ifelse(delta_Z_plus < 0.01, 0.01, ifelse(delta_Z_plus > 0.99, 0.99, delta_Z_plus))
      delta_Z_minus = ifelse(delta_Z_minus < 0.01, 0.01, ifelse(delta_Z_minus > 0.99, 0.99, delta_Z_minus))

      if("censoring" %in% misspecify){

        w.miss_Z_A1 = w.miss_Z_minus_A1 = w.miss_Z_plus_A1 = w.out 
        w.miss_Z_A0 = w.miss_Z_minus_A0 = w.miss_Z_plus_A0 = w.out 
        w.miss_Z_A1$A = rep(1, nrow(w.out));  w.miss_Z_A0$A = rep(0, nrow(w.out))
        w.miss_Z_minus_A1$Z = dat$Z-kappa; w.miss_Z_minus_A1$A = rep(1, nrow(w.out))
        w.miss_Z_plus_A1$Z = dat$Z+kappa; w.miss_Z_plus_A1$A = rep(1, nrow(w.out))
        w.miss_Z_minus_A0$Z = dat$Z-kappa; w.miss_Z_minus_A0$A = rep(0, nrow(w.out))
        w.miss_Z_plus_A0$Z = dat$Z+kappa; w.miss_Z_plus_A0$A = rep(0, nrow(w.out))

         missmod = ranger(R ~., data = cbind(w.miss, R = dat$R)[slong!=split,])
          omega_Z_A1 = predict(missmod, data = w.miss_Z_A1)$predictions
          omega_Z_minus_A1 = predict(missmod, data = w.miss_Z_minus_A1)$predictions
          omega_Z_plus_A1 = predict(missmod, data = w.miss_Z_plus_A1)$predictions
          omega_Z_A0 = predict(missmod, data = w.miss_Z_A0)$predictions
          omega_Z_minus_A0 = predict(missmod, data = w.miss_Z_minus_A0)$predictions
          omega_Z_plus_A0 = predict(missmod, data = w.miss_Z_plus_A0)$predictions
        }else{
           missmod = ranger(R ~., data = cbind(x.miss, R = dat$R)[slong!=split,])
        omega_Z_A1 = predict(missmod, data = x.miss_Z_A1)$predictions
        omega_Z_minus_A1 = predict(missmod, data = x.miss_Z_minus_A1)$predictions
        omega_Z_plus_A1 = predict(missmod, data = x.miss_Z_plus_A1)$predictions
        omega_Z_A0 = predict(missmod, data = x.miss_Z_A0)$predictions
        omega_Z_minus_A0 = predict(missmod, data = x.miss_Z_minus_A0)$predictions
        omega_Z_plus_A0 = predict(missmod, data = x.miss_Z_plus_A0)$predictions
        }
       
     
      ## fit outcome models (time-varying)
      for(j in 1:k){
        dat$Y = as.integer(time.list[j] < dat$obs.Y)
       

        if("outcome" %in% misspecify){
           w.out_Z_A1 = w.out_Z_minus_A1 = w.out_Z_plus_A1 = w.out 
          w.out_Z_A0 = w.out_Z_minus_A0 = w.out_Z_plus_A0 = w.out 
          w.out_Z_A1$A = rep(1, nrow(w.out));  w.out_Z_A0$A = rep(0, nrow(w.out))
          w.out_Z_minus_A1$Z = dat$Z-kappa; w.out_Z_minus_A1$A = rep(1, nrow(w.out))
          w.out_Z_plus_A1$Z = dat$Z+kappa; w.out_Z_plus_A1$A = rep(1, nrow(w.out))
          w.out_Z_minus_A0$Z = dat$Z-kappa; w.out_Z_minus_A0$A = rep(0, nrow(w.out))
          w.out_Z_plus_A0$Z = dat$Z+kappa; w.out_Z_plus_A0$A = rep(0, nrow(w.out))
            outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
            mu_Z_A1 = predict(outmod, data = w.out_Z_A1)$predictions
            mu_Z_minus_A1 = predict(outmod, data = w.out_Z_minus_A1)$predictions
            mu_Z_plus_A1 = predict(outmod, data = w.out_Z_plus_A1)$predictions
            mu_Z_A0 = predict(outmod, data = w.out_Z_A0)$predictions
            mu_Z_minus_A0 = predict(outmod, data = w.out_Z_minus_A0)$predictions
            mu_Z_plus_A0 = predict(outmod, data = w.out_Z_plus_A0)$predictions
          }else{       
            outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
            mu_Z_A1 = predict(outmod, data = x.out_Z_A1)$predictions
            mu_Z_minus_A1 = predict(outmod, data = x.out_Z_minus_A1)$predictions
            mu_Z_plus_A1 = predict(outmod, data = x.out_Z_plus_A1)$predictions
            mu_Z_A0 = predict(outmod, data = x.out_Z_A0)$predictions
            mu_Z_minus_A0 = predict(outmod, data = x.out_Z_minus_A0)$predictions
            mu_Z_plus_A0 = predict(outmod, data = x.out_Z_plus_A0)$predictions
          }
      
      
        phi1_1 = mean((dat$Y*dat$R*dat$A*delta_Z_minus/(omega_Z_A1*pi_Z*delta_Z))[s==split])*
              mean((dat$A*delta_Z_minus/(delta_Z))[s==split]) +
              mean((dat$Y*dat$R*(1-dat$A)*delta_Z_minus/(omega_Z_A0*(1-pi_Z)*delta_Z))[s==split])*
              mean(((1-dat$A)*delta_Z_minus / (delta_Z))[s==split])
          phi1_0 =  mean((dat$Y*dat$R*dat$A*delta_Z_plus/(omega_Z_A1*pi_Z*delta_Z))[s==split])*
              mean((dat$A*delta_Z_plus/ (delta_Z))[s==split] ) +
              mean((dat$Y*dat$R*(1-dat$A)*delta_Z_plus/(omega_Z_A0*(1-pi_Z)*delta_Z))[s==split])*
              mean(((1-dat$A)*delta_Z_plus/ (delta_Z))[s==split])
          phi2_1 = dat$A*delta_Z_minus/delta_Z
          phi2_0 = dat$A*delta_Z_plus/delta_Z

            ## influence-function based estimator:
          ifvals[s==split, j] = mean(phi1_1 - phi1_0) /
                      mean(phi2_1[s==split] - phi2_0[s==split])
      }
    }


  
  # compute estimator (using)
  for(j in 1:k){
    est.eff[j] = mean(ifvals[,j]) # for each delta vaule
  }
    
  return(list(est = est.eff))
              
}

continuous_survival_plugin = function(dat, kappa, x.trt, x.out, x.miss, x.instrument, 
  w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  ## x.trt : covariates for estimating the treatment assignment;
  ## x.out : covariates for estimating the outcome density; 
  ## x.miss : covariates for estimating the missing indicator;
  ## x. instrument : covariates for estimating the instrumental variable indicator
 x.trt_Z = x.trt_Z_minus = x.trt_Z_plus = x.trt
  x.trt_Z_minus$Z = dat$Z - kappa; x.trt_Z_plus$Z = dat$Z + kappa 

  x.out_Z_A1 = x.out_Z_minus_A1 = x.out_Z_plus_A1 = x.out 
  x.out_Z_A0 = x.out_Z_minus_A0 = x.out_Z_plus_A0 = x.out 
  x.out_Z_A1$A = rep(1, nrow(x.out));  x.out_Z_A0$A = rep(0, nrow(x.out))
  x.out_Z_minus_A1$Z = dat$Z-kappa; x.out_Z_minus_A1$A = rep(1, nrow(x.out))
  x.out_Z_plus_A1$Z = dat$Z+kappa; x.out_Z_plus_A1$A = rep(1, nrow(x.out))
  x.out_Z_minus_A0$Z = dat$Z-kappa; x.out_Z_minus_A0$A = rep(0, nrow(x.out))
  x.out_Z_plus_A0$Z = dat$Z+kappa; x.out_Z_plus_A0$A = rep(0, nrow(x.out))

  x.miss_Z_A1 = x.miss_Z_minus_A1 = x.miss_Z_plus_A1 = x.out 
  x.miss_Z_A0 = x.miss_Z_minus_A0 = x.miss_Z_plus_A0 = x.out 
  x.miss_Z_A1$A = rep(1, nrow(x.out));  x.miss_Z_A0$A = rep(0, nrow(x.out))
  x.miss_Z_minus_A1$Z = dat$Z-kappa; x.miss_Z_minus_A1$A = rep(1, nrow(x.out))
  x.miss_Z_plus_A1$Z = dat$Z+kappa; x.miss_Z_plus_A1$A = rep(1, nrow(x.out))
  x.miss_Z_minus_A0$Z = dat$Z-kappa; x.miss_Z_minus_A0$A = rep(0, nrow(x.out))
  x.miss_Z_plus_A0$Z = dat$Z+kappa; x.miss_Z_plus_A0$A = rep(0, nrow(x.out))




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
        w.trt_Z = w.trt_Z_minus = w.trt_Z_plus = w.trt
        w.trt_Z_minus$Z = dat$Z - kappa; w.trt_Z_plus$Z = dat$Z + kappa 
        trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
        pi_Z = predict(trtmod, data = w.trt_Z)$predictions
        pi_Z_minus = predict(trtmod, data = w.trt_Z_minus)$predictions
        pi_Z_plus = predict(trtmod, data = w.trt_Z_plus)$predictions
      }else{
        trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
        pi_Z = predict(trtmod, data = x.trt_Z)$predictions
        pi_Z_minus = predict(trtmod, data = x.trt_Z_minus)$predictions
        pi_Z_plus = predict(trtmod, data = x.trt_Z_plus)$predictions
      }
       pi_Z = ifelse(pi_Z < 0.01, 0.01, ifelse(pi_Z > 0.99, 0.99, pi_Z))
      pi_Z_plus = ifelse(pi_Z_plus < 0.01, 0.01, ifelse(pi_Z_plus > 0.99, 0.99, pi_Z_plus))
      pi_Z_minus = ifelse(pi_Z_minus < 0.01, 0.01, ifelse(pi_Z_minus > 0.99, 0.99, pi_Z_minus))
    

      if("censoring" %in% misspecify){


        w.miss_Z_A1 = w.miss_Z_minus_A1 = w.miss_Z_plus_A1 = w.out 
        w.miss_Z_A0 = w.miss_Z_minus_A0 = w.miss_Z_plus_A0 = w.out 
        w.miss_Z_A1$A = rep(1, nrow(w.out));  w.miss_Z_A0$A = rep(0, nrow(w.out))
        w.miss_Z_minus_A1$Z = dat$Z-kappa; w.miss_Z_minus_A1$A = rep(1, nrow(w.out))
        w.miss_Z_plus_A1$Z = dat$Z+kappa; w.miss_Z_plus_A1$A = rep(1, nrow(w.out))
        w.miss_Z_minus_A0$Z = dat$Z-kappa; w.miss_Z_minus_A0$A = rep(0, nrow(w.out))
        w.miss_Z_plus_A0$Z = dat$Z+kappa; w.miss_Z_plus_A0$A = rep(0, nrow(w.out))

         missmod = ranger(R ~., data = cbind(w.miss, R = dat$R)[slong!=split,])
          omega_Z_A1 = predict(missmod, data = w.miss_Z_A1)$predictions
          omega_Z_minus_A1 = predict(missmod, data = w.miss_Z_minus_A1)$predictions
          omega_Z_plus_A1 = predict(missmod, data = w.miss_Z_plus_A1)$predictions
          omega_Z_A0 = predict(missmod, data = w.miss_Z_A0)$predictions
          omega_Z_minus_A0 = predict(missmod, data = w.miss_Z_minus_A0)$predictions
          omega_Z_plus_A0 = predict(missmod, data = w.miss_Z_plus_A0)$predictions
        }else{
           missmod = ranger(R ~., data = cbind(x.miss, R = dat$R)[slong!=split,])
        omega_Z_A1 = predict(missmod, data = x.miss_Z_A1)$predictions
        omega_Z_minus_A1 = predict(missmod, data = x.miss_Z_minus_A1)$predictions
        omega_Z_plus_A1 = predict(missmod, data = x.miss_Z_plus_A1)$predictions
        omega_Z_A0 = predict(missmod, data = x.miss_Z_A0)$predictions
        omega_Z_minus_A0 = predict(missmod, data = x.miss_Z_minus_A0)$predictions
        omega_Z_plus_A0 = predict(missmod, data = x.miss_Z_plus_A0)$predictions
        }
       
     
      ## fit outcome models (time-varying)
      for(j in 1:k){
        dat$Y = as.integer(time.list[j] < dat$obs.Y)
      

        if("outcome" %in% misspecify){
          w.out_Z_A1 = w.out_Z_minus_A1 = w.out_Z_plus_A1 = w.out 
          w.out_Z_A0 = w.out_Z_minus_A0 = w.out_Z_plus_A0 = w.out 
          w.out_Z_A1$A = rep(1, nrow(w.out));  w.out_Z_A0$A = rep(0, nrow(w.out))
          w.out_Z_minus_A1$Z = dat$Z-kappa; w.out_Z_minus_A1$A = rep(1, nrow(w.out))
          w.out_Z_plus_A1$Z = dat$Z+kappa; w.out_Z_plus_A1$A = rep(1, nrow(w.out))
          w.out_Z_minus_A0$Z = dat$Z-kappa; w.out_Z_minus_A0$A = rep(0, nrow(w.out))
          w.out_Z_plus_A0$Z = dat$Z+kappa; w.out_Z_plus_A0$A = rep(0, nrow(w.out))
            outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
            mu_Z_A1 = predict(outmod, data = w.out_Z_A1)$predictions
            mu_Z_minus_A1 = predict(outmod, data = w.out_Z_minus_A1)$predictions
            mu_Z_plus_A1 = predict(outmod, data = w.out_Z_plus_A1)$predictions
            mu_Z_A0 = predict(outmod, data = w.out_Z_A0)$predictions
            mu_Z_minus_A0 = predict(outmod, data = w.out_Z_minus_A0)$predictions
            mu_Z_plus_A0 = predict(outmod, data = w.out_Z_plus_A0)$predictions
          }else{       
            outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
            mu_Z_A1 = predict(outmod, data = x.out_Z_A1)$predictions
            mu_Z_minus_A1 = predict(outmod, data = x.out_Z_minus_A1)$predictions
            mu_Z_plus_A1 = predict(outmod, data = x.out_Z_plus_A1)$predictions
            mu_Z_A0 = predict(outmod, data = x.out_Z_A0)$predictions
            mu_Z_minus_A0 = predict(outmod, data = x.out_Z_minus_A0)$predictions
            mu_Z_plus_A0 = predict(outmod, data = x.out_Z_plus_A0)$predictions
          }
      
      
        phi1_1 = mu_Z_plus_A1*pi_Z_plus  + mu_Z_plus_A0*(1-pi_Z_plus) 


    
        phi1_0 = mu_Z_minus_A1*pi_Z_minus  + mu_Z_minus_A0*(1-pi_Z_minus) 
    
    
        phi2_1 = pi_Z_plus 
        phi2_0 = pi_Z_minus 


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
  
  #return(list(est = est.eff, sigma = sigma,
  #            ll1 = eff.ll, ul1 = eff.ul, calpha = calpha, ll2 = eff.ll2, ul2 = eff.ul2))
  
}

continuous_survival_if = function(dat, kappa, x.trt, x.out, x.miss, x.instrument, 
  w.trt = NULL, w.out = NULL, w.miss = NULL, w.instrument = NULL, time.list, nsplits, misspecify = NULL){
  ## x.trt : covariates for estimating the treatment assignment;
  ## x.out : covariates for estimating the outcome density; 
  ## x.miss : covariates for estimating the missing indicator;
  ## x. instrument : covariates for estimating the instrumental variable indicator
  x.trt_Z = x.trt_Z_minus = x.trt_Z_plus = x.trt
  x.trt_Z_minus$Z = dat$Z - kappa; x.trt_Z_plus$Z = dat$Z + kappa 

  x.out_Z_A1 = x.out_Z_minus_A1 = x.out_Z_plus_A1 = x.out 
  x.out_Z_A0 = x.out_Z_minus_A0 = x.out_Z_plus_A0 = x.out 
  x.out_Z_A1$A = rep(1, nrow(x.out));  x.out_Z_A0$A = rep(0, nrow(x.out))
  x.out_Z_minus_A1$Z = dat$Z-kappa; x.out_Z_minus_A1$A = rep(1, nrow(x.out))
  x.out_Z_plus_A1$Z = dat$Z+kappa; x.out_Z_plus_A1$A = rep(1, nrow(x.out))
  x.out_Z_minus_A0$Z = dat$Z-kappa; x.out_Z_minus_A0$A = rep(0, nrow(x.out))
  x.out_Z_plus_A0$Z = dat$Z+kappa; x.out_Z_plus_A0$A = rep(0, nrow(x.out))

  x.miss_Z_A1 = x.miss_Z_minus_A1 = x.miss_Z_plus_A1 = x.out 
  x.miss_Z_A0 = x.miss_Z_minus_A0 = x.miss_Z_plus_A0 = x.out 
  x.miss_Z_A1$A = rep(1, nrow(x.out));  x.miss_Z_A0$A = rep(0, nrow(x.out))
  x.miss_Z_minus_A1$Z = dat$Z-kappa; x.miss_Z_minus_A1$A = rep(1, nrow(x.out))
  x.miss_Z_plus_A1$Z = dat$Z+kappa; x.miss_Z_plus_A1$A = rep(1, nrow(x.out))
  x.miss_Z_minus_A0$Z = dat$Z-kappa; x.miss_Z_minus_A0$A = rep(0, nrow(x.out))
  x.miss_Z_plus_A0$Z = dat$Z+kappa; x.miss_Z_plus_A0$A = rep(0, nrow(x.out))



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
        w.trt_Z = w.trt_Z_minus = w.trt_Z_plus = w.trt
        w.trt_Z_minus$Z = dat$Z - kappa; w.trt_Z_plus$Z = dat$Z + kappa 
        trtmod = ranger(A ~ . , dat = cbind(w.trt, A = dat$A)[slong!=split,])
        pi_Z = predict(trtmod, data = w.trt_Z)$predictions
        pi_Z_minus = predict(trtmod, data = w.trt_Z_minus)$predictions
        pi_Z_plus = predict(trtmod, data = w.trt_Z_plus)$predictions
      }else{
        trtmod = ranger(A ~ . , dat = cbind(x.trt, A = dat$A)[slong!=split,])
        pi_Z = predict(trtmod, data = x.trt_Z)$predictions
        pi_Z_minus = predict(trtmod, data = x.trt_Z_minus)$predictions
        pi_Z_plus = predict(trtmod, data = x.trt_Z_plus)$predictions
      }
       pi_Z = ifelse(pi_Z < 0.01, 0.01, ifelse(pi_Z > 0.99, 0.99, pi_Z))
      pi_Z_plus = ifelse(pi_Z_plus < 0.01, 0.01, ifelse(pi_Z_plus > 0.99, 0.99, pi_Z_plus))
      pi_Z_minus = ifelse(pi_Z_minus < 0.01, 0.01, ifelse(pi_Z_minus > 0.99, 0.99, pi_Z_minus))
    

      if("instrument" %in% misspecify){
        instmod = ranger(Z ~ . , dat = cbind(w.instrument, Z = dat$Z)[slong!=split,],
          probability = TRUE) # this takes long time
         Zvalues = dat$Z[slong!=split]
         Z_predict = predict(instmod, data = cbind(w.instrument))$predictions     
       }else{
           instmod = ranger(Z ~ . , dat = cbind(x.instrument, Z = dat$Z)[slong!=split,],
          probability = TRUE) # this takes long time   
          Zvalues = dat$Z[slong!=split]
          Z_predict = predict(instmod, data = cbind(x.instrument))$predictions
       }
      

        delta_Z = delta_Z_plus = delta_Z_minus = c()
       for(i in 1:nrow(dat)){
        tmp.smooth = density(Zvalues, weights = Z_predict[i,])
        delta_Z[i] = tmp.smooth$y[which.min( abs(tmp.smooth$x - dat$Z[i]))]
        delta_Z_plus[i] = tmp.smooth$y[which.min( abs(tmp.smooth$x - (dat$Z[i] + kappa) ))]
        delta_Z_minus[i] = tmp.smooth$y[which.min( abs(tmp.smooth$x - (dat$Z[i] - kappa) ))]
       }
       
      delta_Z = ifelse(delta_Z < 0.01, 0.01, ifelse(delta_Z > 0.99, 0.99, delta_Z))
      delta_Z_plus = ifelse(delta_Z_plus < 0.01, 0.01, ifelse(delta_Z_plus > 0.99, 0.99, delta_Z_plus))
      delta_Z_minus = ifelse(delta_Z_minus < 0.01, 0.01, ifelse(delta_Z_minus > 0.99, 0.99, delta_Z_minus))

      if("censoring" %in% misspecify){


          w.miss_Z_A1 = w.miss_Z_minus_A1 = w.miss_Z_plus_A1 = w.out 
          w.miss_Z_A0 = w.miss_Z_minus_A0 = w.miss_Z_plus_A0 = w.out 
          w.miss_Z_A1$A = rep(1, nrow(w.out));  w.miss_Z_A0$A = rep(0, nrow(w.out))
          w.miss_Z_minus_A1$Z = dat$Z-kappa; w.miss_Z_minus_A1$A = rep(1, nrow(w.out))
          w.miss_Z_plus_A1$Z = dat$Z+kappa; w.miss_Z_plus_A1$A = rep(1, nrow(w.out))
          w.miss_Z_minus_A0$Z = dat$Z-kappa; w.miss_Z_minus_A0$A = rep(0, nrow(w.out))
          w.miss_Z_plus_A0$Z = dat$Z+kappa; w.miss_Z_plus_A0$A = rep(0, nrow(w.out))
         missmod = ranger(R ~., data = cbind(w.miss, R = dat$R)[slong!=split,])
          omega_Z_A1 = predict(missmod, data = w.miss_Z_A1)$predictions
          omega_Z_minus_A1 = predict(missmod, data = w.miss_Z_minus_A1)$predictions
          omega_Z_plus_A1 = predict(missmod, data = w.miss_Z_plus_A1)$predictions
          omega_Z_A0 = predict(missmod, data = w.miss_Z_A0)$predictions
          omega_Z_minus_A0 = predict(missmod, data = w.miss_Z_minus_A0)$predictions
          omega_Z_plus_A0 = predict(missmod, data = w.miss_Z_plus_A0)$predictions
        }else{
           missmod = ranger(R ~., data = cbind(x.miss, R = dat$R)[slong!=split,])
        omega_Z_A1 = predict(missmod, data = x.miss_Z_A1)$predictions
        omega_Z_minus_A1 = predict(missmod, data = x.miss_Z_minus_A1)$predictions
        omega_Z_plus_A1 = predict(missmod, data = x.miss_Z_plus_A1)$predictions
        omega_Z_A0 = predict(missmod, data = x.miss_Z_A0)$predictions
        omega_Z_minus_A0 = predict(missmod, data = x.miss_Z_minus_A0)$predictions
        omega_Z_plus_A0 = predict(missmod, data = x.miss_Z_plus_A0)$predictions
        }
       
     
      ## fit outcome models (time-varying)
      for(j in 1:k){
        dat$Y = as.integer(time.list[j] < dat$obs.Y)
       

        if("outcome" %in% misspecify){
          w.out_Z_A1 = w.out_Z_minus_A1 = w.out_Z_plus_A1 = w.out 
          w.out_Z_A0 = w.out_Z_minus_A0 = w.out_Z_plus_A0 = w.out 
          w.out_Z_A1$A = rep(1, nrow(w.out));  w.out_Z_A0$A = rep(0, nrow(w.out))
          w.out_Z_minus_A1$Z = dat$Z-kappa; w.out_Z_minus_A1$A = rep(1, nrow(w.out))
          w.out_Z_plus_A1$Z = dat$Z+kappa; w.out_Z_plus_A1$A = rep(1, nrow(w.out))
          w.out_Z_minus_A0$Z = dat$Z-kappa; w.out_Z_minus_A0$A = rep(0, nrow(w.out))
          w.out_Z_plus_A0$Z = dat$Z+kappa; w.out_Z_plus_A0$A = rep(0, nrow(w.out))
            outmod = ranger(Y ~., data = cbind(w.out, Y = dat$Y)[slong!=split,])
            mu_Z_A1 = predict(outmod, data = w.out_Z_A1)$predictions
            mu_Z_minus_A1 = predict(outmod, data = w.out_Z_minus_A1)$predictions
            mu_Z_plus_A1 = predict(outmod, data = w.out_Z_plus_A1)$predictions
            mu_Z_A0 = predict(outmod, data = w.out_Z_A0)$predictions
            mu_Z_minus_A0 = predict(outmod, data = w.out_Z_minus_A0)$predictions
            mu_Z_plus_A0 = predict(outmod, data = w.out_Z_plus_A0)$predictions
          }else{       
            outmod = ranger(Y ~., data = cbind(x.out, Y = dat$Y)[slong!=split,])
            mu_Z_A1 = predict(outmod, data = x.out_Z_A1)$predictions
            mu_Z_minus_A1 = predict(outmod, data = x.out_Z_minus_A1)$predictions
            mu_Z_plus_A1 = predict(outmod, data = x.out_Z_plus_A1)$predictions
            mu_Z_A0 = predict(outmod, data = x.out_Z_A0)$predictions
            mu_Z_minus_A0 = predict(outmod, data = x.out_Z_minus_A0)$predictions
            mu_Z_plus_A0 = predict(outmod, data = x.out_Z_plus_A0)$predictions
          }
      
      
        phi1_1 = mu_Z_plus_A1*pi_Z_plus  + mu_Z_plus_A0*(1-pi_Z_plus) + 
           dat$R*delta_Z_minus*(dat$A*dat$Y - dat$A*mu_Z_A1)/(omega_Z_A1*delta_Z) + 
            delta_Z_minus*mu_Z_A1*(dat$A - pi_Z)/delta_Z +
           dat$R*delta_Z_minus*((1-dat$A)*dat$Y - (1-dat$A)*mu_Z_A0)/(omega_Z_A0*delta_Z) + 
           delta_Z_minus*mu_Z_A0*((1-dat$A) - (1-pi_Z))/delta_Z 


    
        phi1_0 = mu_Z_minus_A1*pi_Z_minus  + mu_Z_minus_A0*(1-pi_Z_minus) + 
           dat$R*delta_Z_plus*(dat$A*dat$Y - dat$A*mu_Z_A1)/(omega_Z_A1*delta_Z) + 
            delta_Z_plus*mu_Z_A1*(dat$A - pi_Z)/delta_Z +
           dat$R*delta_Z_plus*((1-dat$A)*dat$Y - (1-dat$A)*mu_Z_A0)/(omega_Z_A0*delta_Z) + 
           delta_Z_plus*mu_Z_A0*((1-dat$A) - (1-pi_Z))/delta_Z 
    
    
        phi2_1 = pi_Z_plus + 
            delta_Z_minus*(dat$A - pi_Z)/delta_Z
        phi2_0 = pi_Z_minus  + 
            delta_Z_plus*(dat$A - pi_Z)/delta_Z


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

