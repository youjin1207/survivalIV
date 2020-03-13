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


hazard_funA = function(I = 1, Z, A, dat){ 
  fit = glm(Y ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + I + A, data = dat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                  colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                  diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  I = rep(1, nrow(dat)),
                                  A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}


censor_funA = function(I = 1, Z, A, dat){ 
  fit = glm(C ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + I + A, data = dat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                  colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                  diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  I = rep(1, nrow(dat)),
                                  A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
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


hazard_fun_naive = function(I = 1, A, dat){ 
  fit = glm(Y ~ A + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + I, data = dat, family = binomial())
  Y.fit = predict(fit, data.frame(A = rep(A, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                  colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                  diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  I = rep(1, nrow(dat))), type = "response")
  return(Y.fit)
}


censor_fun_naive = function(I = 1, A, dat){ 
  fit = glm(C ~ A + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + I, data = dat, family = binomial())
  Y.fit = predict(fit, data.frame(A = rep(A, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                  colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                  diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  I = rep(1, nrow(dat))), type = "response")
  return(Y.fit)
}




binary_survival_para_if_colo = function(dat, time.list){

	n = nrow(dat)
	k = length(time.list) # for treatment density
	ifvals = matrix(0, nrow = n, ncol = k)
	est.eff = rep(NA, k)
	
	# fit a treatment model
  pi_1 = pi_fun(Z = 1, dat = dat)
	pi_0 = rep(0, nrow(dat))

  # fit an instrumental variable model
	delta_1 = delta_fun(1, dat = dat)
	delta_0 = 1-delta_1
	
  # fit a censoring model
	omega_Z1_A1 = w_fun(Z = 1, A = 1, dat = dat)
	omega_Z1_A0 = w_fun(Z = 1, A = 0, dat = dat)
	omega_Z0_A0 = w_fun(Z = 0, A = 0, dat = dat)
		
	## fit outcome models (time-varying)
	for(j in 1:k){
		dat$Y = as.integer(time.list[j] < dat$obs.Y)
		
		mu_Z1_A1 = mu_fun(R = 1, Z = 1, A = 1, dat = dat)
		mu_Z1_A0 = mu_fun(R = 1, Z = 1, A = 0, dat = dat)
		mu_Z0_A0 = mu_fun(R = 1, Z = 0, A = 0, dat = dat)
		
    		 
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
		}
	

		# compute estimator (using)
		for(j in 1:k){
			est.eff[j] = mean(ifvals[,j]) # for each delta vaule
		}

		return(list(est = est.eff))
	}


binary_survival_para_if_hazard = function(dat, time.list){
	
	n = nrow(dat)
	k = length(time.list) # for treatment density
	ifvals = matrix(0, nrow = n, ncol = k)
	est.eff = rep(NA, k)

	# fit a treatment model
	pi_1 = pi_fun(Z = 1, dat = dat)	
	pi_0 = rep(0, nrow(dat))
	
	# fit an instrumental variable model
	delta_1 = delta_fun(1, dat = dat)
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
   
	## fit hazard models (time-varying)
	for(j in 2:k){
		
			dat$Y = as.integer(dat$obs.Y  == time.list[j] & dat$R == 1)
      dat$C = as.integer(dat$obs.Y  == time.list[j] & dat$R == 0)
      dat$I = as.integer(time.list[j-1] < dat$obs.Y)
			
			# fit a hazard function
			h_Z1_A1[j,] = hazard_funA(I = 1, Z = 1, A = 1, dat = dat)
			h_Z1_A0[j,] =	hazard_funA(I = 1, Z = 1, A = 0,  dat = dat)
			h_Z0_A0[j,] =	hazard_funA(I = 1, Z = 0, A = 0, dat = dat)
			

			S_tau_Z1_A1[j,] = S_tau_Z1_A1[j-1,]*(1-h_Z1_A1[j,])
  		S_tau_Z1_A0[j,] = S_tau_Z1_A0[j-1,]*(1-h_Z1_A0[j,])
  		S_tau_Z0_A0[j,] = S_tau_Z0_A0[j-1,]*(1-h_Z0_A0[j,])
			

			S_tau_Z1_A1[j,] = ifelse(S_tau_Z1_A1[j,] < 0.05, 0.05, S_tau_Z1_A1[j,])
			S_tau_Z1_A0[j,] = ifelse(S_tau_Z1_A0[j,] < 0.05, 0.05, S_tau_Z1_A0[j,])
			#S_tau_Z0_A1[j,] = ifelse(S_tau_Z0_A1[j,] < 0.05, 0.05, S_tau_Z0_A1[j,])
			S_tau_Z0_A0[j,] = ifelse(S_tau_Z0_A0[j,] < 0.05, 0.05, S_tau_Z0_A0[j,])

		  # fit a censoring hazard function
			omega_Z1_A1[j,] = censor_funA(I = 1, Z = 1, A = 1, dat = dat)
			omega_Z1_A0[j,] = censor_funA(I = 1, Z = 1, A = 0, dat = dat)
	  	omega_Z0_A0[j,] = censor_funA(I = 1, Z = 0, A = 0, dat = dat)
		

		
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
	}


	return(list(est = est.eff))
}


binary_survival_para_naive = function(dat, time.list){
	
	n = nrow(dat)
	k = length(time.list) 
	ifvals = matrix(0, nrow = n, ncol = k)
	est.eff = rep(NA, k)

	# fit treatment model
	pi_1 = pi_fun_naive(dat = dat)
	pi_1 = ifelse(pi_1 > 0.95, 0.95, ifelse(pi_1 < 0.05, 0.05, pi_1))
	
	
	 S_tau_A1 = S_tau_A0 = G_tau_A1 = G_tau_A0 = matrix(0, length(time.list), nrow(dat))
   h_A1 = h_A0 = omega_A1 = omega_A0 = matrix(0, length(time.list), nrow(dat))
   D_A1 = D_A0 = matrix(0, length(time.list), nrow(dat))
   ## when j = 1
   S_tau_A1[1,] = S_tau_A0[1,] = rep(1, nrow(dat))
   G_tau_A1[1,] = G_tau_A0[1,] = rep(1, nrow(dat))
   D_A1[1,] = D_A0[1,] = rep(1, nrow(dat))
   ifvals[, 1] = mean(D_A1[1,] - D_A0[1,], na.rm = TRUE) 
	

	 ## fit hazard models (time-varying)
	 for(j in 2:k){
  	
  		dat$Y = as.integer(dat$obs.Y == time.list[j] & dat$R == 1)
  		dat$C = as.integer(dat$obs.Y == time.list[j]  & dat$R == 0)
  		dat$I = as.integer(time.list[j-1]  < dat$obs.Y)
		
		
  		h_A1[j,] = hazard_fun_naive(I = 1, A = 1, dat = dat)
  		h_A0[j,] = hazard_fun_naive(I = 1, A = 0, dat = dat)
		  

  		S_tau_A1[j,] = S_tau_A1[j-1,]*(1-h_A1[j,])
  		S_tau_A0[j,] = S_tau_A0[j-1,]*(1-h_A0[j,])
      S_tau_A1[j,] = ifelse(S_tau_A1[j,] < 0.05, 0.05, S_tau_A1[j,])
  		S_tau_A0[j,] = ifelse(S_tau_A0[j,] < 0.05, 0.05, S_tau_A0[j,])


    	
    	omega_A1[j,] = censor_fun_naive(I = 1, A = 1, dat = dat)
    	omega_A0[j,] = censor_fun_naive(I = 1, A = 0, dat = dat)
    	

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
  	}

	return(list(est = est.eff))
}


