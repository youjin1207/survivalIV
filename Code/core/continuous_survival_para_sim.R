source("survIV/nuisance_continuous_survival_para.R")

continuous_survival_para_if = function(dat, kappa, time.list, misspecify = NULL, survival = "Cox"){
	
	n = nrow(dat)
	k = length(time.list) # for treatment density
	ifvals = matrix(nrow = n, ncol = k)
	est.eff = rep(NA, k)
	subdat = dat

	## fit treatment model
	if("treatment" %in% misspecify){
		pi_Z = pi_fun_wrong(Z = dat$Z, subdat = subdat, dat = dat)
  		pi_Z_minus = pi_fun_wrong(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		pi_Z_plus = pi_fun_wrong(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}else{
  		pi_Z = pi_fun(Z = dat$Z, subdat = subdat, dat = dat)
  		pi_Z_minus = pi_fun(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		pi_Z_plus = pi_fun(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}
	
	# adjust too small/large values
	pi_Z = ifelse(pi_Z < 0.01, 0.01, ifelse(pi_Z > 0.99, 0.99, pi_Z))
	pi_Z_minus = ifelse(pi_Z_minus < 0.01, 0.01, ifelse(pi_Z_minus > 0.99, 0.99, pi_Z_minus))
	pi_Z_plus = ifelse(pi_Z_plus < 0.01, 0.01, ifelse(pi_Z_plus > 0.99, 0.99, pi_Z_plus))
	
	## fit instrument model	
	if("instrument" %in% misspecify){
		delta_Z = delta_fun_continuous_wrong(Z = dat$Z, subdat = subdat, dat = dat)
  		delta_Z_minus = delta_fun_continuous_wrong(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		delta_Z_plus = delta_fun_continuous_wrong(Z = dat$Z + kappa, subdat = subdat, dat = dat)			
	}else{
		delta_Z = delta_fun_continuous(Z = dat$Z, subdat = subdat, dat = dat)
  		delta_Z_minus = delta_fun_continuous(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		delta_Z_plus = delta_fun_continuous(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}
	
	# adjust too small/large values
	delta_Z = ifelse(delta_Z < 0.01, 0.01, ifelse(delta_Z > 0.99, 0.99, delta_Z))
	delta_Z_minus = ifelse(delta_Z_minus < 0.01, 0.01, ifelse(delta_Z_minus > 0.99, 0.99, delta_Z_minus))
	delta_Z_plus = ifelse(delta_Z_plus < 0.01, 0.01, ifelse(delta_Z_plus > 0.99, 0.99, delta_Z_plus))
  	
  	## fit censoring model	
	if("censoring" %in% misspecify){
		omega_Z_A1 = w_fun_wrong(Z = dat$Z, A = 1, subdat = subdat, dat = dat)
		omega_Z_minus_A1 = w_fun_wrong(Z = dat$Z - kappa, A = 1, subdat = subdat, dat = dat)
		omega_Z_plus_A1 = w_fun_wrong(Z = dat$Z + kappa, A = 1, subdat = subdat, dat = dat)		  	
		omega_Z_A0 = w_fun_wrong(Z = dat$Z, A = 0, subdat = subdat, dat = dat)
		omega_Z_minus_A0 = w_fun_wrong(Z = dat$Z - kappa, A = 0, subdat = subdat, dat = dat)
		omega_Z_plus_A0 = w_fun_wrong(Z = dat$Z + kappa, A = 0, subdat = subdat, dat = dat)
	}else{
		omega_Z_A1 = w_fun(Z = dat$Z, A = 1, subdat = subdat, dat = dat)
		omega_Z_minus_A1 = w_fun(Z = dat$Z - kappa, A = 1, subdat = subdat, dat = dat)
		omega_Z_plus_A1 = w_fun(Z = dat$Z + kappa, A = 1, subdat = subdat, dat = dat)		  	
		omega_Z_A0 = w_fun(Z = dat$Z, A = 0, subdat = subdat, dat = dat)
		omega_Z_minus_A0 = w_fun(Z = dat$Z - kappa, A = 0, subdat = subdat, dat = dat)
		omega_Z_plus_A0 = w_fun(Z = dat$Z + kappa, A = 0, subdat = subdat, dat = dat)
	}


	## fit outcome model
  	if("outcome" %in% misspecify){
      newdata_Z_A1 = data.frame(Z = dat$Z, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A1 = data.frame(Z = dat$Z - kappa, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A1 = data.frame(Z = dat$Z + kappa, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_A0 = data.frame(Z = dat$Z, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A0 = data.frame(Z = dat$Z - kappa, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A0 = data.frame(Z = dat$Z + kappa, W.1 = dat$W.1, 
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
      newdata_Z_A1 = data.frame(Z = dat$Z, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A1 = data.frame(Z = dat$Z - kappa,X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A1 = data.frame(Z = dat$Z + kappa, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_A0 = data.frame(Z = dat$Z, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A0 = data.frame(Z = dat$Z - kappa, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A0 = data.frame(Z = dat$Z + kappa, X.1 = dat$X.1, 
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
	        mu_Z_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
	        mu_Z_minus_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_minus_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
	        mu_Z_plus_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_plus_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
	        mu_Z_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
	        mu_Z_minus_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_minus_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
	        mu_Z_plus_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_plus_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
		
		}else if(survival == "additive"){
		 	addbaseline = addbase$cumhaz[which.min(abs(addbase$times - time.list[j]))] # baseline
       		mu_Z_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       		mu_Z_minus_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_minus_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       		mu_Z_plus_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_plus_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       		mu_Z_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       		mu_Z_minus_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_minus_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
       		mu_Z_plus_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_plus_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))				
			
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
		ifvals[, j] = mean(phi1_1 - phi1_0) / mean(phi2_1 - phi2_0)
	}
	
	
	for(j in 1:k){
		est.eff[j] = mean(ifvals[,j]) # for each delta vaule
	}

	return(list(est = est.eff))

}


continuous_survival_para_ipw = function(dat, kappa, time.list, misspecify = NULL){
	

	n = nrow(dat)
	k = length(time.list) # for treatment density
	ifvals = matrix(nrow = n, ncol = k)
	est.eff = rep(NA, k)
	subdat = dat

	## fit treatment model
	if("treatment" %in% misspecify){
		pi_Z = pi_fun_wrong(Z = dat$Z, subdat = subdat, dat = dat)
  		pi_Z_minus = pi_fun_wrong(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		pi_Z_plus = pi_fun_wrong(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}else{
  		pi_Z = pi_fun(Z = dat$Z, subdat = subdat, dat = dat)
  		pi_Z_minus = pi_fun(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		pi_Z_plus = pi_fun(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}
	
	# adjust too small/large values
	pi_Z = ifelse(pi_Z < 0.01, 0.01, ifelse(pi_Z > 0.99, 0.99, pi_Z))
	pi_Z_minus = ifelse(pi_Z_minus < 0.01, 0.01, ifelse(pi_Z_minus > 0.99, 0.99, pi_Z_minus))
	pi_Z_plus = ifelse(pi_Z_plus < 0.01, 0.01, ifelse(pi_Z_plus > 0.99, 0.99, pi_Z_plus))
	
	## fit instrument model	
	if("instrument" %in% misspecify){
		delta_Z = delta_fun_continuous_wrong(Z = dat$Z, subdat = subdat, dat = dat)
  		delta_Z_minus = delta_fun_continuous_wrong(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		delta_Z_plus = delta_fun_continuous_wrong(Z = dat$Z + kappa, subdat = subdat, dat = dat)			
	}else{
		delta_Z = delta_fun_continuous(Z = dat$Z, subdat = subdat, dat = dat)
  		delta_Z_minus = delta_fun_continuous(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		delta_Z_plus = delta_fun_continuous(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}
	
	# adjust too small/large values
	delta_Z = ifelse(delta_Z < 0.01, 0.01, ifelse(delta_Z > 0.99, 0.99, delta_Z))
	delta_Z_minus = ifelse(delta_Z_minus < 0.01, 0.01, ifelse(delta_Z_minus > 0.99, 0.99, delta_Z_minus))
	delta_Z_plus = ifelse(delta_Z_plus < 0.01, 0.01, ifelse(delta_Z_plus > 0.99, 0.99, delta_Z_plus))
  	
  	## fit censoring model	
	if("censoring" %in% misspecify){
		omega_Z_A1 = w_fun_wrong(Z = dat$Z, A = 1, subdat = subdat, dat = dat)
		omega_Z_minus_A1 = w_fun_wrong(Z = dat$Z - kappa, A = 1, subdat = subdat, dat = dat)
		omega_Z_plus_A1 = w_fun_wrong(Z = dat$Z + kappa, A = 1, subdat = subdat, dat = dat)		  	
		omega_Z_A0 = w_fun_wrong(Z = dat$Z, A = 0, subdat = subdat, dat = dat)
		omega_Z_minus_A0 = w_fun_wrong(Z = dat$Z - kappa, A = 0, subdat = subdat, dat = dat)
		omega_Z_plus_A0 = w_fun_wrong(Z = dat$Z + kappa, A = 0, subdat = subdat, dat = dat)
	}else{
		omega_Z_A1 = w_fun(Z = dat$Z, A = 1, subdat = subdat, dat = dat)
		omega_Z_minus_A1 = w_fun(Z = dat$Z - kappa, A = 1, subdat = subdat, dat = dat)
		omega_Z_plus_A1 = w_fun(Z = dat$Z + kappa, A = 1, subdat = subdat, dat = dat)		  	
		omega_Z_A0 = w_fun(Z = dat$Z, A = 0, subdat = subdat, dat = dat)
		omega_Z_minus_A0 = w_fun(Z = dat$Z - kappa, A = 0, subdat = subdat, dat = dat)
		omega_Z_plus_A0 = w_fun(Z = dat$Z + kappa, A = 0, subdat = subdat, dat = dat)
	}
  		
  	for(j in 1:k){
  		subdat$Y = as.integer(time.list[j] < subdat$obs.Y)
  		dat$Y = as.integer(time.list[j] < dat$obs.Y)
  			
		phi1_1 = mean((dat$Y*dat$R*dat$A*delta_Z_minus/(omega_Z_A1*pi_Z*delta_Z)))*
    			mean((dat$A*delta_Z_minus/(delta_Z))) + 
    			mean((dat$Y*dat$R*(1-dat$A)*delta_Z_minus/(omega_Z_A0*(1-pi_Z)*delta_Z)))*
    			mean(((1-dat$A)*delta_Z_minus / (delta_Z)))
  
  		phi1_0 =  mean((dat$Y*dat$R*dat$A*delta_Z_plus/(omega_Z_A1*pi_Z*delta_Z)))*
    			mean((dat$A*delta_Z_plus/ (delta_Z))) + 
    			mean((dat$Y*dat$R*(1-dat$A)*delta_Z_plus/(omega_Z_A0*(1-pi_Z)*delta_Z)))*
    			mean(((1-dat$A)*delta_Z_plus/ (delta_Z)))
  
  		phi2_1 = dat$A*delta_Z_minus/delta_Z
  		phi2_0 = dat$A*delta_Z_plus/delta_Z
    	
      	## influence-function based estimator:
  		ifvals[ ,j] = mean(phi1_1 - phi1_0) / 
										mean(phi2_1 - phi2_0) 
	}
		

	for(j in 1:k){
		est.eff[j] = mean(ifvals[,j]) 
	}

	return(list(est = est.eff))
			
}


continuous_survival_para_plugin = function(dat, kappa, time.list, misspecify = NULL, survival = "Cox"){
	
	n = nrow(dat)
	k = length(time.list) # for treatment density
	ifvals = matrix(nrow = n, ncol = k)
	est.eff = rep(NA, k)
	subdat = dat

	## fit treatment model
	if("treatment" %in% misspecify){
		pi_Z = pi_fun_wrong(Z = dat$Z, subdat = subdat, dat = dat)
  		pi_Z_minus = pi_fun_wrong(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		pi_Z_plus = pi_fun_wrong(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}else{
  		pi_Z = pi_fun(Z = dat$Z, subdat = subdat, dat = dat)
  		pi_Z_minus = pi_fun(Z = dat$Z - kappa, subdat = subdat, dat = dat)
  		pi_Z_plus = pi_fun(Z = dat$Z + kappa, subdat = subdat, dat = dat)
	}
	
	# adjust too small/large values
	pi_Z = ifelse(pi_Z < 0.01, 0.01, ifelse(pi_Z > 0.99, 0.99, pi_Z))
	pi_Z_minus = ifelse(pi_Z_minus < 0.01, 0.01, ifelse(pi_Z_minus > 0.99, 0.99, pi_Z_minus))
	pi_Z_plus = ifelse(pi_Z_plus < 0.01, 0.01, ifelse(pi_Z_plus > 0.99, 0.99, pi_Z_plus))
	
	## fit outcome model
  	if("outcome" %in% misspecify){
      newdata_Z_A1 = data.frame(Z = dat$Z, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A1 = data.frame(Z = dat$Z - kappa, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A1 = data.frame(Z = dat$Z + kappa, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_A0 = data.frame(Z = dat$Z, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A0 = data.frame(Z = dat$Z - kappa, W.1 = dat$W.1, 
                           W.2 = dat$W.2, W.3 = dat$W.3, W.4 = dat$W.4, W.5 = dat$W.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A0 = data.frame(Z = dat$Z + kappa, W.1 = dat$W.1, 
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
      newdata_Z_A1 = data.frame(Z = dat$Z, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A1 = data.frame(Z = dat$Z - kappa,X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A1 = data.frame(Z = dat$Z + kappa, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(1, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_A0 = data.frame(Z = dat$Z, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_minus_A0 = data.frame(Z = dat$Z - kappa, X.1 = dat$X.1, 
                         X.2 = dat$X.2, X.3 = dat$X.3, X.4 = dat$X.4, X.5 = dat$X.5,
                           A = rep(0, nrow(dat)),  R = rep(1, nrow(dat)))
      newdata_Z_plus_A0 = data.frame(Z = dat$Z + kappa, X.1 = dat$X.1, 
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
        mu_Z_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z_minus_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_minus_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z_plus_A1 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_plus_A1) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z_minus_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_minus_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
        mu_Z_plus_A0 = exp(-baseline*exp(as.numeric(as.matrix(newdata_Z_plus_A0) %*% as.matrix(summary(coxfit)$coefficients[1:8,1]))))
	
	}else if(survival == "additive"){
	 	addbaseline = addbase$cumhaz[which.min(abs(addbase$times - time.list[j]))] # baseline
   		mu_Z_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
   		mu_Z_minus_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_minus_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
   		mu_Z_plus_A1 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_plus_A1) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
   		mu_Z_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
   		mu_Z_minus_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_minus_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))
   		mu_Z_plus_A0 = exp(-addbaseline - time.list[j]*as.numeric(as.matrix(newdata_Z_plus_A0) %*% as.matrix(summary(addfit)$coefficients[1:8,1])))				
		
	}
		 
		phi1_1 = mu_Z_plus_A1*pi_Z_plus  + mu_Z_plus_A0*(1-pi_Z_plus) 

	    phi1_0 = mu_Z_minus_A1*pi_Z_minus  + mu_Z_minus_A0*(1-pi_Z_minus) 
	    
	    phi2_1 = pi_Z_plus 
	    phi2_0 = pi_Z_minus  

	  	## influence-function based estimator:
		ifvals[, j] = mean(phi1_1 - phi1_0) /
	              mean(phi2_1 - phi2_0)
		
	}

	# compute estimator (using)
	for(j in 1:k){
		est.eff[j] = mean(ifvals[,j]) 
	}


	return(list(est = est.eff))
}
