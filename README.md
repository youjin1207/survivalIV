# Nonparametric Instrumental Variable Estimator for Survival Outcomes

## Data

Our application study is based on the Prostate, Lung, Colorectal and Ovarian (PLCO) Cancer Screening Trial (https://cdas.cancer.gov/plco/). The PLCO data is available to the general scientific community upon request. 


### Sample data

We uploaded the sample data ``Data/sample.csv`` that has the same data structure as the PLCO data we analyzed but was randomly manipulated from the original data, but this data is provided just for illustrative purpose not for reproducing the results.


## Code

## Core functions used for simulation

- Code/core/binary_survival_para_sim.R

- Code/core/binary_survival_ranger_sim.R 

- Code/core/continuous_survival_para_sim.R 

- Code/core/continuous_survival_rander_sim.R


## Core functions used for real data analysis

- Code/core/binary_survival_para_colo.R : provides three influence function-based estimators described in Application Section with parametric estimation.

- Code/core/binary_survival_ranger_colo.R : provides three influence function-based estimators described in Application Section with nonparametric estimation using R package ranger.

## Implementation functions

- Code/cox_binary.R : simulation for comparing three estimators (IF, IPW, Plug-in) under the covariate-dependent censoring  with survival times generated from a Cox proportional hazard model.

- Code/cox_additive.R : simulation for comparing three estimators (IF, IPW, Plug-in) under the covariate-dependent censoring with survival times generated from an additive hazards model.

- Code/scenario1.R : simulation for comparing three influence function-based estimators under scenario 1 (no unmeasured confounding and outcome-dependent censoring)

- Code/scenario2.R : simulation for comparing three influence function-based estimators under scenario 2 (unmeasured confounding and random censoring)

- Code/scenario3.R : simulation for comparing three influence function-based estimators under scenario 3 (no unmeasured confounding and administrative censoring)

- Code/realanalysis.R : code for real data analysis both using parametric and nonparametric estimation.

- Code/additive_continuous.R : simulation for comparing three estimators (IF, IPW, Plug-in) under non-ignorable censoring with survival times generated from an additive hazards model and a continuous IV.

## Instructions for general use

For parametric estimation, replace the covariates set in a regression function with user-specific covariates. For example, in ``mu_fun`` for a regression function in Code/core/binary_survival_para_colo.R, you can replace the covariates with your own covariates X1, X2, and X3.

```
## PLCO data

mu_fun = function(R = 1, Z, A, subdat, dat){ 

  fit = glm(Y ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + A + R, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                    colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                    diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  R = rep(1, nrow(dat)), A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}

## User-specific covariates of X1, X2, and X3

dat = data.frame(Y = Y, Z = Z, A = A, X1 = X1, X2 = X2, X3 = X3)
mu_fun = function(R = 1, Z, A, subdat, dat){ 

  fit = glm(Y ~ Z + age + sex + fh_cancer_yes + fh_cancer_no +
              colo_fh_yes + colo_fh_no + polyps_f_yes + polyps_f_no + 
              diabetes_f_yes + diabetes_f_no + A + R, data = subdat, family = binomial())
  Y.fit = predict(fit, data.frame(Z = rep(Z, nrow(dat)), age = dat$age, sex = dat$sex, 
                                  fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                                    colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                                  polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                                    diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                                  R = rep(1, nrow(dat)), A = rep(A, nrow(dat))), type = "response")
  return(Y.fit)
}
```

For nonparametric estimation, we can replace training data set of `data` used in `ranger` containing the covariates to be adjusted. For example, the training data set for the outcome model, `x.out`, requires the baseline covariates, instrument (`Z`), and treatment (`A`) variable to be included in `data`. 

```
## PLCO data
x.out =  data.frame(age = dat$age, sex = dat$sex, 
                fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                Z = dat$Z, A = dat$A)

## User-specific covariates of X1, X2, and X3
dat = data.frame(Y = Y, Z = Z, A = A, X1 = X1, X2 = X2, X3 = X3)
x.out = data.frame(X1 = dat$X1, X2 = dat$X2, X3 = dat$X3,
                Z = dat$Z, A = dat$A)

```