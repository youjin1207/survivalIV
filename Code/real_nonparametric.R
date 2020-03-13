library(MASS)
library(ranger)

source("Code/core/binary_survival_ranger_colo.R")
dat = read.table("Data/sample.csv", sep="," , header = TRUE)

time.list = seq(0, 8000, 100)
x.trt = data.frame(age = dat$age, sex = dat$sex, 
               fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
               colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
               polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
               diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
               Z = dat$Z)
x.out =  data.frame(age = dat$age, sex = dat$sex, 
                fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                Z = dat$Z, A = dat$A)
x.miss =  data.frame(age = dat$age, sex = dat$sex, 
                 fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                 colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                 polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                 diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                 Z = dat$Z, A = dat$A)
x.instrument =  data.frame(age = dat$age, sex = dat$sex, 
                       fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                       colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                       polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                       diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no)


x.trt.naive = data.frame(age = dat$age, sex = dat$sex, 
               fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
               colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
               polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
               diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no)
x.out.naive =  data.frame(age = dat$age, sex = dat$sex, 
                fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no, A = dat$A)
x.miss.naive = data.frame(age = dat$age, sex = dat$sex, 
                 fh_cancer_yes = dat$fh_cancer_yes, fh_cancer_no = dat$fh_cancer_no,
                 colo_fh_yes = dat$colo_fh_yes, colo_fh_no = dat$colo_fh_no, 
                 polyps_f_yes = dat$polyps_f_yes, polyps_f_no = dat$polyps_f_no,
                 diabetes_f_yes = dat$diabetes_f_yes, diabetes_f_no = dat$diabetes_f_no,
                 A = dat$A)



if_estimator = binary_survival_ranger_if_colo(dat = dat, x.trt = x.trt, x.out = x.out, 
x.miss = x.miss, x.instrument = x.instrument,
time.list = time.list, nsplits = 2)

## discretize the observed failure times
dat_2 = dat
dat_2$obs.Y = trunc(dat$obs.Y/100)*100
if_hazard_estimator = binary_survival_ranger_if_hazard_colo(dat = dat_2,  x.trt = x.trt, x.out = x.out, 
x.miss = x.miss, x.instrument = x.instrument, time.list = time.list, nsplits = 2)

naive_estimator = binary_survival_ranger_naive_colo(dat = dat_2,  x.trt = x.trt.naive, x.out = x.out.naive, 
x.miss = x.miss.naive, x.instrument = NULL, time.list = time.list, nsplits = 2)
