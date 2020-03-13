library(MASS)
source("Code/core/binary_survival_para_colo.R")
dat = read.table("Data/sample.csv", sep="," , header = TRUE)


time.list = seq(0, 8000, 100)
if_estimator = binary_survival_para_if_colo(dat = dat, time.list = time.list)

## discretize the observed failure times
dat_2 = dat
dat_2$obs.Y = trunc(dat$obs.Y/100)*100
if_hazard_estimator = binary_survival_para_if_hazard(dat = dat_2, time.list = time.list)
naive_estimator = binary_survival_para_naive(dat = dat_2, time.list = time.list)


