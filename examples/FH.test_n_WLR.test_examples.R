# Examples for FH.test and WLR.test
set.seed(12345)
data_temp<- nphsim(nsim=1,lambdaC=log(2)/6, lambdaE = c(log(2)/6,log(2)/6*0.7), ssC=250, intervals = c(2),ssE=250, gamma=500/14, R=14, eta=1e-5, fixEnrollTime = TRUE)$simd
data_final<-data.trim.d(100,data_temp)[[1]]
rho=1
gamma=0
# compare the 3 different ways below:
#library(survival)
sqrt(survdiff(Surv(survival,delta)~trt, data =data_final,rho=rho)$chisq)
FH.test(survival=data_final$survival,delta=data_final$delta,trt=data_final$trt,rho=rho,gamma=gamma)
WLR.test(survival=data_final$survival,delta=data_final$delta,trt=data_final$trt,w=function(...){survKM_minus(...)^rho*(1-survKM_minus(...))^gamma})

