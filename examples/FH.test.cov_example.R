# Examples for FH.test.cov and FH.test.cor

#install.packages("devtools")
# library(devtools)
# install_github("keaven/nphsim")
library(nphsim)
eps<-2 # delayed effect
p<-0.5 #treatment assignment
b<-30 # an intrinsic parameter to decide the number of intervals per time unit
tau<- 18 # end of the study
R<-14 # accrual period [0,R]
omega<- tau-R
lambda<-log(2)/6 # control group risk hazard
theta<-0.7 # hazard ratio
lambda.trt<- lambda*theta #hazard after the change point for the treatment arm
rho<- 0 # parameter for the weights
gamma<-1 #parameter for the weights
alpha<-0.025 #type 1 error
beta<-0.1 #type 2 error
# First we decide the sample size:
size_FH <- sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
n_FH <-size_FH$n
n_event_FH<-size_FH$n_event
d_fixed<-ceiling(0.6*n_event_FH)
accrual.rt<-n_FH/R # the needed arrual rate
#Generate data accordingly, use eta=1e-5 to inhibit censoring
data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt),
                    ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),
                    gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd
# Example 1 for WLR.test.cov and WLR.test.cor: Fleming-Harrington family Weights
# I will let w1 be the default 1
# define a WLRT for w2 accodring to the rho and gamma defined above.
w2<-function(...){survKM_minus(...)^rho*(1-survKM_minus(...))^gamma}
data_interim<-data.trim.d(d_fixed,data_temp)[[1]] #data trimmed at the interim stage, the second enry on the list is the interim time, refer to function data.trim.d for details.
WLR.test.cov(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2)
WLR.test.cor(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2)
# The variance should be identical to the output from I_t and correlation is 1 if two weights are identical.
WLR.test.cov(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w1=w2,w2=w2)
I_t.2(data_interim,data_interim,rho,gamma) # Use the exact at-risk Set other than using the fixed esitmated treatment assignment probability used in Hasegawa (2016), which is implemented in I_t
I_t(data_interim,data_interim,rho,gamma) # Different from I_t.2, since it simplify the at-risk process by repplacing it with a fixed estimated treatment assignment probability.
WLR.test.cor(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w1=w2,w2=w2)

#Example 2 for WLR.test.cov and WLR.test.cor: any Weights
w2_2<-function(v,...){1-exp(-v*0.25)}
WLR.test.cov(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2_2)
WLR.test.cor(survival=data_interim$survival,delta=data_interim$delta,trt=data_interim$trt,w2=w2_2)




