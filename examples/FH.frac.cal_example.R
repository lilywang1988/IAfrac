# Examples for FH.frac.cal
# install.packages("devtools")
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
accrual.rt<-n_FH/R # the needed arrual rate
#Generate data accordingly, use eta=1e-5 to inhibit censoring
data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt),
                    ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),
                    gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd

# Example 1 for FH.frac.cal: Set trimmed=F and work on the crude dataset from nphsim()
inf_frac_vec1<-FH.frac.cal(data_temp,c(10,15,18),I_denom,rho,gamma,trimmed=F)
inf_frac_vec1

# Example 2 for FH.frac.cal: First trim the data before inputting into FH.frac.cal() setting trimmed=T, and obtain the whole spectrum.
I_denom<-I.1(rho, gamma,lambda,theta,eps,R,p,t.star=tau)*n_FH
tau.star=21 #in case the ratio=1 when t>tau
#Trim the data
data_temp2 <-data.trim(tau.star,data_temp)
t_seq <- seq(0.1,tau_star,0.1) # the time series to check the information fraction
inf_frac_vec2<-FH.frac.cal(data_temp2,t_seq,I_denom,rho,gamma,trimmed=T)
# WLRT at the interim
interim_index<- which.min(abs(inf_frac_vec2-0.6))
interim_time<-t_seq[interim_index]
interim_frac<-inf_frac_vec2[interim_index]
# WLRT at the final
final_index<- which.min(abs(inf_frac_vec2-1))
final_time<-t_seq[final_index]
final_frac<-inf_frac_vec2[final_index]
