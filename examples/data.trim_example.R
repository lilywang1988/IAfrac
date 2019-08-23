# Examples for data.trim and data.trim.d
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

#Obtain the full information at the final stage based on the generated data
#Trim the data upto the final stage when n_event_FH events have been observed
data_final<-data.trim.d(n_event_FH,data_temp)[[1]]
data_final.time<-data.trim.d(n_event_FH,data_temp)[[2]]

data_final2<-data.trim(tau,data_temp)
