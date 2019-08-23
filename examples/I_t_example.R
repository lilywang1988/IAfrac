# Examples for I_t
# install.packages("devtools")
# library(devtools)
# install_github("keaven/nphsim")
library(nphsim)
eps<-2 # delayed effect
p<-0.5 #treatment assignment
b<-30 # an intrinsic parameter to decide the number of intervals per time unit
tau<- 20 # end of the study
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
data_temp1 <-data.trim.d(n_event_FH,data_temp)[[1]]
I_t(data_temp1,data_temp1,rho,gamma) # the estimated information at the final stage
#Trim the data upto certain event numbers at the interim stage when 60% of the events have been observed. Have been trimed once to get data_temp1, no need to add additional variables, thus set the third argument to be F.

I_t.2(data_temp1,data_temp1,rho,gamma) # If we consider the change of the at-risk set, which is not necessary to be a fixed probability.

data_temp2 <- data.trim.d(ceiling(0.6*n_event_FH),data_temp1,F)[[1]]
I_t(data_temp1,data_temp2,rho,gamma) # Use the full dataset data_temp to provide the survival function, and check the estimated information for the trimmed data set data_temp2 with only 60% of the planned events have been observed.
I_t.2(data_temp1,data_temp2,rho,gamma) # If we consider the change of the at-risk set, which is not necessary to be a fixed probability.
