# Examples for approx.I
eps<-2 # delayed effect
p<-0.5 #treatment assignment
tau<-18 # end of the study
R<-14 # accrual period [0,R]
lambda<-log(2)/6 # control group risk hazard
theta<-0.7 # hazard ratio
lambda.trt<- lambda*theta
rho<- 0 # parameter for the weights
gamma<-1 #parameter for the weights
S1<-function(x){
   ifelse(x>eps,exp(-theta*lambda*x)*getc(theta,lambda,eps),exp(-lambda*x))
}
S0<-function(x){
    exp(-lambda*x)
}
S_pool<-function(x){
  p*S1(x)+(1-p)*S0(x)
}
func<-function(x){
  min((tau-x)/R,1)*(S_pool(x)^rho*(1-S_pool(x))^gamma)^2
}

approx.I(t.star=tau,p,S1=S1,S0=S0,fun=func,n.length=1e6)
I.1(rho,gamma,lambda,theta,eps,R,p,tau)

func2<-function(x){
  min((10-x)/R,1)*(S_pool(x)^rho*(1-S_pool(x))^gamma)^2
}


approx.I(t.star=10,p,S1=S1,S0=S0,fun=func2,n.length=1e6)
I.1(rho,gamma,lambda,theta,eps,R,p,t.star=10)

# Covariance approximation for two weights: 1 and G(0,1)
rho1=rho2=0
gamma1=0
gamma2=1
func3<-function(x){
  min((10-x)/R,1)*(S_pool(x)^rho1*(1-S_pool(x))^gamma1)*(S_pool(x)^rho2*(1-S_pool(x))^gamma2)
}
approx.I(t.star=10,p,S1=S1,S0=S0,fun=func3,n.length=1e6)
I.1.cov(rho1,gamma1,rho2,gamma2,lambda,theta,eps,R,p,t.star=10)
