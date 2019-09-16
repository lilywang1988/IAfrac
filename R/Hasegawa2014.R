######## Fleming and Harrington #########
######## Ref. paper T. Hasegawa 2014 pharmaceutical statistics ##########
sample.size_FH<-function(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta){
  # every sequence starts from 0
  # b<-30
  z.alpha<-qnorm(p=alpha,lower.tail = FALSE)
  z.beta<-qnorm(p=beta,lower.tail = FALSE)
  n_sub<-floor(b*tau)
  t<-c(0,seq(1,n_sub)/b)
  h_1<-rep(lambda,(n_sub+1)) #control
  h_2<-c(rep(lambda,(eps*b)),rep(lambda.trt,n_sub-eps*b+1)) #treatment
  N_1<-rep((1-p),(n_sub+1))
  N_2<-rep(p,(n_sub+1))
  for(i in 1:(n_sub-1)){
    N_1[i+1]<-N_1[i]*(1-h_1[i]/b-(t[i]>omega)/b/(tau-t[i]))
    N_2[i+1]<-N_2[i]*(1-h_2[i]/b-(t[i]>omega)/b/(tau-t[i]))
  }
  N_1[n_sub+1]<-N_2[n_sub+1]<-0

  f_S_1<-function(x) exp(-lambda*x)
  f_S_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-(lambda*eps+lambda.trt*(x-eps)))
  #f_S_2_2<-function(x) (x<eps)*exp(-lambda*x)+(x>=eps)*exp(-eps*lambda.trt*(1/theta-1))*exp(-lambda.trt*x)
  S_1<-f_S_1(t)
  S_2<-f_S_2(t)
  S<-(1-p)*S_1+p*S_2
  D<-(h_1*N_1+h_2*N_2)/b
  theta_seq<-h_2/h_1
  phi<-N_2/N_1
  r<-S^rho*(1-S)^gamma
  num_vec<-D*r*(phi*theta_seq/(1+phi*theta_seq)-phi/(1+phi))
  den_vec<-D*r^2*phi/(1+phi)^2
  E.star_num<-sum(num_vec[1:n_sub])
  E.star_den<-sqrt(sum(den_vec[1:n_sub]))
  E.star<-E.star_num/E.star_den
  n<-(z.alpha+z.beta)^2/E.star^2
  n_event<-sum(D)*n
  return(list(n=ceiling(ceiling(n*p)/p), n_event= ceiling(ceiling(n_event*p)/p),E.star=E.star,sum_D=sum(D),den_vec=den_vec[1:n_sub],num_vec=num_vec[1:n_sub],time_vec=seq(1,n_sub)/b))
}
if(F){
  # Example 1 in the paper Hasegawa(2014)
  p<-2/3
  tau<-66
  omega<-18
  eps<-6
  #theta<-0.79 #percent of full trt effect
  #lambda<-log(2)/6 # median survival time 10 months
  #lambda.trt<-lambda*theta #full treatment effect
  m1=21.7
  m2=25.8
  lambda<-log(2)/m1#log(2)/21.7 # median survival time 10 months
  lambda.trt<-log(2)*(m1-eps)/(m2-eps)/m1 #full treatment effect
  (theta=lambda.trt/lambda)
  alpha<-0.025
  beta<-0.1
  rho=0
  gamma=1
  b=30
  res=sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
  length(res$den_vec)
}





#' The average hazard ratio calculation according to Hasegawa(2016) Phamaceutical Statistics paper
avg.haz<-function(theta, eps,lambda,p=1/2){
  term1<-1/2/p*(1-exp(-2*p*lambda*eps))+theta/p/(1+theta)*exp(-lambda*eps*(1+p-theta*(1-p)))
  term2<-1/2/p*(1-exp(-2*p*lambda*eps))+1/p/(1+theta)*exp(-lambda*eps*(1+p-theta*(1-p)))
  term1/term2
}
# test
#lambda=log(2)/6
#theta=0.7
#eps=2
#avg.haz(theta,eps,lambda)
