Fleming-Harrington class weighted log-rank tests and interim analysis
================
[Lili Wang](mailto:lilywang@umich.edu)
2019-08-23

Purpose
-------

This R package is to implement important contributions from Dr. Hasegawa (2014 and 2016) to calculate sample sizes and information fractions (IF) for Fleming-Harrington class weighted log-rank tests (FH-WLRT) in interim analysis (IA).

This R package will be improved gradually and we hope to release it to CRAN eventually when it is well established. Please let us know if you find anything to be corrected. We highly appreciate your feedback.

Model assumptions
-----------------

Let the survival function is *S*<sub>0</sub>(*t*) for both treatment arms (treatment and control) under the null, but *S*<sub>1</sub>(*t*) for the treatment in comparison with *S*<sub>0</sub>(*t*) for the control arm under the alternative.

We consider a simple case in current version of the R package, with a delayed (*ϵ*) treatment effect and can be described with a piece-wise exponential distribution is proposed.

$$\\begin{array}{cc}
 S\_0(t)= \\exp(-\\lambda t),   
& 
S\_1(t)=\\left\\{\\begin{array}{c l}
exp(-\\lambda t) &for \\ t\\leq \\epsilon;\\\\
c \\exp(-\\theta\\lambda t)& for \\ t&gt; \\epsilon. 
\\end{array}  \\right.
\\end{array}$$

Note that
*c* = exp(−(1 − *θ*)*λ**ϵ*)
, and the corresponding density functions are

$$\\begin{array}{cc}
 f\_0(t)= \\lambda\\exp(-\\lambda t),   
& 
f\_1(t)=\\left\\{\\begin{array}{c l}
\\lambda\\exp(-\\lambda t)=\\lambda S\_1(t) &for \\ t\\leq \\epsilon;\\\\
\\theta\\lambda c\\exp(-\\theta\\lambda t)=\\theta\\lambda S\_1(t)& for \\ t&gt; \\epsilon}. 
\\end{array}  \\right.
\\end{array}$$

Install the package
-------------------

To install the R package from Github, you will need to install another R package "devtools". Please uncomment the codes to install them.

``` r
# install.packages("devtools")
# library(devtools)
# install_github("lilywang1988/IAfrac")
library(IAfrac)
#> Loading required package: data.table
#> Loading required package: survival
```

Vignette 1: sample size caulcation for FH-WLRT
----------------------------------------------

``` r
# Sample size calculation using an Example 1 from Hasegawa (2014)
 p<-2/3
  tau<-66 # The end of the study
  omega<-18 # tau-R, and R is the end of the arrual time, thus omega is can be viewed as the time between the end of the study and the last accrual. 
  eps<-6 # The change point, where under the H0, the HR between the treatment and the control group changes from 1 to theta
  m1=21.7  # median survival time for the control group
  m2=25.8  # median survival time for the treatment group
  lambda<-log(2)/m1 # hazard for the concontrol group and the treatment group before eps, the change point. 
  lambda.trt<-log(2)*(m1-eps)/(m2-eps)/m1 # hazaed for the treatment group after eps, the change point. 
  theta=lambda.trt/lambda # the hazard ratio after the change point
  alpha<-0.025 # the type I error under control, which is one-sided here
  beta<-0.1 # the type II error under control
  rho=0 # the first parameter for the Fleming-Harrington weight: S(t^-)^\rho (1-S(t^-))^\gamma. 
  gamma=1 # the second parameter for the Fleming-Harrington weight: S(t^-)^\rho (1-S(t^-))^\gamma. 
  b=30 # an intrinsic parameter to decide the number of intervals per time unit
  sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)$n
#> [1] 1974
 #1974, identical to the paper's report
```

Vignette 2: data manipulation and information fraction
------------------------------------------------------

``` r

# install.packages("devtools")
# library(devtools)
# install_github("keaven/nphsim")
library(nphsim)
#> Loading required package: survMisc
#> Loading required package: mvtnorm
#> Loading required package: Matrix
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:data.table':
#> 
#>     between, first, last
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: survminer
#> Loading required package: ggplot2
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:survMisc':
#> 
#>     autoplot
#> Loading required package: ggpubr
#> Loading required package: magrittr
#> Loading required package: survRM2
#> Warning: replacing previous import 'ggplot2::autoplot' by
#> 'survMisc::autoplot' when loading 'nphsim'
set.seed(123456)
eps<-2 # delayed effect
p<-0.5 #treatment assignment
b<-30 # an intrinsic parameter to decide the number of intervals per time unit
tau<- 18 # end of the study
R<-14 # accrual period [0,R]
omega<- tau-R # tau-R, and R is the end of the arrual time, thus omega is can be viewed as the time between the end of the study and the last accrual. 
lambda<-log(2)/6 # control group risk hazard
theta<-0.7 # hazard ratio
lambda.trt<- lambda*theta #hazard after the change point for the treatment arm
rho<- 0 # parameter for the weights
gamma<-1 #parameter for the weights
alpha<-0.025 #type 1 error
beta<-0.1 #type 2 error
# First we decide the sample size:
size_FH <- sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
n_FH <-size_FH$n # number of patients needed
n_event_FH<-size_FH$n_event # number of events needed, n_event_FH~n_FH*size_FH$sum_D
accrual.rt<-n_FH/R # the needed arrual rate
#Generate data accordingly

# use eta=1e-5 to inhibit censoring when generating the data
data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p), gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd
# First  under H1 obtain the maximum information at time tau, the end of the study
I_denom1<-I.1(rho, gamma,lambda,theta,eps,R,p,t.star=tau)*n_FH

# Second under H0, obtain the maximum information at time tau, the end of the study
I_denom0<-I.0(rho, gamma,lambda,R,p,t.star=tau)*n_FH

# Example 1: Set trimmed=F and work on the crude dataset from nphsim()
# Pick 3 time points, 10, 15, 18 to check their information fractions under the two hypothesis
inf_frac_vec11<-FH.frac.cal(data_temp,c(10,15,18),I_denom1,rho,gamma,trimmed=F)
inf_frac_vec11
#> [1] 0.1836210 0.6139667 1.0497077

inf_frac_vec10<-FH.frac.cal(data_temp,c(10,15,18),I_denom0,rho,gamma,trimmed=F)
inf_frac_vec10
#> [1] 0.1501642 0.5020986 0.8584453

# Example 2: First trim the data before inputting into FH.frac.cal() setting trimmed=T, and obtain the whole spectrum.
tau.star=21 # set a pseudo ending point in case IF=1 when t>tau
#Trim the data according by time tau.star
data_temp2 <-data.trim(tau.star,data_temp)
t_seq <- seq(0.1,tau.star,0.1) # the time series to check the information fraction
# under H1
inf_frac_vec21<-FH.frac.cal(data_temp2,t_seq,I_denom1,rho,gamma,trimmed=T)
# under H0
inf_frac_vec20<-FH.frac.cal(data_temp2,t_seq,I_denom0,rho,gamma,trimmed=T)
# WLRT at the interim under H1
interim_index1<- which.min(abs(inf_frac_vec21-0.6))
(interim_time1<-t_seq[interim_index1])
#> [1] 14.9
(interim_frac1<-inf_frac_vec21[interim_index1]) # should be close to 0.6
#> [1] 0.5888858
# WLRT at the interim under H0
interim_index0<- which.min(abs(inf_frac_vec20-0.6))
(interim_time0<-t_seq[interim_index0])
#> [1] 15.9
(interim_frac0<-inf_frac_vec20[interim_index0]) # should be close to 0.6
#> [1] 0.5972703
# WLRT at the final stage under H1
final_index1<- which.min(abs(inf_frac_vec21-1)) 
(final_time1<-t_seq[final_index1])
#> [1] 17.6
(final_frac1<-inf_frac_vec21[final_index1])  # should be close to 1
#> [1] 0.995822

# WLRT at the final stage under H0
final_index0<- which.min(abs(inf_frac_vec20-1))
(final_time0<-t_seq[final_index0])
#> [1] 19.3
(final_frac0<-inf_frac_vec20[final_index0])  # should be close to 1
#> [1] 0.9859277

# Example 3: Stop by event count d, note that data.trim.d returns a list of two components, the first is the dataset, the second is the stopping time
data_temp3.ls<-data.trim.d(d=200,data_temp)
data_temp3<-data_temp3.ls[[1]]
(interim_time_d<-data_temp3.ls[[2]]) # the interim time with d observed events
#> [1] 8.800698
```

Vignette 3: information prediction and estimation
-------------------------------------------------

``` r
# install.packages("devtools")
# library(devtools)
# install_github("keaven/nphsim")
library(nphsim)
set.seed(123456)
eps<-2 # delayed effect
p<-0.5 #treatment assignment
b<-30 # an intrinsic parameter to decide the number of intervals per time unit
tau<- 18 # end of the study
R<-14 # accrual period [0,R]
omega<- tau-R # tau-R, and R is the end of the arrual time, thus omega is can be viewed as the time between the end of the study and the last accrual. 
lambda<-log(2)/6 # control group risk hazard
theta<-0.7 # hazard ratio
lambda.trt<- lambda*theta #hazard after the change point for the treatment arm
rho<- 0 # parameter for the weights
gamma<-1 #parameter for the weights
alpha<-0.025 #type 1 error
beta<-0.1 #type 2 error
# First we decide the sample size:
size_FH <- sample.size_FH(eps,p,b,tau,omega,lambda,lambda.trt,rho, gamma,alpha,beta)
n_FH <-size_FH$n # number of patients needed
n_event_FH<-size_FH$n_event # number of events needed, n_event_FH~n_FH*size_FH$sum_D
accrual.rt<-n_FH/R # the needed arrual rate
#Generate data accordingly

# use eta=1e-5 to inhibit censoring when generating the data
data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd

#Obtain the full information at the final stage based on the generated data
#Trim the data up to the final stage when n_event_FH events have been observed
data_temp1.ls <-data.trim.d(n_event_FH,data_temp)
data_temp1 <-data_temp1.ls[[1]]
I_t(data_temp1,data_temp1,rho,gamma) # the estimated information of data_Temp1 at the final stage using data_temp1 to obtain the survival function estimates
#> [1] 22.67557
(data_temp1_t<-data_temp1.ls[[2]]) # the true stopping time based on the data data_temp1
#> [1] 17.76395

# Now if we would like to predict the stopping time
t_ls<-seq(0.1,25,0.1)
(pred_t1<-t_ls[which.min(abs(size_FH$sum_D-sapply(t_ls,function(t) {h.tilde(0,lambda,theta,eps,R,p,t)})))]) #predicted data_temp1_t, they should be similar
#> [1] 18

I.1(rho,gamma,lambda,theta,eps,R,p,pred_t1)*n_FH # predicted information, should be close to I_t(data_temp1,data_temp1,rho,gamma) above.
#> [1] 21.73793

#Generate data under H0, use eta=1e-5 to inhibit censoring
data_temp0 <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd

#Trim the data upto the final stage when n_event_FH events have been observed
data_temp0.ls <-data.trim.d(n_event_FH,data_temp0)
data_temp0 <-data_temp0.ls[[1]]
I_t(data_temp0,data_temp0,rho,gamma) # estimated information at the final stage
#> [1] 22.24485
(data_temp0_t<-data_temp0.ls[[2]]) # the true stopping time based on the data data_temp1
#> [1] 17.05209

# Now if we would like to predict the stopping time
t_ls<-seq(0.1,25,0.1)
(pred_t0<-t_ls[which.min(abs(size_FH$sum_D-sapply(t_ls,function(t) {h.tilde(0,lambda,1,eps,R,p,t)})))]) #predicted data_temp1_t, they should be similar
#> [1] 16.8
I.0(rho,gamma,lambda,R,p,pred_t0)*n_FH # predicted information, should be close to I_t(data_temp0,data_temp0,rho,gamma) above
#> [1] 22.75341
```

References
----------

1.  Kalbfleisch, J. D., & Prentice, R. L. (2011). The statistical analysis of failure time data (Vol. 360). John Wiley & Sons.

2.  Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.

3.  Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.

4.  Ye, T., & Yu, M. (2018). A robust approach to sample size calculation in cancer immunotherapy trials with delayed treatment effect. Biometrics, 74(4), 1292-1300.

5.  Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019). Design and monitoring of survival trials in complex scenarios. Statistics in medicine, 38(2), 192-209.
