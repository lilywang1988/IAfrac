Fleming-Harrington class weighted log-rank tests and interim analysis
================
[Lili Wang](mailto:lilywang@umich.edu)
2019-08-22

Purpose
-------

This R package is to implement important works from Dr. Hasegawa (2014 and 2016) for Fleming-Harrington class weighted log-rank tests (FH-WLRT) sample size calculation and the information fraction monitoring for interim analysis (IA).

This R package will be improved and upgraded in the near future.

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
  tau<-66
  omega<-18
  eps<-6
  m1=21.7  #median survival time for placebo group
  m2=25.8  # median survival time for treatment group
  lambda<-log(2)/m1
  lambda.trt<-log(2)*(m1-eps)/(m2-eps)/m1
  theta=lambda.trt/lambda
  alpha<-0.025
  beta<-0.1
  rho=0
  gamma=1
  b=30
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
data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p), gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd
# First  under H1 obtain the maximum information at time tau, the end of the study
I_denom1<-I.1(rho, gamma,lambda,theta,eps,R,p,t.star=tau)*n_FH

# Second under H0, obtain the maximum information at time tau, the end of the study
I_denom0<-I.0(rho, gamma,lambda,R,p,t.star=tau)*n_FH

# Example 1: Set trimmed=F and work on the crude dataset from nphsim()
inf_frac_vec11<-FH.frac.cal(data_temp,c(10,15,18),I_denom1,rho,gamma,trimmed=F)
inf_frac_vec11
#> [1] 0.1836210 0.6139667 1.0497077

inf_frac_vec10<-FH.frac.cal(data_temp,c(10,15,18),I_denom0,rho,gamma,trimmed=F)
inf_frac_vec10
#> [1] 0.1501642 0.5020986 0.8584453

# Example 2: First trim the data before inputting into FH.frac.cal() setting trimmed=T, and obtain the whole spectrum.
tau.star=21 #in case the ratio=1 when t>tau
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
(interim_frac1<-inf_frac_vec21[interim_index1])
#> [1] 0.5888858
# WLRT at the interim under H0
interim_index0<- which.min(abs(inf_frac_vec21-0.6))
(interim_time0<-t_seq[interim_index0])
#> [1] 14.9
(interim_frac0<-inf_frac_vec20[interim_index0])
#> [1] 0.4815877
# WLRT at the final stage under H1
final_index1<- which.min(abs(inf_frac_vec21-1))
(final_time1<-t_seq[final_index1])
#> [1] 17.6
(final_frac1<-inf_frac_vec21[final_index1])
#> [1] 0.995822

# WLRT at the final stage under H0
final_index0<- which.min(abs(inf_frac_vec20-1))
(final_time0<-t_seq[final_index0])
#> [1] 19.3
(final_frac0<-inf_frac_vec20[final_index0])
#> [1] 0.9859277

# Example 3: Stop by event count d
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
#Generate data under H1, use eta=1e-5 to inhibit censoring
data_temp <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda.trt), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd

#Obtain the full information at the final stage based on the generated data
#Trim the data upto the final stage when n_event_FH events have been observed
data_temp1.ls <-data.trim.d(n_event_FH,data_temp)
data_temp1 <-data_temp1.ls[[1]]
I_t(data_temp1,data_temp1,rho,gamma) # the estimated information at the final stage
#> [1] 22.67557
I.1(rho,gamma,lambda,theta,eps,R,p,data_temp1.ls[[2]])*n_FH # predicted information
#> [1] 21.09362

#Generate data under H0, use eta=1e-5 to inhibit censoring
data_temp2 <- nphsim(nsim=1,lambdaC=lambda, lambdaE = c(lambda,lambda), ssC=ceiling(n_FH*(1-p)),intervals = c(eps),ssE=ceiling(n_FH*p),gamma=accrual.rt, R=R, eta=1e-5, fixEnrollTime = TRUE)$simd

#Trim the data upto the final stage when n_event_FH events have been observed
data_temp2.ls <-data.trim.d(n_event_FH,data_temp2)
data_temp2 <-data_temp2.ls[[1]]
I_t(data_temp2,data_temp2,rho,gamma) # estimated information at the final stage
#> [1] 22.24485
I.0(rho,gamma,lambda,R,p,data_temp2.ls[[2]])*n_FH # predicted information
#> [1] 23.55037
```

References
----------

1.  Kalbfleisch, J. D., & Prentice, R. L. (2011). The statistical analysis of failure time data (Vol. 360). John Wiley & Sons.

2.  Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.

3.  Hasegawa, T. (2016). Group sequential monitoring based on the weighted log‐rank test statistic with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 15(5), 412-419.

4.  Ye, T., & Yu, M. (2018). A robust approach to sample size calculation in cancer immunotherapy trials with delayed treatment effect. Biometrics, 74(4), 1292-1300.

5.  Luo, X., Mao, X., Chen, X., Qiu, J., Bai, S., & Quan, H. (2019). Design and monitoring of survival trials in complex scenarios. Statistics in medicine, 38(2), 192-209.
