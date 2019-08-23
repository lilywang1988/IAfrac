\name{avg.haz}
\alias{avg.haz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calcualte the average hazard ratios
}
\description{
Calculate the average hazard ratios according to Kalbfleisch and Prentice (1981) or in the paper Hasegawa (2014) for piece-wise exponential survival functions (only one change point \code{eps}).
}
\usage{
avg.haz(theta, eps, lambda, p = 1/2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
hazard ratio after eps between the treatment and the control group, assuming that the hazard rato is 1 before eps.
}
  \item{eps}{
The change point, before which, the hazard ratio is 1, and after which, the hazard ratio is theta.
}
  \item{lambda}{
The constant hazard for the control arm.
}
  \item{p}{
Treatment assignment probability.
}
}
\details{
This is to estimate the average hazard ratios when there exist non-proportional hazard ratios.
 \deqn{AHR(x) = \frac{\frac{1}{2p}(1-e^{-2p\lambda\epsilon})+\frac{\theta}{p(1+theta)}e^{-\lambda\epsilon(1+p-\theta(1-p))}}{\frac{1}{2p}(1-e^{-2p\lambda\epsilon})+\frac{1}{p(1+theta)}e^{-\lambda\epsilon(1+p-\theta(1-p))}}}{%
 [1/2/p*(1-exp(-2*p*lambda*eps))+theta/p/(1+theta)*exp(-lambda*eps*(1+p-theta*(1-p)))]/[1/2/p*(1-exp(-2*p*lambda*eps))+1/p/(1+theta)*exp(-lambda*eps*(1+p-theta*(1-p)))]
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kalbfleisch, J. D., & Prentice, R. L. (2011). The statistical analysis of failure time data (Vol. 360). John Wiley & Sons.

Hasegawa, T. (2014). Sample size determination for the weighted log‐rank test with the Fleming–Harrington class of weights in cancer vaccine studies. Pharmaceutical statistics, 13(2), 128-135.

}
\author{
Lili Wang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
# test
lambda=log(2)/6
theta=0.7
eps=2
avg.haz(theta,eps,lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
