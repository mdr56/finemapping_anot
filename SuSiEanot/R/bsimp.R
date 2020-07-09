




library(glmnet)

#' Simple Bayesian linear regression model
#'
#' This function allows you to express your love of cats.
#' @param y vector of phenotypes standardized and centered.
#' @param X SNP matrix of dimenstion n by p (centered and standardized)
#' @param sigma_0 prior probablity for beta
#' @param sigma prior sigma
#' @examples
#' Bayes_simp()

Bayes_simp<-function(y,x,sigma,sigma_0){
  
  b_hat<-solve(t(x)%*%x)%*%t(x)%*%y
  s<-sigma/(t(x)%*%x)
  sigma_1<-(1/s+1/sigma_0)^(-1)
  mu_1<-(sigma_1/s)*b_hat
  z<-b_hat/sqrt(s)
  
  BF<-sqrt(s/(s+sigma_0))*exp((z^2/2)*(sigma_0/(sigma_0+s)))
  
  list('post_sigma'=sigma_1,'b_hat'=b_hat,'post_mu'=mu_1,'BF'=BF)
  
}
