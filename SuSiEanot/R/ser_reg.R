

library(glmnet)

#' Single effect Bayesian regression model
#'
#' This function allows you to express your love of cats.
#' @param y vector of phenotypes standardized and centered.
#' @param X SNP matrix of dimenstion n by p (centered and standardized)
#' @param sigma_0 prior probablity for beta
#' @param sigma prior sigma
#' @param pi vector of length p represent prior probability of the SNPs
#' @param credible_prob probability of the credible set containing atleast 1 causal SNP
#' @examples
#' SER_fun()

SER_fun<-function(y,X,sigma,sigma_0,pi=NULL,credible_prob=0.95){
  
  if(is.null(pi)){ pi=rep(1/dim(X)[2],dim(X)[2])}
  p<-dim(X)[2]
  
  BF_vector<-mu<-sigma1_vector<-c()
  
  for(i in 1:p){
    model<-Bayes_simp(y=y,x=X[,i],sigma=sigma,sigma_0=sigma_0)
    mu[i]<-model$post_mu
    BF_vector[i]<-model$BF
    sigma1_vector[i]<-model$post_sigma
    
  }
  alpha<-BF_vector*pi/(sum(BF_vector*pi))
  ind<-order(alpha,decreasing = T)
  sorted_alpha<-alpha[ind]
  cum_alpha<-cumsum(sorted_alpha)
  sel_genes<-min(which(cum_alpha>credible_prob))
  CS<-ind[1:sel_genes]
  list('mu_post'=mu,'sigma_post'=sigma1_vector,'alpha'=alpha,'CS'=CS)
  
}