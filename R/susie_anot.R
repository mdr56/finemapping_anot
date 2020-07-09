


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
#' @param L Number of credible set
#' @param z_cov matrix of annotation
#' @examples
#' anot_susie()



anot_susie<-function(y,X,L,pi=NULL,credible_prob=0.95,sigma=1,sigma_0=1,z_cov){
  c_p<-credible_prob
  sigma<-sigma
  sigma_0<-sigma_0
  #beta_prev<-rep(0,dim(z_cov)[2]+1)
  beta_prev<-rep(0,2)
  iter<-1
  repeat{
    result_IBSS<-IBSS_function_a(y,X,L=L,pi=pi,credible_prob=c_p)
    prob_response<-result_IBSS$PIP
    
    model_prob<-glm(prob_response~z_cov,family=binomial(link='logit'))
    if(max(abs(model_prob$coefficients-beta_prev))<0.1) break;
    pi<-predict(model_prob,type='response')
    beta_prev<-model_prob$coefficients
    iter<-iter+1
    if(iter>10) break;
  }
  list('result_Susie'=result_IBSS,'pi'=pi,'model_logit'=model_prob)
  
}
