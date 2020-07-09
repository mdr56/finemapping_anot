


library(glmnet)

#' Susie_anot_penal
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




anot_susie_penal<-function(y,X,L,pi=NULL,credible_prob=0.95,sigma=1,sigma_0=1,z_cov){
  c_p<-credible_prob
  sigma<-sigma
  sigma_0<-sigma_0
  beta_prev<-rep(0,dim(z_cov)[2]+1)
  iter<-0
  repeat{
    result_IBSS<-IBSS_function_a(y,X,L=L,pi=pi,credible_prob=c_p)
    prob_response<-result_IBSS$PIP
    pred_caus<-rbinom(length(prob_response),size=1,prob_response)
    model_prob<-cv.glmnet(z_cov,pred_caus,family='binomial')
    if(max(abs(coef(model_prob)-beta_prev))<0.1) break;
    print(max(abs(coef(model_prob)-beta_prev)))
    pi<-predict(model_prob,z_cov,type='response')
    beta_prev<-coef(model_prob)
    print(beta_prev)
    iter<-iter+1
    if(iter>10) break;
  }
  list('result_Susie'=result_IBSS,'pi'=pi,'model_logit'=model_prob)
  
}
