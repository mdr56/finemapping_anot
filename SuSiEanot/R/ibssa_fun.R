
library(glmnet)

#' IBSS function annotation
#'
#' This function allows you to express your love of cats.
#' @param y vector of phenotypes standardized and centered.
#' @param X SNP matrix of dimenstion n by p (centered and standardized)
#' @param sigma_0 prior probablity for beta
#' @param sigma prior sigma
#' @param pi vector of length p represent prior probability of the SNPs
#' @param credible_prob probability of the credible set containing atleast 1 causal SNP
#' @examples
#' IBSS_function()



IBSS_function_a<-function(y,X,L,pi=NULL,credible_prob=0.95,sigma=1,sigma_0=1){
  
  p<-dim(X)[2]
  
  mu_matrix<-b_matrix<-alpha_matrix<-sigma1_matrix<-matrix(0,nrow=p,ncol=L)
  init_b_matrix<-b_matrix
  repeat{
    for(l in 1:L){
      pred_val<-X%*%b_matrix
      if(L>2){
        r<-y-apply(pred_val[,-l],1,sum)}
      else{  r<-y-pred_val[,-l]}
      model<-SER_fun(y=r,X=X,sigma=sigma,sigma_0=sigma_0,pi=pi)
      mu_matrix[,l]<-model$mu_post
      alpha_matrix[,l]<-model$alpha
      sigma1_matrix[,l]<-model$sigma_post
      b_matrix[,l]<-alpha_matrix[,l]*mu_matrix[,l]
    }
    if(sum(abs(b_matrix-init_b_matrix))<0.05) break;
    init_b_matrix<-b_matrix
  }
  
  CS<-apply(alpha_matrix,2, function(x) cal_CS(x,credible_prob = credible_prob))
  PIP<-1-apply(1-alpha_matrix,1,prod)
  list('alpha_matrix'=alpha_matrix,'mu_matrix'=mu_matrix,'sigma1_matrix'=sigma1_matrix,'CS'=CS,'PIP'=PIP)
  
  
}