cal_CS<-function(alpha,credible_prob=0.95){
  ind<-order(alpha,decreasing = T)
  sorted_alpha<-alpha[ind]
  cum_alpha<-cumsum(sorted_alpha)
  sel_genes<-min(which(cum_alpha>=credible_prob))
  CS<-ind[1:sel_genes]
  return(CS)
}