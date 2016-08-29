.Rfunc <- function(alpha1,alpha2,L){
  t(alpha1)%*%solve(diag(L))%*%alpha2
}