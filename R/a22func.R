.a22func <- function(a11,a12,R2,R12,l1, l2){
  bb <- R12*a12/a11^2 * l1/l2
  root <- sqrt( R12^2*a12^2/a11^4 *l1^2/l2^2 - (1/l2 + l1/l2*a12^2/a11^2)*(l1/a11^2*R12^2 - R2))
  denom <- (1/l2+l1/l2*a12^2/a11^2)
  w1 <- (bb - root)/denom
  w2 <- (bb + root)/denom
  return(c(w1, w2))
}