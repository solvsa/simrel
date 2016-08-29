simrelplot <- function(obj, ncomp=min(obj$p,obj$n,20), ask=TRUE, print.cov=FALSE){

  def.par <- par(no.readonly = TRUE)
  if(!ask){
    layout(matrix(c(1,1,2,3),2,2, byrow=TRUE))
  }
  par(mar=c(5.1, 4.1, 4.1, 4.1))
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  dev.hold()
  plot(obj$beta[,1], type="h",col=2, lwd=1, xlab="Variable number", ylab=expression(beta),
       main="True regression coefficients",ylim=range(obj$beta),xlim=c(1,obj$p*1.1))
  abline(h=0, col=1)
  if(obj$type=="bivariate"){
    points(obj$beta[,2], type="h", lwd=1, col=3)
    legend(x="bottomright",legend=c(expression(beta[1]),expression(beta[2])),lty=1,col=c(2,3), bty="n",text.width=0.5)
    abline(h=0, col=1)
  }
  dev.flush()
  
  if(obj$type=="univariate"){
    covs <- abs(obj$Sigma[2:obj$p,1])
    covs.sc <- covs/max(covs)
  }
  else if(obj$type=="bivariate"){
    covs <- abs(obj$Sigma[3:obj$p,1])
    covs.sc <- covs/max(covs)
    covs2 <- abs(obj$Sigma[3:obj$p,2])
    covs2.sc <- covs2/max(covs2)
  }
  dev.hold()
    plot(obj$lambda[1:ncomp], type="h", lwd=2, col=1, 
       main="Relevant components plot",
       xlab="Components", ylab="Eigenvalue", axes=F, ylim=c(0,1), cex.lab=1)
    points(1:ncomp, covs.sc[1:ncomp], type="p", pch=20, cex=2, col=2)
  if(obj$type=="bivariate"){
    points(1:ncomp, covs2.sc[1:ncomp], type="p", pch=20, cex=2, col=3)
    legend(x="topright",legend=c(expression(Y[1]),expression(Y[2])),pch=20,pt.cex=2,col=c(2,3), bty="n",text.width=1)
  }
    axis(1)
    axis(2,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
    axis(4,at=seq(0,1,0.1), labels=as.character(round(seq(0,max(covs),length.out=length(seq(0,1,0.1))),2)))
    mtext("Covariance (absolute value)",side=4, line=3, cex=0.85)
    box()
  dev.flush()
  
  X <- scale(obj$X, center=TRUE, scale=FALSE)
  Y <- scale(obj$Y[,1], center=TRUE, scale=FALSE)
  if(obj$type=="bivariate"){
    Y2 <- scale(obj$Y[,2], center=TRUE, scale=FALSE)
  }
  svdres <- svd(X)
  eigval <- (svdres$d^2)/(obj$n-1)#/(svdres$d^2)[1]
  eigval.sc <- eigval/eigval[1]
  Z <- X%*%svdres$v
  covs <- abs(cov(Y, Z))
  covs.sc <- covs/max(abs(covs))
  if(obj$type=="bivariate"){
    covs2 <- abs(cov(Y2, Z))
    covs2.sc <- covs2/max(abs(covs2))
  }
  dev.hold()
    plot(eigval.sc[1:ncomp], type="h", lwd=2, col=1, 
       main="Estimated relevant components plot",
       xlab="Components", ylab="Eigenvalue", axes=F, ylim=c(0,1), cex.lab=1)
    points(1:ncomp, covs.sc[1:ncomp], type="p", pch=20, cex=2, col=2)
  if(obj$type=="bivariate"){
    points(1:ncomp, covs2.sc[1:ncomp], type="p", pch=20, cex=2, col=3)
    legend(x="topright",legend=c(expression(Y[1]),expression(Y[2])),pch=20,pt.cex=2,col=c(2,3), bty="n",text.width=1)
  }
    axis(1)
    #axis(2,at=seq(0,1,0.1), labels=as.character(seq(0,1,0.1)))
    axis(2,at=seq(0,max(eigval.sc),0.1), labels=as.character(round(seq(0,max(eigval),length.out=length(seq(0,max(eigval.sc),0.1))),2)))
    axis(4,at=seq(0,max(covs.sc),0.1), labels=as.character(round(seq(0,max(covs),length.out=length(seq(0,max(covs.sc),0.1))),2)))
    mtext("Covariance (absolute value)",side=4, line=3, cex=0.85)
    box()
  dev.flush() 
  if(print.cov){
    covs <- covs[1,1:min(ncomp, obj$p)]
    cat("Absolute values of estimated covariances\n")
    names(covs) <- paste("Component", 1:min(ncomp,obj$p),sep="")
    print(abs(round(covs,3)))
  }
  par(def.par)
}