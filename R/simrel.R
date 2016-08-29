simrel <- function (n, p, q, relpos, gamma, R2, ntest = NULL, muY = NULL, 
          muX = NULL, lambda.min =.Machine$double.eps, sim = NULL) 
{
  if (!is.null(sim)) {
    betaX <- sim$beta
    beta0 <- sim$beta0
    muY <- sim$muY
    muX <- sim$muX
    qpos <- sim$relpred
    p <- sim$p
    q <- sim$q
    gamma <- sim$gamma
    lambdas <- sim$lambda
    R2 <- sim$R2
    relpos <- sim$relpos
    minerror <- sim$minerror
    Sigma <- sim$Sigma
    R <- sim$Rotation
    warning(paste("All parameters are collected from the supplied 'sim' object. \n"))
  }
  m <- length(relpos)
  if (q < m) 
    stop(paste("the number of relevant predictors must at least be equal to", 
               m, "\n"))
  if (is.null(sim)) {
    irrelpos <- (1:p)[-relpos]
    extradim <- q - m
    qpos <- sort(c(relpos, sample(irrelpos, extradim, replace = F)))
    nu <- lambda.min*exp(-gamma)/(1-lambda.min)
    if(lambda.min<0 || lambda.min>=1){stop("Parameter lambda.min must be in the interval [0,1]\n")}
    lambdas <- (exp(-gamma*(1:p))+nu)/(exp(-gamma)+nu)    
    SigmaZ <- diag(lambdas)
    SigmaZinv <- diag(1/lambdas)
    Sigmazy <- matrix(0, p, 1)
    r <- runif(m, 0, 1) * sample(c(-1, 1), m, replace = TRUE)
    Sigmazy[relpos, ] <- sign(r) * sqrt(R2 * abs(r)/sum(abs(r)) * 
                                          lambdas[relpos])
    SigmaY <- 1
    Sigma <- rbind(c(SigmaY, t(Sigmazy)), cbind(Sigmazy, 
                                                SigmaZ))
    Q <- matrix(rnorm(q^2), q)
    Q <- scale(Q, scale = F)
    Rq <- qr.Q(qr(Q))
    R <- diag(p)
    R[qpos, qpos] <- Rq
    if (q < (p - 1)) {
      Q <- matrix(rnorm((p - q)^2), (p - q))
      Q <- scale(Q, scale = F)
      Rnq <- qr.Q(qr(Q))
      R[(1:p)[-qpos], (1:p)[-qpos]] <- Rnq
    }
    betaZ <- SigmaZinv %*% Sigmazy
    betaX <- R %*% betaZ
    beta0 <- 0
    if (!(is.null(muY))) {
      beta0 <- beta0 + muY
    }
    if (!(is.null(muX))) {
      beta0 <- beta0 - t(betaX) %*% muX
    }
    R2 <- t(Sigmazy) %*% betaZ
    minerror <- SigmaY - R2
  }
  pd <- all(eigen(Sigma)$values > 0)
  if (pd) {
    Sigmarot <- chol(Sigma)
    Ucal <- matrix(rnorm(n * (p + 1), 0, 1), nrow = n)
    U1cal <- Ucal %*% Sigmarot
    Y <- U1cal[, 1, drop = F]
    if (!(is.null(muY))) {
      Y <- Y + rep(muY, n)
    }
    Z <- U1cal[, 2:(p + 1)]
    X <- Z %*% t(R)
    if (!(is.null(muX))) {
      X <- sweep(X, 2, muX, "+")
    }
    colnames(X) <- as.character(1:p)
    if (!is.null(ntest)) {
      Utest <- matrix(rnorm(ntest * (p + 1), 0, 1), nrow = ntest)
      U1test <- Utest %*% Sigmarot
      TESTY <- U1test[, 1, drop = F]
      if (!(is.null(muY))) 
        TESTY <- TESTY + rep(muY, ntest)
      TESTZ <- U1test[, 2:(p + 1)]
      TESTX <- TESTZ %*% t(R)
      if (!(is.null(muX))) {
        TESTX <- sweep(TESTX, 2, muX, "+")
      }
      colnames(TESTX) <- as.character(1:p)
    }
    else {
      TESTX <- NULL
      TESTY <- NULL
    }
  }
  else {
    stop("Correlation matrix is not positive definit \n")
  }
  res <- list()
  res$call <- match.call()
  res$X <- X
  res$Y <- Y
  res$beta <- betaX
  res$beta0 <- beta0
  res$muY <- muY
  res$muX <- muX
  res$relpred <- qpos
  res$TESTX <- TESTX
  res$TESTY <- TESTY
  res$n <- n
  res$p <- p
  res$m <- m
  res$q <- q
  res$gamma <- gamma
  res$lambda <- lambdas
  res$R2 <- drop(R2)
  res$relpos <- relpos
  res$minerror <- minerror
  res$Sigma <- Sigma
  res$Rotation <- R
  res$type = "univariate"
  class(res) <- "simrel"
  return(res)
}