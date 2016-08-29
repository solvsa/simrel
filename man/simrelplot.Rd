\name{simrelplot}
\alias{simrelplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for plotting of an object of class simrel as output from a call to simrel().
}
\description{
The plotting function produces three plots:
1) A barplot of the true regression coefficents in the linear regression model.
2) Scree-plot of true eigenvalues (barplot) with true covariances (absolute values) between components and response overlayed (red dots). The covariances are scaled by the largest covariance (in absolute value).
3) Scree-plot of eigenvalues estimated from the simulated data (barplot) with estimated covariances (absolute values) between components 
and simulated response overlayed (red dots). The estimated covariances are scaled by the largest estimated covariance (in absolute value).
}
\usage{
  simrelplot(obj, ncomp = min(obj$p, obj$n, 20), ask = TRUE, print.cov=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ A simrel object. Output of call from function \code{simrel()}.
}
  \item{ncomp}{
The number of desired components to plot.
}
  \item{ask}{
Logical; if TRUE, the user is asked before displaying each plot separately, see \code{par(ask=.)}.
}
  \item{print.cov}{
Logical; if TRUE, then the absolute values of the estimated covariances between the latent components and the response are printed (non-scaled).
}
}

\details{
A plot of true eigenvalues of the covariance matrix with covariances between principal components and the response overlayed (plot number 2) 
displays important properties of the predictor variables and their relevance to the response variable. This is discussed in Sæbø et al. (2015).
The relevant components (indexed by \code{obj$relpos}) have non-zero covariances with the response. The size of the eigenvalues (variances) of the
relevant components is tightly connected to the performance of predictors. As discussed in Helland and Almøy (1994) prediction is relatively easy
if the relevant components also have large variances, whereas the opposite is typically the case if the relevant components have small variances (eigenvalues).
A corresponding plot of eigenvalues and covariances estimated from the simulated data (plot number 3) is included for comparison. This kind of plot
may in a real data setting give information about the number of relevant components and their positions (when ordered according to declining eigenvalues).
}

\references{
Helland, I. S. and Almøy, T., 1994, Comparison of prediction methods when only a few components are relevant, \emph{J. Amer. Statist. Ass.}, \bold{89}(426), 583 -- 591.

Sæbø, S., Almøy, T. and Helland, I. S., 2015, simrel - A versatile tool for linear model data simulation based on the concept of a relevant subspace and relevant predictors, \emph{Chemometr. Intell. Lab.}(in press),doi:10.1016/j.chemolab.2015.05.012.
}

\author{
Solve Sæbø
}

\examples{
mydata <- simrel(n=250, p=20, q=5, relpos=c(2,4), gamma=0.25, R2=0.75 )
simrelplot(mydata)

#Displaying single plots:
simrelplot(mydata, ask=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ simulations }
\keyword{ plot }% __ONLY ONE__ keyword per line
