\name{mbrdsim}
\alias{mbrdsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to set up a design for a given set of factors with their specific levels
using the MBR-design method.
}
\description{
The multi-level binary replacement (MBR) design approach is used here in order to facilitate the investigation of the effects of 
the data properties on the performance of estimation/prediction methods. The mbrdsim function
takes as input a list containing a set of factors with their levels. The output is an MBR-design
with the combinations of the factor levels to be run.
}
\usage{
	mbrdsim(simlist, fraction, gen=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simlist}{ A named list containing the levels of a set of (multi-level) factors.}
  \item{fraction}{Design fraction at bit-level. Full design: fraction=0, half-fraction: fraction=1, and so on.}
  \item{gen}{Generators for the fractioning at the bit level. Default is \code{NULL} for which the generators are chosen
automatically by the \code{FrF2} function. See documentation of \code{FrF2} for details on how to set the generators.}
}

\value{
  \item{BitDesign }{The design at bit-factor level. The object is of class design, as output from FrF2. Function design.info()
  can be used to get extra design info of the bit-design. The bit-factors are named.numbered if the input factor list is named.}
  \item{Design }{The design at original factor level, non-randomized. The factors are named if the input factor list is named.}%% ...
}
\references{
Martens, H., M�ge, I., T�ndel, K., Isaeva, J., H�y, M. and S�b��, S., 2010, Multi-level binary replacement (MBR) design for computer experiments in high-dimensional nonlinear systems, \emph{J, Chemom}, \bold{24}, 748--756.
}
\author{
Solve S�b�
}

\examples{
#Input: A list of factors with their levels (number of levels must be a multiple of 2).
simlist <- list(R2=c(0.5,0.9), 
                pos = c(12,45),
                gamma = c(0.1, 0.9),
                comp = 1:8,
                alpha = c(0.01, 0.05, 0.1, 0.25))

# A 1/8 design
des <- mbrdsim(simlist, fraction=3)
#Setting generators manually as interactions between base bit-level factors.
des <- mbrdsim(simlist, fraction=3, gen=list(c(1,2,4,5),c(2,3,5),c(3,4,5)))

#The MBRD-design at original factor level
des$Design

#Info about the bit-design including bit-level aliasing (and resolution if \code{gen = NULL})
#Not run
#library(DoE.base)
#design.info(des$BitDesign)
}

\keyword{ MBRD }
\keyword{ Design }% __ONLY ONE__ keyword per line
