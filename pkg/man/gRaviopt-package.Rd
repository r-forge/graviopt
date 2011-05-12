\name{gRaviopt-package}
\alias{gRaviopt-package}
\alias{gRaviopt}
\docType{package}
\title{
gRaviopt
}
\description{
Implementation of an stochastic agent-based optimisation algorithm using the laws of gravity and motion.
Loosely based on the CSS algorithm of A. Kaveh and S. Talatahari described in
"A novel heuristic optimization method: charged system search", Acta Mech 213, 267-289 (2010).
}
\details{
\tabular{ll}{
Package: \tab gRaviopt\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-06\cr
License: \tab GNU \cr
LazyLoad: \tab yes\cr
}
}
\author{
Peter Kehler
Maintainer: <peter.kehler.jr@googlemail.com>
}
\references{
A. Kaveh and S. Talatahari: 
\emph{A novel heuristic optimization method: charged system search},
Acta Mech 213, 267--289 (2010)
}
\examples{
## Rastrigin02 function
## The function has a global maximum f(x) = 0 at the point (0,0).  
## gRaviopt searches for maxima of the objective function between
## lower and upper bounds on each parameter to be optimized. 

Rastrigin02  <- function(X){
  -((X[,1]^2 - 10*cos(2*pi*X[,1]^2) + 10) + (X[,2]^2 - 10*cos(2*pi*X[,2]^2) + 10))
}

## This version of the function is needed for gRaviopt.Plot
Rastrigin02.2d  <- function(x,y){
  -((x*x - 10*cos(2*pi*x) + 10) + (y*y - 10*cos(2*pi*y) + 10))
}

# optimization process of Rastrigin02
Rast02 <- gRaviopt(fn= Rastrigin02, Par=2, n=20, lower.limits = -3, upper.limits = 3,man.scaling=TRUE,alpha=0.05)

# the best solutions found
Rast02$Memory

# the movements of the particles during the optimization process
gRaviopt.Plot(fn= Rastrigin02.2d, gRaviopt.Result=Rast02, Par=2, iterations=200, n=20, lower.limits = -3, upper.limits = 3, Movements=TRUE,man.scaling=TRUE,alpha=0.1,Nice=FALSE)
}
\keyword{optimize}