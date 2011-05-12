\name{gRaviopt.Plot}
\alias{gRaviopt.Plot}
\title{
graphical output for functions with two parameters
}
\description{
visualisises the motion of the agents during the optimization process
}
\usage{
gRaviopt.Plot(fn, gRaviopt.Result, Par = 2, iterations = 200, n = 20, lower.limits = -10, upper.limits = 10, Movements = FALSE, Nice = TRUE, man.scaling = F, alpha = 0.1)
}
\arguments{
  \item{fn}{
needs a 2D version of the function that is optimized with \code{gRaviopt}
}
  \item{gRaviopt.Result}{
graphical output
}
  \item{Par}{
number of parameters used. \code{Par} should always be two.
}
  \item{iterations}{
number of iterations used for the optimization process
}
  \item{n}{
number of agents used for the optimization
}
\item{lower.limits}{
the lower limits of the parameters used in "fn".
}
  \item{upper.limits}{
the upper limits of the parameters used in "fn".
}
  \item{Movements}{
should the movements be visualised or a static graphic of the optimized function 
}
  \item{Nice}{
if \code{FALSE} no velocity, force vectors or new positions are drawn
}
  \item{man.scaling}{
the user can specicify wether to choose the the radius of local optimization manually or automatically.
}
  \item{alpha}{
radius that divises between local and global optimization.
}
}
\value{
either graphical output of the optimized function (usind \code{persp} and \code{image}) or the movements of the agents during the optimization process
}
\references{
A. Kaveh and S. Talatahari: 
\emph{A novel heuristic optimization method: charged system search},
Acta Mech 213, 267--289 (2010)
}
\author{
Peter Kehler \email{peter.kehler.jr@googlemail.com}
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
