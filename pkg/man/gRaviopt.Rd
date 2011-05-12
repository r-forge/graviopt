\name{gRaviopt}
\alias{gRaviopt}
\title{
Agent-based stochastic global optimization.
}
\concept{maximization}
\description{
Performs stochastic global optimization via Newtons laws of gravity and motion.
}
\usage{
gRaviopt(fn, Par, lower.limits = -10, upper.limits = 10, n = 20, m = 20, iterations = 200, man.scaling = FALSE, alpha = 0.1)
}
\arguments{
  \item{fn}{
  the function to be optimized (maximised). \code{NA} and \code{NaN} values are not allowed.
}
  \item{Par}{
number of parameters used in "fn".
}
  \item{lower.limits}{
the lower limits of the parameters used in "fn".
}
  \item{upper.limits}{
the upper limits of the parameters used in "fn".
}
  \item{n}{
number of agents that should be used for the optimization process.
}
  \item{m}{
number of best solutions to be stored in the memory.
}
  \item{iterations}{
number of iterations used for optimization process.
}
  \item{man.scaling}{
the user can specicify wether to choose the the radius of local optimization manually or automatically.
}
  \item{alpha}{
radius that divises between local and global optimization.
}
\value{
The output of the function \code{gRaviopt} contains a matrix \code{GMemory} with the \code{m} best solutions and an array \code{GP} of all calculations.
}
}
\references{
A. Kaveh and S. Talatahari: 
\emph{A novel heuristic optimization method: charged system search},
Acta Mech 213, 267--289 (2010)
}
\author{
Peter Kehler \email{peter.kehler.jr@googlemail.com}
}
\seealso{
\code{\link{gRaviopt.Plot}} for graphical output.
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

