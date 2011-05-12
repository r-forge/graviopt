gRaviopt.Plot <-
function(fn, gRaviopt.Result, Par=2, iterations=200, n=20, lower.limits = -10, upper.limits = 10, Movements=FALSE,Nice=TRUE,man.scaling=F, alpha=0.1){

if(Par!=2) stop("fn has ",Par," parameters instead of 2.")

GP <- gRaviopt.Result$GP
GM <- gRaviopt.Result$Memory

circles <- function(x, y, r, col = rep(0, length(x)),border = rep(1, length(x)), ...) {
    
    circle <- function(x, y, r, ...) {
    ang <- seq(0, 2*pi, length = 100)
    xx <- x + r * cos(ang)
    yy <- y + r * sin(ang)
    polygon(xx, yy, ...)
	}

    for(i in 1:length(x)) {
        circle(x[i], y[i], r[i], col = col[i], border = border[i], ...)
    }
}


# radius of charges: "a" in the article
Radius <- function(GP,k){
	if(man.scaling == TRUE){
		alpha*max(dist(GP[,1:Par,k]))
	} else {
		0.5*(1 - k/iterations)*max(dist(GP[,1:Par,k]))
	}
}

x <- seq(lower.limits,upper.limits,length.out=50)
y <- seq(lower.limits,upper.limits,length.out=50)
z <- outer(x,y,fn)

Ncol <- 50       	        # the number of colors to use
zlim <- c(min(z),max(z)) 	# limits in z coordinates
nlevels <- 25       		# see option nlevels in contour
theta <- 30         		# see option theta in persp
phi <- 30           		# see option phi in persp
		      

nrz <- nrow(z)
ncz <- ncol(z)

z1 <- z + max(abs(z))
Grenze <- max(z1) - min(z1)



couleurs <- tail(topo.colors(trunc(Ncol)),Ncol)
fcol <- couleurs[trunc(z1/Grenze*(Ncol-1))+1]
dim(fcol) <- c(nrz,ncz)
fcol <- fcol[-nrz,-ncz]

if(Movements==FALSE){
par(mfrow=c(1,2))
persp(x,y,z,col= fcol,zlim=zlim,theta=theta,phi=phi,ticktype = "detailed",xlab = "x", ylab = "y", zlab = "f(x,y)")

par(mar=c(2,2,2,2))
image(x,y,z,col= topo.colors(25))  # couleurs
contour(x,y,z,add=TRUE,nlevels=nlevels)
box()
}	


if(Movements==TRUE){
layout(matrix(c(1,1,1,2), 1, 4))
	for (i in 1:iterations){
		Sys.sleep(0.3)
		#par(ask=T)
		image(x,y,z,col=topo.colors(25),main = paste("Iteration",i,"of",iterations,"iterations."))
		contour(x,y,z,add=TRUE,nlevels=25)	

		points(GP[,c("x1","x2"),i],col="red",pch=21,bg="red")
		points(GM[,c("x1","x2")],col="black",pch=22)
		text(
		GP[,c("x1","x2"),i],
		labels=c(1:n),
		cex=1,
		adj=c(0,-0.3)
		)
		Rad <- Radius(GP=GP,k=i)
		if(Nice==TRUE){
		circles(x=GP[,"x1",i], y=GP[,"x2",i],r=rep(Rad,n))
		arrows(x0=GP[,"x1",i], y0=GP[,"x2",i], x1=(GP[,1,i]+GP[,"F_x1",i]), y1=(GP[,2,i]+GP[,"F_x2",i]),lty=2,col="blue") # force
		arrows(x0=GP[,"x1",i], y0=GP[,"x2",i], x1=GP[,1,i]+GP[,"v_x1",i], y1=GP[,2,i]+GP[,"v_x2",i],lty=2,col="red") # vel
		arrows(x0=GP[,"x1",i], y0=GP[,"x2",i], x1=GP[,"x1",i+1], y1=GP[,"x2",i+1],lty=2,col="green") # new Pos
		}
		plot.new()
	   	plot.window(c(2,2),c(2,2))
		legend("topleft",title="Legend",c(expression(sum(F[i], i==1, n)),expression(V[i]),expression(X[i+1]),"best results"),lty=c(2,2,2,NA),pch=c(NA,NA,NA,22),col=c("blue","red","green","black"),bg="white")
		mtext(paste("current best result:\n",round(GP[1,Par+1,i],6),"\nOverall best result\n",round(GM[1,Par+1])),side=3, line=-20)
		}
	}
}

