gRaviopt <-
function(
					fn, 
					Par, 
					lower.limits = -10, 
					upper.limits = 10, 
					n=20, m=20, 
					iterations=200,
					man.scaling=FALSE,
					alpha=0.1)
{


Names <- c(paste("x",1:Par,sep=""),"fn_xi",paste("v_x",1:Par,sep=""),paste("F_x",1:Par,sep=""))

# ---------------------------------------------------------------------------
# Initialization
# starting points
X <- matrix(runif(Par*n,lower.limits, upper.limits),ncol=Par)

fn_xi <- rep(0,n)  
v_xi <- matrix(0,nrow=n,ncol=Par)
F_xi <- matrix(0,nrow=n,ncol=Par)

GP <- array(NA,c(n,3*Par + 1,iterations))
GP[,,1] <- cbind(X,fn_xi,v_xi,F_xi)
dimnames(GP) <- list(c(paste(1:n)),Names)

GMemory <- matrix(NA,ncol=3*Par + 1,nrow=n)
colnames(GMemory) <- Names


# ---------------------------------------------------------------------------
# handling outliers
Handling <- function(GP,k){
  for (i in 1:Par){
		out.left <- which(GP[,i,k] < lower.limits)
		out.right <- which(GP[,i,k] > upper.limits)
		
		if(sum(out.left,out.right) >= 1){
			if(length(out.left) >= 1){
			GP[out.left,i,k] <<- runif(length(out.left),lower.limits, upper.limits)
			}
		if(length(out.right) >= 1){
			GP[out.right,i,k] <<- runif(length(out.right),lower.limits, upper.limits)
			} 
		}
	}	
}	



#----------------------------------------------------------------------------
# sorting solutions
SortSolutions <- function(GP,k){  
	sortet <- sort(GP[,"fn_xi",k],decreasing=TRUE,index.return=TRUE)$ix 
	GP[,,k] <<- GP[sortet,,k]
}


# ---------------------------------------------------------------------------
# gravity memory with m/n best solutions
GM <- function(GP,k){
   
	 if(k == 1){
	 GMemory <<- GP[1:m,,1]	 
	 } else {
	 
	Aux <- rbind(GMemory,GP[,,k])
	sortet <- sort(Aux[,"fn_xi"],decreasing=TRUE,index.return=TRUE)$ix[1:m]
	GMemory <<- Aux[sortet,]
	
	}	
}


#---------------------------------------------------------------------------- Forces start
# radius (separating global/local optimisation); search strategies 
Radius <- function(GP,k){
	if(man.scaling == TRUE){
		alpha*max(dist(GP[,1:Par,k]))
	} else {
		0.5*(1 - k/iterations)*max(dist(GP[,1:Par,k]))
	}
}


#----------------------------------------------------------------------------
# probabilities for attraction
p_ij <- function(GP,k,i,j){	
	if((GP[j,"fn_xi",k] > GMemory[i,"fn_xi"]) & (runif(1) > 0.1)){    
    	return(0)
    	}else {return(1)}
}


# ---------------------------------------------------------------------------
# calculate distances between particles 
r_ij <- function(GP,k,i,j){  
	result <- sqrt(sum((GMemory[i,1:Par] - GP[j,1:Par,k])^2))
	return(result)		
}


# resulting force on a particle j ---------------------------------------------
F_j <- function(GP,k){

a <- Radius(GP=GP,k)
summand <- matrix(NA,ncol=Par,nrow=m)

for (j in 1:n){
	for (i in 1:m){
  
  	r <- r_ij(GP,k,i,j) # distance between particles i and j
 
# gravitational fields ----------------------------------------------------------
		if (r < a){ # particle inside inner gravitational radius
			summand[i,] <- abs(GMemory[i,"fn_xi"])/{a*a*a}*(GMemory[i,1:Par] - GP[j,1:Par,k])*p_ij(GP,k,i,j) 
		} else{  	# particle inside outside gravitational radius
			summand[i,] <- abs(GMemory[i,"fn_xi"])/{r*r*r}*(GMemory[i,1:Par] - GP[j,1:Par,k])*p_ij(GP,k,i,j)
		}	
	}
	GP[j,c(paste("F_x",1:Par,sep="")),k] <<- apply(summand,2,sum) 
	}
}


#---------------------------------------------------------------------------- Forces end
# calculate new position ---------------------------------------------------------------
k_f <- function(k){0.5*(1 + k/iterations)} # Force, -> 1
k_v <- function(k){0.5*(1 - k/iterations)} # Velocity, -> 0

X_new <- function(GP,k){
  GP[,1:Par,k+1] <<- runif(n)*k_f(k)*GP[,c(paste("F_x",1:Par,sep="")),k] + runif(n)*k_v(k)*GP[,c(paste("v_x",1:Par,sep="")),k] + GP[,1:Par,k]
	}


# calculate new velocity --------------------------------------------------------------
V_new <- function(GP,k){
	GP[,c(paste("v_x",1:Par,sep="")),k+1] <<- GP[,1:Par,k+1] - GP[,1:Par,k]
	}


#----------------------------------------------------------------------------
# ---- calculations
pb <- txtProgressBar(min=1,max=iterations-1,style=3)
for (k in 1:(iterations-1)){
	
	setTxtProgressBar(pb, k)

 	Handling(GP=GP,k)
	GP[,"fn_xi",k] <- fn(X=GP[,1:Par,k])
  	SortSolutions(GP=GP,k)
	GM(GP=GP,k)
	F_j(GP=GP,k)
	X_new(GP=GP,k)
	V_new(GP=GP,k)
	
 	if(k == (iterations-1)){
 		cat("\nProcess finished.\n")
 		close(pb)
 	}
}

#------------------------------------------ end of function  ----------------------
return(list(GP = GP, Memory = GMemory[,1:(Par+1)]))
}

