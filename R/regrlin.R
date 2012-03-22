
### implementation by Gianluca Bontempi

`.regrlin` <- function(X,Y,X.ts=NULL,lambda=1e-3){
#require(MASS)

	n <- NCOL(X) # number input variables
	p <- n+1
	N <- NROW(X) # number training data 
  
	XX <- cbind(array(1,c(N,1)),X)

	H1 <- MASS::ginv(t(XX)%*%XX+lambda*diag(p))
	beta.hat <- H1%*%t(XX)%*%Y
	H <- XX%*%H1%*%t(XX)
	Y.hat <- XX%*%beta.hat
	e <- Y-Y.hat
	var.hat.w <- (t(e)%*%e)/(N-p)
	e.loo <- e/(1-diag(H))

  return(e.loo)
}
