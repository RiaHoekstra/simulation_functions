genGVAR <- function(nNode, propPositiveCon = 0.8, mean = 0.3, sd = 0.1, ...){
  
  # Generate contemporaneous network using GGM:
  omega <- bootnet::genGGM(nNode, propPositive = propPositiveCon, ...) # partial correlation matrix
  
  # Get covariance matrix sigma: 
  # sigma <- cov2cor(solve(diag(ncol(omega)) - omega)) 
  
  # Transpose from partial correlation matrix to a precision matrix:
  n <- ncol(omega)
  kappa <- solve(cov2cor(solve(diag(n) - omega)))
  kappa[abs(omega) < sqrt(.Machine$double.eps) & !diag(n)==1] <- 0
  
  # Generate temporal network:
  beta <- diag(1, nNode) 
  
  for (i in 1:nNode){
    beta[(i+(1-1)) %% nNode+1, i] <- sample(c(-1, 1), 1, 0.5)
  }
  
  beta <- beta * rnorm(nNode^2, mean, sd)
  
  # Transpose beta to PDC: 
  labmda <- solve(kappa)
  PDC <- t(beta / sqrt(diag(labmda) %o% diag(kappa) + beta^2))
  
  # Check eigenvalues:
  if(any(eigen(kappa)$values < 0)){
    warning("Precision matrix is not positive semi definite")
  }
  
  # Return:
  return(list(kappa = kappa,
              PCC = omega,
              beta = beta, 
              PDC = PDC))
}
