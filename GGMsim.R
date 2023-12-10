GGMsim <- function(n, PCC, missing = 0){
  
  # check if matrix is positive semi-definite
  if (any(eigen(diag(ncol(PCC)) - PCC)$values < 0)){
    stop("matrix is not positive semi-definite") 
  }
  
  # Get covariance matrix sigma: 
  sigma <- cov2cor(solve(diag(ncol(PCC)) - PCC))  
  
  # Generate data:
  data <- mvtnorm::rmvnorm(n, sigma = sigma)
  
  # Add missing:
  if (missing > 0){
    for (i in 1:ncol(data)){
      data[runif(nrow(data)) < missing, i] <- NA
    }
  }
  
  return(data)
}

