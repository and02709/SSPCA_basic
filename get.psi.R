#' This function accepts a data matrix and the response kernel matrix
#' and generates an nxp matrix that has supervision
#' 
#' @param X nxp matrix of predictors
#' @param L nxn response kernel
#' @param n an integer storing the number response observations
#' @param p an integer storing the number of predictors 
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

get.psi <- function(X,L,n,p){
  # Centering matrix
  H <- diag(x=1, nrow=n) - 1/n*rep(1,n)%*%t(rep(1,n))
  # Eigendecomposition of L
  Eigendecomp <- eigen(L)
  # Need to retain eigenvectors to reconstruct Delta
  U <- Eigendecomp$vectors
  # Eigenvalues of kernel response matrix L
  EV <- Eigendecomp$values
  # Generation of diagonal matrix of square root eigenvalues of L
  Sigmat <- diag(sqrt(zapsmall(EV)))
  # Generation of Delta such that t(Delta)%*%Delta=L
  Delta <- Rfast::Tcrossprod(Rfast::mat.mult(U, Sigmat), U)
  # Generation of Psi matrix
  return(Rfast::Crossprod(Delta, Rfast::mat.mult(H,X)))
}