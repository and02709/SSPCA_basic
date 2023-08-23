#' This function is the parent function for performing SPCA
#' @param X nxp matrix of predictors
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param kernel designation for either linear or delta kernel
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

SPCA <- function(X, Y, npc, kernel=c("linear", "delta")){
  # Convert X to a matrix
  X <- as.matrix(X)
  # Convert Y to a matrix
  Y <- as.matrix(Y)
  # Stop if number of observations doesn't match between X and Y
  if(dim(X)[1] != dim(Y)[1]) stop("Number of observations do not match")
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  
  # Generate response kernel
  L <- SSPCA::respkernel(Y=Y, n=n, kernel=kernel)
  
  if(p > n){
    return(SSPCA::dual.decomp(X=X,L=L, npc=npc, n=n, p=p))
  }
  else{
    return(SSPCA::reg.decomp(X=X, L=L, npc=npc, n=n, p=p))
  }
}
