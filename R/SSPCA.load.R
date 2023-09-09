#' This function performs sparse supervised PCA using number of nonzero
#' loadings for the shrinkage mechanism
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param nonzero.loadings desired number of nonzero loadings
#' @param kernel specification of the response kernel
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

SSPCA.load <- function(X, Y, npc, 
                       nonzero.loadings=NULL, 
                       kernel=c("linear", "delta"), niter=50, trace=F){
  
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  
  # Check nonzero.loadings
  if(is.null(nonzero.loadings)) nonzero.loadings <- p
  if(nonzero.loadings < 1) stop("Please provide a number of nonzero loadings between 1 and p")
  if(nonzero.loadings > p) stop("Please provide a number of nonzero loadings between 1 and p")
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  
  # Generate response kernel
  L <- SSPCA::respkernel(Y=Y, n=n, kernel=kernel)
  # Calculate Psi matrix
  Psi <- SSPCA::get.psi(X=X, L=L, n=n, p=p)
  # Perform Penalized matrix decomposition (PMD)
  if(nonzero.loadings==p) return(SPCA.svd(X=Psi, npc=npc))
  out <- SSPCA::PMDvL1.load(X=Psi, npc=npc, n=n, p=p, nonzero.loadings=nonzero.loadings, niter=niter, trace=trace)
  
  return(out)
}