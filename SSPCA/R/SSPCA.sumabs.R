#' This function performs sparse supervised PCA using sum of absolute values
#' of loadings for the v vector as the shrinkage mechanism
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param sumabsv shrinkage using sum of absolute value of loadings
#' @param kernel specification of the response kernel
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

SSPCA.sumabs <- function(X, Y, npc, 
                         sumabsv=NULL, 
                         kernel=c("linear", "delta"), niter=50, trace=F){
  
  # Find number of observations
  n <- dim(X)[1]
  # Find number of predictors
  p <- dim(X)[2]
  
  # Check nonzero.loadings
  if(is.null(sumabsv)) sumabsv <- sqrt(p)
  if(sumabsv < 1) stop("Please provide a number of nonzero loadings between 1 and sqrt(p)")
  if(sumabsv > sqrt(p)) stop("Please provide a number of nonzero loadings between 1 and sqrt(p)")
  
  # Check names of predictors
  if(is.null(colnames(X))) colnames(X) <- paste0("X",seq(1:p))
  
  # Generate response kernel
  L <- respkernel(Y=Y, n=n, kernel=kernel)
  # Calculate Psi matrix
  Psi <- get.psi(X=X, L=L, n=n, p=p)
  # Perform Penalized matrix decomposition (PMD)
  if(sumabsv==sqrt(p)) return(SPCA.svd(X=Psi, npc=npc))
  out <- PMDvL1.sumabs(X=Psi, npc=npc, n=n, p=p, sumabsv=sumabsv, niter=niter, trace=trace)
  
  return(out)
}