#' Parent function or sparse supervised PCA
#' @param X nxp predictor matrix
#' @param Y nx1 response vector
#' @param npc number of desired principal components
#' @param sparsity.type specifies mechanism for vector shrinkage
#' @param nonzero.loadings specifies the desired number of nonzero loadings
#' @param sumabsv shrinkage using sum of absolute value of loadings
#' @param kernel specification of the response kernel
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

SSPCA <- function(X, Y, npc, sparsity.type=c("loadings", "sumabs"), 
                  nonzero.loadings=NULL, sumabsv=NULL, 
                  kernel=c("linear", "delta"), niter=50, trace=F){
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
  
  if(sparsity.type != "loadings" && sparsity.type!= "sumabs") stop("Must specify correct sparse penalty")
  if(sparsity.type=="loadings"){
    return(SSPCA.load(X=X,Y=Y,npc=npc,nonzero.loadings = nonzero.loadings,
                      kernel=kernel, niter=niter,trace=trace))
  } 
  if(sparsity.type=="sumabs"){
    return(SSPCA.sumabs(X=X,Y=Y,npc=npc,sumabsv=sumabsv,
                        kernel=kernel, niter=niter,trace=trace))
  }
  
  stop("SOMETHING GONE WRONG")
  
}