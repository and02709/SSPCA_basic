#' This function performs supervised PCA on the psi matrix
#' @param X nxp matrix to undergo supervised PCA
#' @param npc number of desired principal components
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

SPCA.svd <- function(X, npc){
  temp <- RSpectra::svds(X, npc)
  return(list(d=temp$d, u=temp$u, v=temp$v))
}