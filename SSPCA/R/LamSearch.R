#' This function peforms identifies the desired lambda penalty resulting
#' in the desired number of nonzero loadings in the v matrix
#' @param v vector to undergo shrinkage
#' @param nzl number of desired nonzero loadings
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples

LamSearch <- function(v,nzl){
  vorder <- v[order(abs(v), decreasing=T)]
  v1 <- abs(vorder[nzl])
  v2 <- abs(vorder[nzl+1])
  return((v1+v2)/2)
}