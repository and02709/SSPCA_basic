#' This function peforms penalized matrix decomposition for one pair of
#' u and v vectors.  This uses orthogonilization for the previous U vectors.
#' Uses sum of absolute value of loadings as the shrinkage mechanism
#' @param X matrix to undergo PMD
#' @param d supplied singular value
#' @param u supplied left vector
#' @param v supplied right vector
#' @param n number of observations
#' @param p number of predictors
#' @param sumasv sum of absolute value of loadings of the v vector
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @param npc number of desired principal components
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples

SMDorth.sumabs <- function(X, us, d, u, v, n, p, sumabsv, niter, trace, npc){
  oldv <- rnorm(p,0,1)
  oldu <- rnorm(n,0,1)
  cat("Vector ", npc, ": ")
  for(iter in 1:niter){
    if((sum(abs(oldv-v)) < 1e-7) && (sum(abs(oldu-u)) < 1e-7)) break
    oldv <- v
    oldu <- u
    if(trace) cat(iter," ", fill=F)
    # update u
    argu <- Rfast::mat.mult(X,v)
    numer <- lsfit(y=argu, x=us, intercept=FALSE)$res
    u <- matrix(numer/l2n(numer),ncol=1)
    # update v
    argv <- Rfast::Crossprod(u,X)
    v <- matrix(argv/l2n(argv),ncol=1)
    # Find appropriate shrinkage 
    lamv <- BinarySearch(argu=argv,sumabs=sumabsv)
    # soft threshold v
    sv <- soft(argv,lamv)
    v <- matrix(sv/l2n(sv),ncol=1)
  }
  cat("\n")
  d <- as.numeric(t(u)%*%(X%*%v))
  return(list(d=d, u=u, v=v))
}