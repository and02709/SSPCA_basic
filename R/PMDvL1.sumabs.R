#' This function peforms penalized matrix decomposition on a supplied matrix
#' and uses sum of absolute value of loadings as the shrinkage mechanism
#' @param X matrix to undergo PMD
#' @param npc number of desired principal components
#' @param n number of observations
#' @param p number of predictors
#' @param sumasv sum of absolute value of loadings of the v vector
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples

PMDvL1.sumabs <- function(X, npc, n, p, sumabsv, niter, trace){
  # Initialize the list of singular values and orthogonal matrices
  d <- rep(0, npc)
  u <- matrix(0, nrow=n, ncol=npc)
  v <- matrix(0, nrow=p, ncol=npc)
  
  # Initialize X matrix to be used
  Xuse <- X
  
  for(i in 1:npc){
    if(sumabsv==sqrt(p)){
      temp <- RSpectra::svds(Xuse, k=1)
      d[i] <- temp$d
      u[,i] <- temp$u
      v[,i] <- temp$v
    }
    else{
      temp <- RSpectra::svds(Xuse, k=1)
      if(i==1){
        temp <- SSPCA::SMD.sumabs(X=Xuse, d=temp$d, u=temp$u, v=temp$v, n=n, p=p, 
                                  sumabsv=sumabsv, niter=niter, trace=trace,
                                  npc=i)
      }
      else{
        temp <- SSPCA::SMDorth.sumabs(X=Xuse, us=u[,1:(i-1)], d=temp$d, u=temp$u, v=temp$v, n=n, p=p, 
                                      sumabsv=sumabsv, niter=niter, trace=trace,
                                      npc=i)
      }
      d[i] <- temp$d
      u[,i] <- temp$u
      v[,i] <- temp$v
    }
    if(i==npc) break
    Xuse <- Xuse - d[1]*u[,1]%*%t(v[,1])
  }
  
  return(list(d=d,u=u,v=v))
  
}