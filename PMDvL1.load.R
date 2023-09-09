#' This function peforms penalized matrix decomposition on a supplied matrix
#' and uses nonzero loadings as the shrinkage mechanism
#' @param X matrix to undergo PMD
#' @param npc number of desired principal components
#' @param n number of observations
#' @param p number of predictors
#' @param nonzero.loadings desired number of nonzero loadings
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples

PMDvL1.load <- function(X, npc, n, p, nonzero.loadings, niter, trace){
  # Initialize the list of singular values and orthogonal matrices
  d <- rep(0, npc)
  u <- matrix(0, nrow=n, ncol=npc)
  v <- matrix(0, nrow=p, ncol=npc)
  
  # Initialize X matrix to be used
  Xuse <- X
  
  for(i in 1:npc){
    if(nonzero.loadings==p){
      temp <- RSpectra::svds(Xuse, k=1)
      d[i] <- temp$d
      u[,i] <- temp$u
      v[,i] <- temp$v
    }
    else{
      temp <- RSpectra::svds(Xuse, k=1)
      if(i==1){
        temp <- SSPCA::SMD.load(X=Xuse, d=temp$d, u=temp$u, v=temp$v, n=n, p=p, 
                                nonzero.loadings=nonzero.loadings, niter=niter, trace=trace,
                                npc=i)
      }
      else{
        temp <- SSPCA::SMDorth.load(X=Xuse, us=u[,1:(i-1)], d=temp$d, u=temp$u, v=temp$v, n=n, p=p, 
                                    nonzero.loadings=nonzero.loadings, niter=niter, trace=trace,
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
