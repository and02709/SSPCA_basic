#' This function accepts a response vector and converts into either a linear or
#' a delta kernel
#' 
#' @param Y initial response vector
#' @param n an integer storing the number response observations
#' @param kernel stores they type of kernel we wish to build, either linear or 
#' delta
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 
  

respkernel <- function(Y, n, kernel){
  
  if (kernel=="linear"){
    L <- Rfast::Tcrossprod(Y,Y)
  } else if (kernel=="delta"){
    yf <- as.factor(Y)
    L<-matrix(0,n,n)
    for (i in levels(yf)){
      tmp<-yf==i
      L[tmp,tmp]<-1
    } 
  } else {
    stop("Please select a valid kernel, linear kernel or delta kernel")
  } 
  return(L)
}