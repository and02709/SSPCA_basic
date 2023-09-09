#' This function performs soft thresholding on a vector.
#' @param x vector to undergo soft thresholding
#' @param d value for absolute value loading shrinkage
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}