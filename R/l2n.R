#' This function calculates the l2 norm of a vector
#' @param vec input vector
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  return(a)
}