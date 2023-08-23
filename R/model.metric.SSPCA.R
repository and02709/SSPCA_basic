#' cross validation function for sparse supervised PCA
#' @param Y vector of response 
#' @param X predictor matrix
#' @param V pca vector
#' @param mod holds the learned model
#' @param kernel flag for kernel type
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

model.metric.SSPCA <- function(Y,X,V,mod,kernel=c("linear","delta")){
  Z <- X%*%V
  df <- data.frame(y=Y,z=Z)
  n <- nrow(df)
  if(kernel=="linear"){
    yhat <- predict(mod, df)
    metric <- sqrt(sum((Y-yhat)^2)/n)
  }
  else{
    yhat <- predict(mod,df) |> as.matrix()
    measure.yhat <- (Y==yhat)
    acc <- mean(measure.yhat)
    metric <- 1-acc
  }
  return(metric)
}