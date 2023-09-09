#' cross validation function for sparse supervised PCA
#' @param arg.sparse vector for folds and sparse parameter
#' @param df.partition list of folds
#' @param npc number of desired principal components
#' @param n.folds number of folds to perform cross validation
#' @param sparsity.type specifies mechanism for vector shrinkage
#' @param nonzero.loadings specifies the desired number of nonzero loadings
#' @param sumabsv shrinkage using sum of absolute value of loadings
#' @param kernel specification of the response kernel
#' @param niter number of iterations in the SMD algorithm
#' @param trace display algorithm progress
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples 

cv.partition.SSPCA <- function(arg.sparse, df.partition, npc, n.folds, sparsity.type=c("loadings", "sumabs"),
                               nonzero.loadings=NULL,sumabsv=NULL,kernel=c("linear", "delta"), niter,trace){
  test.index <- arg.sparse[[1]]
  if(sparsity.type=="loadings") nonzero.loadings <- arg.sparse[[2]]
  else sumabsv <- arg.sparse[[2]]
  
  xtrain <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(-c(y,.folds))
  xtest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(-c(y,.folds))
  ytrain <- df.partition |> dplyr::filter(.folds!=test.index) |> dplyr::ungroup() |> dplyr::select(y)
  ytest <- df.partition |> dplyr::filter(.folds==test.index) |> dplyr::ungroup() |> dplyr::select(y)
  
  xmeans <- colMeans(xtrain)
  xtrain <- t(apply(xtrain, 1, function(x) x-xmeans))
  xtest <- t(apply(xtest, 1, function(x) x-xmeans))
  
  
  if(sparsity.type=="loadings") sspca.out <- SSPCA::SSPCA(X=xtrain,Y=ytrain,npc=npc,
                                                          sparsity.type=sparsity.type,nonzero.loadings=nonzero.loadings,
                                                          sumabsv=NULL,kernel=kernel,niter=niter,trace=trace)
  else sspca.out <- SSPCA::SSPCA(X=xtrain,Y=ytrain,npc=npc,
                                 sparsity.type=sparsity.type,nonzero.loadings=NULL,
                                 sumabsv=sumabsv,kernel=kernel,niter=niter,trace=trace)
  V <- sspca.out$v
  
  mod <- SSPCA::model.build.SSPCA(Y=as.matrix(ytrain),X=as.matrix(xtrain),V=V,kernel=kernel)
  metric <- SSPCA::model.metric.SSPCA(Y=as.matrix(ytest),X=as.matrix(xtest),V=V,mod=mod,kernel=kernel)
  return(metric)
}