
#' plsDA function is an internal function of package easyPLSDA used to perform PLS-DA.
#'
#' @param X A matrix containing the explanatory variables of the model.
#' @param Y A dummy matrix of the response classes corresponding to the X observations.
#' @param ncomp The required number of components to compute.
#' @param method The method that will be used to compute the different matrices used in regression.
#'  Default "classic" is for PLS2 algorithm, "SIMPLS" is the other method based on matrix SVD.
#' @param auto.select.var Should the function automatically select variables using VIP model ? Default: True.
#' @param threshold The threshold for selecting the important variable.
#'  Variables with weighted VIP under this threshold will be removed if auto.select.var is TRUE. Default: 0.8.
#' @param tol The convergence threshold for latent scores computation, only for "classic" method. Default: 10^-9.
#'
#' @return S3 object of class plsDA.
#' \item{selected.var}{Names of explanatory variables selected.}
#' \item{X}{Matrix of explanatory variables.}
#' \item{scores}{Scores matrices of X and Y variables projected into the new space.}
#' \item{weights}{Matrices of weights for both X and Y variables.}
#' \item{loadings}{One (or two for SIMPLS method) matrix of loadings.}
#' \item{mode}{The mode used to compute the matrices.}
#' \item{explained.var}{The explained variance by latent vectors of X and Y.}
#' \item{VIP}{Matrix of Variable Importance in Projection.}
#'
#'
#'
plsDA <- function(X,Y,ncomp=2,method="classic",auto.select.var=FALSE,threshold=0.8,tol=10^-9){

  instance <- list()
  var.names.X <- colnames(X)

  #Test on rank of matrix X
  if(min(Matrix::rankMatrix(X)[1],nrow(X))<ncomp){
    stop("Rank of the matrix is lower than asked number of components")
  }
  computedPls <- pls2(X,Y,ncomp=ncomp,method=method,tol=tol)

  #Computation of explained variance and VIP (taken from https://github.com/gastonstat/):
  R2.X <- colMeans(cor(X,computedPls$scores$X)^2)
  R2.Y <- colMeans(cor(Y,computedPls$scores$X)^2)
  Rd.mat = matrix(0,ncomp, ncomp)
  for (j in 1:ncomp){
    Rd.mat[1:j,j] = R2.Y[1:j]
  }

  #VIP computation
  VIP = sqrt((computedPls$weights$X^2) %*% Rd.mat %*% diag(ncol(X)/cumsum(R2.Y), ncomp, ncomp))

  ### Weighted VIP for entire fitted model
  weighted.vip <- matrix(0, nrow = ncol(X), ncol=ncomp)
  for (i in 1:length(R2.Y)){
    weighted.vip[,i] <- R2.Y[i] * (computedPls$weights$X[,i]^2)
  }
  VIP.weighted <- sqrt(rowSums(weighted.vip) * (ncol(X)/sum(R2.Y)))

  #Autoslection of variables if asked
  if(auto.select.var==T && sum(VIP.weighted<threshold)>0){
    var.to.exclude <- VIP.weighted<threshold
    instance$selected.var <- var.names.X[VIP.weighted>=threshold]
    X <- X[,colnames(X)[!var.to.exclude]]

    var.names.X <- colnames(X)


    computedPls <- pls2(X,Y,ncomp=ncomp,method=method,tol=tol)

    #Computation of explained variance (taken from https://github.com/gastonstat/):
    R2.X <- colMeans(cor(X,computedPls$scores$X)^2)
    R2.Y <- colMeans(cor(Y,computedPls$scores$X)^2)
    Rd.mat = matrix(0,ncomp, ncomp)
    for (j in 1:ncomp){
      Rd.mat[1:j,j] = R2.Y[1:j]
    }
    #VIP computation
    VIP  <- sqrt((computedPls$weights$X^2) %*% Rd.mat %*% diag(ncol(X)/cumsum(R2.Y), ncomp, ncomp))

    ### Weighted VIP for entire fitted model
    weighted.vip <- matrix(0, nrow = ncol(X), ncol=ncomp)


  }else{
    instance$selected.var <- var.names.X
  }

  # combine weighted VIPs with individual component VIPs
  for (i in 1:length(R2.Y)){
    weighted.vip[,i] <- R2.Y[i] * (computedPls$weights$X[,i]^2)
  }
  VIP.weighted <- sqrt(rowSums(weighted.vip) * (ncol(X)/sum(R2.Y)))
  VIP <- cbind(VIP, VIP.weighted)
  dimnames(VIP) = list(var.names.X, c(paste(rep("Component ",ncomp),1:ncomp,sep=""),"Model VIP"))




  #Creation of output
  instance$X <- X
  instance <- append(instance,computedPls)
  instance$explained.var$X <- R2.X
  instance$explained.var$Y <- R2.Y
  instance$VIP <- VIP

  class(instance) <- "PLSDA"

  return(instance)
}
