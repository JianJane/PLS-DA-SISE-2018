plsDA <- function(X,Y,ncomp=2,method="classic",auto.select.var=TRUE,threshold=0.8,tol=10^-9){

  instance <- list()

  var.names.X <- colnames(X)
  #Test on rank of matrix X
  if(min(Matrix::rankMatrix(X)[1],nrow(X))<ncomp){
    stop("Rank of the matrix is lower than number of composantes")
  }
  computedPls <- pls2(X,Y,ncomp=ncomp,method=method,tol=tol)

  #Computation of explained variance (taken from https://github.com/gastonstat/):
  R2.X <- colMeans(cor(X,computedPls$scores$X)^2)
  R2.Y <- colMeans(cor(Y,computedPls$scores$X)^2)
  Rd.mat = matrix(0,ncomp, ncomp)
  for (j in 1:ncomp){
    Rd.mat[1:j,j] = R2.Y[1:j]
  }
  # variable importance computation
  VIP = sqrt((computedPls$weights$X^2) %*% Rd.mat %*% diag(ncol(X)/cumsum(R2.Y), ncomp, ncomp))

  ### Weighted VIP for entire fitted model
  weighted.vip <- matrix(0, nrow = ncol(X), ncol=ncomp)
  for (i in 1:length(R2.Y)){
    weighted.vip[,i] <- R2.Y[i] * (computedPls$weights$X[,i]^2)
  }
  VIP.weighted <- sqrt(rowSums(weighted.vip) * (ncol(X)/sum(R2.Y)))

  #Autoslection of variables if needed
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
    # variable importance computation
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
  instance$Y <- Y
  instance$levels <- levels(Y)
  instance <- append(instance,computedPls)
  instance$explained.var$X <- R2.X
  instance$explained.var$Y <- R2.Y
  instance$VIP <- VIP

  class(instance) <- "PLSDA"

  return(instance)
}
