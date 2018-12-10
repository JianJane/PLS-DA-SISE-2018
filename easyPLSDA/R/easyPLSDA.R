plsDA <- function(formula,data=NULL,ncomp=2,method="classic",auto.select.var=TRUE,threshold=0.8,tol=10^-9,scale=TRUE){

  instance <- list()
  extractedDF <- model.frame(formula,data)
  X <- extractedDF[,-1]
  var.names.X <- colnames(X)
  Y <- extractedDF[,1]

  #Test on rank of matrix X
  if(Matrix::rankMatrix(X)[1]<ncomp){
    stop("Rank of the matrix is lower than number of composantes")
  }
  isfactor <- is.factor(Y)
  if(isfactor){

    Y <- dummies::dummy(Y)
    instance$Y.dummy <- Y
  }else{
    stop("The response variable seems to not be a factor, please provide a categorical factor.")
  }



  #Scale X if needed
  if(scale==TRUE){
    Xp <- scale(X)
    instance$scaled <- TRUE

  }else{
    Xp <- X
    instance$scaled <- FALSE

  }
  computedPls <- pls2(Xp,Y,ncomp=ncomp,method=method,tol=tol)

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
    #Scale X if needed
    if(scale==TRUE){
      Xp <- scale(X)
      instance$scaled <- TRUE

    }else{
      Xp <- X
      instance$scaled <- FALSE

    }

    print(tol)

    computedPls <- pls2(Xp,Y,ncomp=ncomp,method=method,tol=)

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
