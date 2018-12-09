library(Matrix)
library(dummies)
plsDA <- function(formula,data=NULL,ncomp=2,method="classic",tol=10^-9,scale=TRUE){

  instance <- list()
  extractedDF <- model.frame(formula,data)
  X <- extractedDF[,-1]
  Y <- extractedDF[,1]
  instance$X <- X
  instance$Y <- Y
  instance$levels <- levels(Y)
  #Test on rank of matrix X
  if(rankMatrix(X)[1]<ncomp){
    stop("Rank of the matrix is lower than number of composantes")
  }
  isfactor <- is.factor(Y)
  if(isfactor){

    Y <- dummy(Y)
    instance$Y.dummy <- Y
  }



  #
  if(scale==TRUE){
    X=scale(X)
    if(!isfactor){
      Y=scale(Y)
    }

  }

  computedPls <- pls2(X,Y,ncomp=ncomp,method=method,tol=tol)

  instance <- append(instance,computedPls)

  class(instance) <- "PLSDA"

  return(instance)
}
