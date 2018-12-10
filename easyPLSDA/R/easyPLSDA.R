#' easy Partial Least Square Discriminant Analysis class generating function
#'
#' plsDA is used to generate a class objet based on PLS-2 Discriminant analysis algorithms.
#' Its class objects can be used to carry out PLS discriminant analysis, to predict response matrix on unseen data,
#' to calculate the variable importance in projection "VIP", to optimise the prediction by filtering out the less relevant principal components
#' to further optimise the prediction by performing VIP or Variable Importance in Projection analysis
#'
#' @usage plsDA(formula, data=NULL, ncomp=2, method="classic", tol=10^-9, scale=TRUE)
#'
#'
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param ncomp= or ncomp=2 by default, denoting the number of principle components to be returned by the function
#' @param method= or method= "classic' by default, "SIMPLS" is the other method based on matrix diagnolisation SVD
#' @param tol= or tol=10^-9 by default, as the convergence threshold for latent scores itteration
#' @param scale or scale=TRUE by default, after being column centered, X and Y can be further normalised with the standard deviation of each column
#'
#'
#'@details The 'pls2' method function of the class extracts the numerical and geometric features from the independet variable matrix X and its linearly dependent response matrix Y.
#' Features of interest include the latent scores "T.scores" and "U.scores" from X and Y respectively, their respective loadings vectors "P.loadings" and "Q.loadings", along with two sets of column vectors,
#' named "X.weights" and "Y.weights" that are needed to linearly combine the columns of X and Y into their respective latent scores. The itterative method is set as default. Both methods return the features in a list.
#'
#'@details The 'predict' method function of the class predicts the response matrix from an unseen set of data
#'
#'
#'@import Matrix
#'@import dummies
#
#'@example
#'


#library(Matrix)   commented out by J.BI 11/Dec/2018
#library(dummies)  commented out by J.BI 11/Dec/2018

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
