#' easy Partial Least Square Discriminant Analysis class constructor function
#'
#' easyPLSDA class constructor is to generate a class objet that performs PLS-2 Discriminant analysis on subject data sets.
#' It is a class object contains functions for different aspects of analysis and a class method for prediction.
#' All of the analysis is performed upon instantiating the class constructor along with the target data set, the results returned as lists.
#' The prediction is realised separately by calling the method function 'predict', results are returned in a list.
#'
#' Aspects of analysis imbeded in the constructor include:
#' Note: Two methods are made available for this extraction step. The iterative method based on NIPALS PLS2 is refered to as the 'classic' method.
#' The "SIMPLS" method uses SVD decomposition function in R.
#' 1.The extraction of latent score matrices of X and Y, refered to as T.scores and U.scores.
#' Their respective loadings matrices, refered to as P.loadings and Q.loadings. Together with two sets of column vectors,
#' named "X.weights" and "Y.weights" that are needed to linearly combine the columns of X and Y into their respective latent scores.
#' 2.The optimal number of scores or components analysis is performed by tuning the accumulated variance threshold, see 'threshold.comp' in Arguments.
#' 3.Variable importance in projection (VIP) scores are calculated and an option of variable auto selection based on VIP is availabe, see 'auto.select.var' in Arguments.
#'
#'
#' @usage easyPLSDA(formula,data=NULL,ncomp=2,method="classic",auto.select.var=TRUE,threshold=0.8,threshold.comp=0.95,maxi.comp=10,tol=10^-9,scale=TRUE)
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data data = NULL by default, a data matrix or data frame object containing observations data matrix X and response data matrix Y
#' @param ncomp ncomp=2 by default, denoting the number of principle components to be returned by the function
#' @param method method= "classic' by default, "SIMPLS" is the other method based on matrix diagnolisation SVD
#' @param auto.select.var Option for automatically selecting the variables explaining the most amount of accumulated variance in response matrix Y
#' @param threshold threshold = 0.8 by default, high pass threshold for selecting important variables in terms of their explained variance in response matrix Y
#' @param threshold.comp threshold.comp = 0.95 by default high pass threshold for selecting principal components explaining most of the variance
#' @param maxi.comp maxi.comp = 10 by default the maximum number of scores to be included in cross validation for finding the optimal number of scores
#' @param tol tol=10^-9 by default, as the convergence threshold for latent scores itteration
#' @param scale scale=TRUE by default, after being column centered, X and Y can be further normalised with the standard deviation of each column
#'
#'
#'
#'@details The 'predict' method function of the class predicts the response matrix from an unseen set of data
#'
#'
#'@return easyPLSDA returns an object of class "easyPLSDA"
#'The function summary is used to print a summary of the analysis results.
#'An object of class 'easyPLSDA' returns a list containing the following components:
#'
#' \item{scores}{Scores matrices of X and Y variables projected into the new space}
#' \item{weights}{Matrices of weights for both X and Y variables}
#' \item{loadings}{One (or two for SIMPLS method) matrix of loadings}
#' \item{B.mat}{Only for "classic" mode. Matrix of "b" coefficients, result of the scalar product between X latent vectors and Y latent vectors}
#' \item{mode}{The mode used to compute the matrices}
#' \item{selected.var}{Names of explanatory variables selected}
#' \item{X}{Matrix of explanatory variables}
#' \item{scores}{Scores matrices of X and Y variables projected into the new space}
#' \item{weights}{Matrices of weights for both X and Y variables}
#' \item{loadings}{One (or two for SIMPLS method) matrix of loadings}
#' \item{mode}{The mode used to compute the matrices}
#' \item{explained.var}{The explained variance by latent vectors of X and Y}
#' \item{VIP}{Matrix of Variable Importance in Projection}
#'
#'
#'
#'
#'@example
#'## Wine evaluation of a set of 5 wines.## Features: price, sugar, alcohol, acidity.
#' ##Rating features: likeability, compatible with meat, compatible with dessert.
#'X_wine<-matrix(c(7, 7, 13, 7, 4, 3, 14, 7, 10, 5, 12, 5, 16, 7, 11, 3, 13, 3, 10, 3 ),5,4, byrow=TRUE)
#'Y_wine<-matrix(c(14, 7, 8, 10, 7, 6, 8, 5, 5, 2, 4, 7, 6, 2, 4 ),5,3, byrow=TRUE)
#'mX<-as.matrix(X_wine)
#'mY<-as.matrix(Y_wine)
#'
#'obj<-plsDA(formula,data=NULL,ncomp=2,method="classic",tol=10^-9,scale=TRUE)
#'fit.results <- fit(obj)
#'predict.results<- predict(obj)
#'vip.results<-vip(obj)
#'
#'
#'@import Matrix
#'@import dummies

easyPLSDA <- function(formula,data=NULL,ncomp=2,method="classic",auto.select.var=FALSE,threshold=0.8,threshold.comp=0.95,maxi.comp=10,tol=10^-9,scale=TRUE){


  instance <- list()

  #Extraction of the formula variables in two part (explanatory varaibles X and response Y)
  extractedDF <- model.frame(formula,data)
  X <- extractedDF[,-1]
  Y <- extractedDF[,1]

  #Verification step
  isfactor <- is.factor(Y)
  if(isfactor){
    instance$levels <- levels(Y)
    Y <- dummies::dummy(Y)
    colnames(Y) <- instance$levels
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


  if(is.null(ncomp)){
    #Automatic selection by Leave One Out cross validation

    #maxi.comp (10) model by default to avoid Out of Memory Exception
    max.comp <- min(Matrix::rankMatrix(as.matrix(X))[1],nrow(X),maxi.comp)

    press.test <- vector(length = max.comp)
    press = matrix(nrow=nrow(X), ncol=ncol(Y))

    #For the first component:
    for (i in 1:nrow(X)) {
      computedPLSDA <- plsDA(Xp[-i,],Y[-i,],ncomp=1,method=method,auto.select.var=F,threshold=threshold,tol=tol)
      test <- Xp[i,]
      dim(test) <- c(1,length(Xp[i,]))
      colnames(test) <- colnames(X)

      prediction <- predict(computedPLSDA,test,scale=FALSE)
      press[i,] <- (Y[i,]-prediction$pred)^2
    }

    press.test[1] <- sum(press)

    #For the next components:
    for (k in 2:max.comp) {
      press = matrix(nrow=nrow(X), ncol=ncol(Y))
      for (i in 1:nrow(X)) {
        computedPLSDA <- plsDA(Xp[-i,],Y[-i,],ncomp=k,method=method,auto.select.var=F,threshold=threshold,tol=tol)
        test <- Xp[i,]
        dim(test) <- c(1,length(Xp[i,]))
        colnames(test) <- colnames(X)
        prediction <- predict(computedPLSDA,test,scale=FALSE)
        press[i,] <- (Y[i,]-prediction$pred)^2
      }

      print(k)
      press.test[k] <- sum(press)


    }

    #We use the Rk statistic to select the optimal number of components
    Rk <- vector(length = max.comp-1)
    for (p in 2:length(press.test)) {
      Rk[p-1] <- press.test[p]/press.test[p-1]

    }
    ncomp.selected <- which(Rk>threshold.comp)[1]+1
    if(auto.select.var==T){
      computedPLSDAfin <- plsDA(Xp,Y,ncomp=ncomp.selected,method=method,auto.select.var=T,threshold=threshold,tol=tol)
      instance <- append(instance,computedPLSDAfin)
      instance$comp.selected <- ncomp.selected
      instance$Rk <- Rk
    }else{
      computedPLSDAfin <- plsDA(Xp,Y,ncomp=ncomp.selected,method=method,auto.select.var=F,threshold=threshold,tol=tol)
      instance <- append(instance,computedPLSDAfin)
      instance$comp.selected <- ncomp.selected
      instance$Rk <- Rk
    }

  }else{

    #For manual selection of components number
    instance <- append(instance,plsDA(Xp,Y,ncomp=ncomp,method=method,auto.select.var=auto.select.var,threshold=threshold,tol=tol))
    instance$comp.selected <- ncomp
  }

  class(instance) <- "PLSDA"

  return(instance)
}
