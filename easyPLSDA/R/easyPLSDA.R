#'
#' easy Partial Least Square Discriminant Analysis class constructor function
#'
#' easyPLSDA class constructor is used to generate a class object to perform PLS-2 Discriminant analysis on subject data sets.
#' All of the analysis is performed upon instantiating the class constructor along with the target data set, the results are returned as a list.
#' The prediction is realised separately by calling the function 'predict', results are returned in a list. A summary method is also available.
#'
#' @usage easyPLSDA(formula,
#' data=NULL,
#' ncomp=2,
#' method="classic",
#' auto.select.var=FALSE,
#' threshold=0.8,
#' threshold.comp=0.95,
#' maxi.models=10,
#' tol=10^-9,
#' scale=TRUE)
#'
#'
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Similar to lm function.
#' @param data A data matrix or data frame object containing observations data matrix X and response data matrix Y, NULL by default.
#' @param ncomp Positive non zero integer, denoting the number of principale components to be returned by the function, 2 by default.
#' If NULL, the function will compute the optimal number of components by Leave One Out Cross Validation and using the Rk Criterium.
#' @param method Charactor string, "classic' by default, "SIMPLS" is the other method based on matrix SVD.
#' @param auto.select.var Logical. Option to enable automatic selection of the variables explaining the most amount of accumulated variance. FALSE by default.
#' @param threshold Numeric, high pass threshold for selecting important variables in terms of their explained variance in response matrix Y.
#'  0.8 by default.
#' @param threshold.comp Numeric, Threshold used for R of Wold criterium in automatic selection of components.
#' 0.95 by default, corresponding to the adjusted R criterium.
#' @param maxi.models Numeric. The maximum number of models to be tested in cross validation for finding the optimal number of components. 10 by default.
#' If the train dataset allows less than 10 components, then the function test only the available number of components.
#' @param tol  Numeric, the convergence threshold for latent scores computation, only for "classic" method. Default: 10^-9.
#' @param scale Logical,if TRUE (the default) after being column centered, X and Y can be further normalised with the standard deviation of each column.
#'
#'
#'@details The 'pls2' method function of the class extracts the numerical and geometric features from the independet variable matrix X and its linearly dependent response matrix Y.
#' Features of interest include the latent scores "T.scores" and "U.scores" from X and Y respectively, their respective loadings vectors "P.loadings" and "Q.loadings", along with two sets of column vectors,
#' named "X.weights" and "Y.weights" that are needed to linearly combine the columns of X and Y into their respective latent scores. The itterative method is set as default. Both methods return the features in a list.
#' Aspects of analysis included in the returned object include:
#' 1.The extraction of latent score matrices of X and Y, refered to as scores$X and scores$Y.
#' 2.Their respective loadings matrices, refered to as P.loadings (loadings$X) and Q.loadings (only for SIMPLS, loadings$Y). Together with two sets of column vectors,
#' named "X.weights" and "Y.weights" that are needed to linearly combine the columns of X and Y into their respective latent scores.
#'
#' The optimal number of scores or components analysis is performed by tuning the accumulated variance threshold, see 'threshold.comp' in Arguments.
#' Moreover, Variable importance in projection (VIP) scores are calculated and an option of variable auto selection based on VIP is availabe, see 'auto.select.var' in arguments.
#'
#'
#' Note: Two methods are made available for this extraction step. The iterative method based on NIPALS PLS2 is refered to as the 'classic' method.
#' The "SIMPLS" method uses SVD decomposition function in R.
#'
#'@return easyPLSDA returns an object of class "easyPLSDA"
#'An object of class 'easyPLSDA' returns a list containing the following elements:
#' \item{Y}{Factor.The factor of the categorical response variable (Y).}
#' \item{levels}{Character vector of classes names.}
#' \item{Y.dummy}{Numeric matrix.The dummy matrix of the response (Y).}
#' \item{scaled}{Logical, TRUE if the explanatory variables matrix X was scaled (centered and normalized).}
#' \item{selected.var}{Character vector. The names of the selected explanatory variables in the PLS regression computation.}
#' \item{X}{Numeric matrix. The matrix of explanatory variables.}
#' \item{scores}{List of two numeric matrices. Scores matrices of X and Y variables projected into the new space.}
#' \item{weights}{Matrices of weights for both X and Y variables.}
#' \item{loadings}{One (or two for SIMPLS method) matrix of loadings.}
#' \item{B.mat}{Only for "classic" mode. Matrix of "b" coefficients, result of the scalar product between X latent vectors and Y latent vectors.}
#' \item{mode}{The mode used to compute the matrices.}
#' \item{explained.var}{Explained variance by latent vectors of X and Y.}
#' \item{VIP}{Numeric matrix. Matrix of Variable Importance in Projection.}
#' \item{comp.selected}{Numeric. The number of components of the fitted model.}
#'
#'@export
#'@examples
#'
#' #Iris example
#'
#' train <- sample.int(nrow(iris),100)
#' PLSobj <- easyPLSDA(Species~.,iris[train,],ncomp=3)
#' PLSpred <- predict(PLSobj,iris[-train,1:4])
#'
#' table(iris[-train,5],PLSpred$majority.vote)
#'
#'
#'
#'

easyPLSDA <- function(formula,data=NULL,ncomp=2,method="classic",auto.select.var=FALSE,threshold=0.8,threshold.comp=0.95,maxi.models=10,tol=10^-9,scale=TRUE){


  instance <- list()

  #Extraction of the formula variables in two part (explanatory varaibles X and response Y)
  extractedDF <- model.frame(formula,data)
  X <- extractedDF[,-1]
  Y <- extractedDF[,1]
  instance$Y <- Y
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

    #maxi.models (10) model by default to avoid Out of Memory Exception
    max.comp <- min(Matrix::rankMatrix(as.matrix(X))[1],nrow(X),maxi.models)

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

      cat(paste(k,"models tested.\n"))
      press.test[k] <- sum(press)


    }

    #We use the Rk statistic to select the optimal number of components
    Rk <- vector(length = max.comp-1)
    for (p in 2:length(press.test)) {
      Rk[p-1] <- press.test[p]/press.test[p-1]

    }
    ncomp.selected <- which(Rk>threshold.comp)[1]+1
    cat(paste(ncomp.selected,"components selected.\n"))
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
