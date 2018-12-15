#' easy Partial Least Square Discriminant Analysis class constructor function
#'
#' easyPLSDA class constructor is used to generate a class object to perform PLS-2 Discriminant analysis on subject data sets.
#' The returned object contains method to predict response matrix on unseen data : predict and  a summary() method.
#'
#' @usage easyPLSDA(formula,data=NULL,ncomp=2,method="classic",auto.select.var=FALSE,threshold=0.8,threshold.comp=0.95,maxi.comp=10,tol=10^-9,scale=TRUE)
#'
#'
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Similar to lm function.
#' @param data the dataset corresponding to given formula, must be coerced to a data.frame.
#' @param ncomp denoting the number of principal components to be used in the regression.
#' If NULL, the function choose the optimal number of components. Default: 2.
#' @param method The method that will be used to compute the different matrices used in regression.
#' Default "classic" is for PLS2 algorithm, "SIMPLS" is the other method based on matrix SVD.
#' @param tol The convergence threshold for latent scores computation, only for "classic" method. Default: 10^-9.
#' @param scale If TRUE (the default) after being column centered, X and Y can be further normalised with the standard deviation of each column.
#'
#'
#'@details The 'pls2' method function of the class extracts the numerical and geometric features from the independet variable matrix X and its linearly dependent response matrix Y.
#' Features of interest include the latent scores "T.scores" and "U.scores" from X and Y respectively, their respective loadings vectors "P.loadings" and "Q.loadings", along with two sets of column vectors,
#' named "X.weights" and "Y.weights" that are needed to linearly combine the columns of X and Y into their respective latent scores. The itterative method is set as default. Both methods return the features in a list.
#'
#'@details The 'predict' method function of the class predicts the response matrix from an unseen set of data
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
