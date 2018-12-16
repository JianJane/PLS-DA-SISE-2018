#' Prediction method for PLSDA object
#'
#' @param pls2Object Object of class PLSDA.
#' @param newdata data.frame or matrix, New observations for prediction.
#' @param scale default : TRUE, if the data has to be scaled before prediction computation.
#'
#' @return S3 object of class PLSDAprediction
#' @export
#'
#'
predict.PLSDA <- function(pls2Object,newdata,scale=TRUE){
  if(is.null(dim(newdata))){
    stop("New data must be a matrix or a data.frame")
  }
  if(scale==TRUE){
    #We scale new data to not compute the intercept
    newdata <- as.matrix(scale(newdata))

  }else{
    newdata <- newdata
  }
  instance <- list()
  instance$new.data <- newdata

  #Two ways of calculate the prediction depending on the PLS method used
  if(pls2Object$mode=="classic" ){

    #Compute the matrix of regression coefficents
    Bpls <- pls2Object$weights$X%*%solve(t(pls2Object$loadings)%*%pls2Object$weights$X)%*%t(pls2Object$weights$Y)

    #Compute the response
    pred <- newdata[,pls2Object$selected.var]%*%Bpls


    instance$B.hat <- Bpls
    instance$pred <- pred

    #Determine the classes
    majorityvote <- apply(pred, 1, which.max)
    majorityvote <- sapply(majorityvote, function(x){return(pls2Object$levels[x])})
    instance$majority.vote <- majorityvote
    class(instance) <- "PLSDAprediction"

    return(instance)
  }else{

    if(pls2Object$mode=="SIMPLS"){

      inv_pw<-solve(t(pls2Object$loadings$X)%*%pls2Object$weights$X)
      R.loadings <- pls2Object$weights$X%*%inv_pw

      #Compute the matrix of regression coefficents
      Bpls <- R.loadings %*% t(pls2Object$weights$Y)

      #Compute the response
      pred <- newdata[,pls2Object$selected.var]%*%Bpls

      instance$B.hat <- Bpls
      instance$pred <- pred

      #Determine the classes
      majorityvote <- apply(pred, 1, which.max)
      majorityvote <- sapply(majorityvote, function(x){return(pls2Object$levels[x])})
      instance$majority.vote <- majorityvote
      class(instance) <- "PLSDAprediction"
      return(instance)
    }else{

      stop("Error with prediction object, check prediction mode.")


    }
  }


}
