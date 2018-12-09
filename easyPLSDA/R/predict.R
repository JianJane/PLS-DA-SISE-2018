predict.PLSDA <- function(pls2Object,newdata){
  newdata <- as.matrix(scale(newdata))
  instance <- list()
  if(pls2Object$mode=="classic" || pls2Object$mode=="test" ){

    Bpls <- pls2Object$weights$X%*%solve(t(pls2Object$loadings)%*%pls2Object$weights$X)%*%t(pls2Object$weights$Y)
    pred <- newdata%*%Bpls


    instance$B.hat <- Bpls
    instance$pred <- pred


    class(instance) <- "predict"

    return(instance)
  }else{

    if(pls2Object$mode=="SVD"){
      inv_pw<-inv(t(pls2Object$loadings$X)%*%pls2Object$weights$X)

      R.loadings <- pls2Object$weights$X%*%inv_pw

      Bpls <- R.loadings %*% t(pls2Object$weights$Y)

      pred <- newdata%*%Bpls

      instance$B.hat <- Bpls
      instance$pred <- pred


      class(instance) <- "predict"
      return(instance)
    }else{

      stop("Error with prediction object, check prediction mode.")


    }
  }


}
