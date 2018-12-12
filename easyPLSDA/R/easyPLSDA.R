easyPLSDA <- function(formula,data=NULL,ncomp=2,method="classic",auto.select.var=TRUE,threshold=0.8,threshold.comp=0.95,tol=10^-9,scale=TRUE){

  instance <- list()
  extractedDF <- model.frame(formula,data)
  X <- extractedDF[,-1]

  Y <- extractedDF[,1]


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

  if(is.null(ncomp)){
    #Test du Leave one Out
    max.comp <- min(Matrix::rankMatrix(X)[1],nrow(X))

    press.test <- vector(length = max.comp)

    press = matrix(nrow=nrow(X), ncol=ncol(Y))
    #For the first composante :
    for (i in 1:nrow(X)) {
      computedPLSDA <- plsDA(Xp[-i,],Y[-i,],ncomp=1,method=method,auto.select.var=F,threshold=threshold,tol=tol)
      test <- Xp[i,]
      dim(test) <- c(1,length(Xp[i,]))
      colnames(test) <- colnames(X)

      prediction <- predict(computedPLSDA,test,scale=FALSE)
      press[i,] <- (Y[i,]-prediction$pred)^2
    }

    press.test[1] <- sum(press)

    plslist <- list()
    #For the next composantes
    for (k in 2:max.comp) {
      press = matrix(nrow=nrow(X), ncol=ncol(Y))
      for (i in 1:nrow(X)) {
        computedPLSDA <- plsDA(Xp[-i,],Y[-i,],ncomp=k,method=method,auto.select.var=F,threshold=threshold,tol=tol)
        plslist <- append(plslist,computedPLSDA)
        test <- Xp[i,]
        dim(test) <- c(1,length(Xp[i,]))
        colnames(test) <- colnames(X)
        prediction <- predict(computedPLSDA,test,scale=FALSE)
        press[i,] <- (Y[i,]-prediction$pred)^2
      }


      press.test[k] <- sum(press)


    }
    Rk <- vector(length = max.comp-1)
    for (p in 2:length(press.test)) {
      Rk[p-1] <- press.test[p]/press.test[p-1]

    }
    ncomp.selected <- which(valeurs>threshold.comp)[1]+1
    if(auto.select.var==T){
      computedPLSDA <- plsDA(Xp,Y,ncomp=ncomp.selected,method=method,auto.select.var=T,threshold=threshold,tol=tol)
      instance <- append(instance,computedPLSDA)
      instance$comp.selected <- ncomp.selected
      instance$Rk <- Rk
    }else{
      instance <- append(instance,plslist[which(valeurs>threshold.comp)[1]])
      instance$comp.selected <- ncomp.selected
      instance$Rk <- Rk
    }

  }else{
    instance <- plsDA(Xp,Y,ncomp=ncomp,method=method,auto.select.var=auto.select.var,threshold=threshold,tol=tol)
  }

  class(instance) <- "PLSDA"

  return(instance)
}
