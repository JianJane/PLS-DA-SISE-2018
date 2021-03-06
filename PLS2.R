library(Matrix)
library(dummies)
pls2 <- function(X,Y,ncomp=2,method="classic",tol=10^-9,scale=TRUE){
  
  instance <- list()
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
  
  if(method=="SIMPLS"){
    # ---------------------------------PLS-2 Algorithm using SVD-------------------------------
    # -----------------------------------------------------------------------------------------
    
    E0 <- X
    F0 <- Y
    
    # corss product matrix S 
    M0 <- t(E0)%*%F0
    #UDV <- svd(M0)
    
    # 
    T.scores = matrix(ncol = ncomp,nrow = nrow(E0))
    U.scores = matrix(ncol = ncomp,nrow = nrow(F0))
    
    X.weights = matrix(ncol=ncomp,nrow=ncol(E0))
    Y.weights = matrix(ncol=ncomp,nrow=ncol(F0))
    
    P.loadings = matrix(nrow = ncol(E0), ncol = ncomp)
    Q.loadings = matrix(nrow = ncol(F0), ncol=ncomp)
    
    
    
    
    for(i in 1:ncomp){
      
      UDV <- svd(M0)
      
      X.weights[,i] <-UDV$u[,1]
      Y.weights[,i] <- UDV$v[,1]
      
      tau<- E0 %*% UDV$u[,1]
      tau_norm <- tau/sqrt(sum(tau^2))
      
      T.scores[,i] <-tau_norm
      
      u <- F0 %*% UDV$v[,1]
      u_norm <- u/sqrt(sum(u^2))
      
      U.scores[,i] <- u_norm
      
      # X and Y loadings are obtained by regressing against the same vector t 
      p <- t(E0)%*%tau_norm
      q <- t(F0)%*%tau_norm
      
      P.loadings[,i] <- p
      Q.loadings[,i] <- q
      
      E_ <- E0 - tau_norm%*%t(p)
      F_ <- F0 - tau_norm%*%t(q)
      
      
      M_ <- t(E_)%*%F_
      M0<-M_
      E0<-E_
      F0<-F_
      
      
      
    }
    instance$scores$X <- T.scores
    instance$scores$Y <- U.scores
    instance$weights$X <- X.weights
    instance$weights$Y <- Y.weights
    instance$loadings$X <- P.loadings
    instance$loadings$Y <- Q.loadings
    instance$mode <- "SVD"
    
    
  }else{
    if(method=="classic"){
      complatentscores <- function(X,Y){
        
        #Intialisation
        u <- Y[,1]
        x_weights <- NA
        x_weightsold <- NA
        
        while(is.na(x_weights)|| is.na(x_weightsold) || sum((x_weights-x_weightsold)^2)>(tol)){
          
          x_weightsold <- x_weights
          
          #Computation of X weights
          x_weights <- t(X)%*%u/ sum(u^2)
          x_weights <- x_weights/sqrt(sum(x_weights^2))
          
          #Computation of T latent vector
          t <- X%*%x_weights
          
          #Computation of Y weights
          y_weights <- t(Y)%*%t/sum(t^2)
          
          
          #Computation of u latent vector
          u <- Y%*%y_weights/sum(y_weights^2)
          
          
        }
        result <- list(t_scores=t,u_scores=u,x_weights=x_weights,y_weights=y_weights)
        return(result)
      }
      
      #Matrix intilisation
      E0=X
      F0=Y
      T.scores = matrix(ncol = ncomp,nrow = nrow(E0))
      U.scores = matrix(ncol = ncomp,nrow = nrow(F0))
      P.loadings = matrix(ncol=ncomp,nrow = ncol(E0))
      X.weights = matrix(ncol=ncomp,nrow=ncol(E0))
      Y.weights = matrix(ncol=ncomp,nrow=ncol(F0))
      B = matrix(0,ncol=ncomp,nrow=ncomp)
      
      for (i in 1:ncomp) {
        latentlist <- complatentscores(E0,F0)
        
        T.scores[,i] <- latentlist$t_scores
        U.scores[,i] <- latentlist$u_scores
        X.weights[,i] <-latentlist$x_weights
        Y.weights[,i] <- latentlist$y_weights
        
        b <- as.numeric(t(latentlist$u_scores)%*%latentlist$t_scores)
        B[i,i] <- b
        
        p <- t(E0)%*%latentlist$t_scores/sum(latentlist$t_scores^2)
        P.loadings[,i] <-p
        
        E0 <- E0-(latentlist$t_scores%*%t(p))
        F0 <- F0-(latentlist$t_scores%*%t(latentlist$y_weights))
        
      }
      
      
      
      instance$scores$X <- T.scores
      instance$scores$Y <- U.scores
      instance$weights$X <- X.weights
      instance$weights$Y <- Y.weights
      instance$loadings <- P.loadings
      instance$B.mat <- B
      instance$mode <- "classic"
    }else{
      
      stop("Error, choose a proper PLS method.")
      
      
    }
    
  }
  
  
  
  class(instance) <- "PLS2"
  
  return(instance)
}

predict.PLS2 <- function(pls2Object,newdata){
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