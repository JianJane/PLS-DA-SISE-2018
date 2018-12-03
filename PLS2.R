library(Matrix)

pls2 <- function(X,Y,ncomp=2,tol=10^-9,scale=TRUE){
  
  if(rankMatrix(X)[1]<ncomp){
    stop("Rank of the matrix is lower than number of composantes")
  }
  
  #
  if(scale==TRUE){
    X=scale(X)
    Y=scale(Y)
  }
  
  complatentscores <- function(X,Y){
    
    #Intialisation
    u <- Y[,1]
    t <- NA
    told <- NA
    
    while(is.na(t)|| is.na(told) || sum((t-told)^2)>(tol)){
      
      told <- t
      
      #Computation of X weights
      x_weights <- t(X)%*%u
      x_weights <- x_weights/as.numeric(sqrt(t(x_weights)%*%x_weights))
      
      #Computation of T latent vector
      t <- X%*%x_weights
      t <- t/as.numeric(sqrt(t(t)%*%t))
      
      #Computation of Y weights
      y_weights <- t(Y)%*%t
      y_weights <- y_weights/as.numeric(sqrt(t(y_weights)%*%y_weights))
      
      #Computation of u latent vector
      u <- Y%*%y_weights
      
      
    }
    result <- list(t_scores=t,u_scores=u,x_weights=x_weights,y_weights=y_weights)
    return(result)
  }
  
  #Matrix intilisation
  E0=X
  F0=Y
  T.scores = matrix(ncol = ncomp,nrow = nrow(E0))
  U.scores = matrix(ncol = ncomp,nrow = nrow(E0))
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
    
    b <- as.numeric(t(latentlist$t_scores)%*%latentlist$u_scores)
    B[i,i] <- b
    
    p <- t(E0)%*%latentlist$t_scores
    P.loadings[,i] <-p
    
    E0 <- E0-(latentlist$t_scores%*%t(p))
    F0 <- F0-b*(latentlist$t_scores%*%t(latentlist$y_weights))
    
  }
  
  instance <- list()
  
  instance$scores$X <- T.scores
  instance$scores$Y <- U.scores
  instance$weights$X <- X.weights
  instance$weights$Y <- Y.weights
  instance$loadings <- P.loadings
  instance$B.mat <- B
  
  class(instance) <- "PLS2"
  
  return(instance)
}

predict.PLS2 <- function(pls2Object,newdata){
  instance <- list()
  
  Bpls <- pls2Object$weights$X%*%solve(t(pls2Object$loadings)%*%pls2Object$weights$X)%*%t(pls2Object$weights$Y)
  pred <- newdata%*%Bpls
    
   
  instance$B.hat <- Bpls
  instance$pred <- pred
  
  
  class(instance) <- "predict"
  
  return(instance)
  
}