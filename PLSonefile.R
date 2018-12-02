X_wine<-matrix(c(7, 7, 13, 7, 4, 3, 14, 7, 10, 5, 12, 5, 16, 7, 11, 3, 13, 3, 10, 3 ),5,4, byrow=TRUE)
Y_wine<-matrix(c(14, 7, 8, 10, 7, 6, 8, 5, 5, 2, 4, 7, 6, 2, 4 ),5,3, byrow=TRUE)

mX<-as.matrix(mtcars[,1:5])
mY<-as.matrix(mtcars[,6:7])

mX<-scale(mX, center=TRUE, scale=TRUE)
mY<-scale(mY, center=TRUE, scale=TRUE)


complatentscores <- function(X,Y){
  
  u <- mY[,1]
  
  t <- NA
  told <- NA
  compteur <- 1
  while(is.na(t)|| is.na(told) || sum((t-told)^2)>(10^(-7))){
    
    told <- t
    
    x_weights <- t(X)%*%u
    x_weights <- x_weights/as.numeric(sqrt(t(x_weights)%*%x_weights))
    t <- X%*%x_weights
    t <- t/as.numeric(sqrt(t(t)%*%t))

    y_weights <- t(Y)%*%t
    y_weights <- y_weights/as.numeric(sqrt(t(y_weights)%*%y_weights))
    
    u <- Y%*%y_weights
    
   
    print(compteur)
    compteur <- compteur +1
  }
  result <- list(t_scores=t,u_scores=u,x_weights=x_weights,y_weights=y_weights)
  return(result)
  
}

#Number of comp to compute:
comp.max=3
#Matrix intilisation
T.scores = matrix(ncol = comp.max,nrow = nrow(mX))
U.scores = matrix(ncol = comp.max,nrow = nrow(mX))
P.loadings = matrix(ncol=comp.max,nrow = ncol(mX))
X.weights = matrix(ncol=comp.max,nrow=ncol(mX))
Y.weights = matrix(ncol=comp.max,nrow=ncol(mY))
B = matrix(0,ncol=comp.max,nrow=comp.max)

E0=mX
F0=mY
for (i in 1:comp.max) {
  latentlist <- complatentscores(E0,F0)
  
  T.scores[,i] <- latentlist$t_scores
  U.scores[,i] <- latentlist$u_scores
  X.weights[,i] <-latentlist$x_weights
  Y.weights[,i] <- latentlist$y_weights
  
  b <- as.numeric(t(latentlist$t_scores)%*%latentlist$u_scores)
  B[i,i] <- b
  
  p <- t(E0)%*%latentlist$t_scores
  P.loadings[,i] <-p
  
  E0 <- E0-latentlist$t_scores%*%t(p)
  F0 <- F0-b*(latentlist$t_scores%*%t(latentlist$y_weights))
  
  
  
}


print("T matrix :")
print(T.scores)
print("U matrix :")
print(U.scores)
print("P loadings matrix :")
print(P.loadings)
print("X weights W matrix :")
print(X.weights)
print("Y weights C matrix :")
print(Y.weights)
print(" B vector")
print(diag(B))

#Need library MASS for pseudo inverse
# predict1 <- function(Xnew,P.loadings,B,Y.weights){
#   Bpls <- ginv(t(P.loadings))%*%B%*%t(Y.weights)
#   print(Bpls)
#   return(Xnew%*%Bpls)
# }

predict2 <- function(Xnew,P.loadings,X.weights,Y.weights){
  Bpls <- X.weights%*%solve(t(P.loadings)%*%X.weights)%*%t(Y.weights)
  return(Xnew%*%Bpls)
}