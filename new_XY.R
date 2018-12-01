# function for computing the factor loading with converged 't'

new_XY<-function(l,mX, mY){

# initialising matrices
#-------------------------------  
  X<- mX
  Y<- mY
  
  # accessing items from list object 'l'
  #-------------------------------
  tau <- as.matrix(l[[1]][,1])
  u <- as.matrix(l[[1]][,2])
  w<- as.matrix(l[[2]])
  c_<- as.matrix(l[[3]]) 
  
  # calculating 'b' and 'p'
  #------------------------
  b <- t(tau)%*%u
  p <- t(X)%*%tau
    
  # subtracting the part of spanned by 't' from feature and response matrix
  # -----------------------------------------------------------------------
  X_ <- X - tau%*%t(p) 
  Y_ <- Y - as.numeric(b)*tau%*%t(c_)  # c is a column vector, transposed to row vector
  
  #LV: I changed the return to add the b scalar
  return(list(X_,Y_,p,b))
}