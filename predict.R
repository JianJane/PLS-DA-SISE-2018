# function that predict the Y with new X observations

#corpcor package is used to compute pseudoinverse matrix
library(corpcor)

predict <- function(Xnew, y_weights,p_loads,B){
  p<-p_loads[,-1]
  print(p)
  B<-B[2:nrow(B),2:ncol(B)]
  
  
  Ypred <- as.matrix(Xnew)%*%Bpls
  return(Ypred)
  
}

