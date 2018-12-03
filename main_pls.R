# ---------------------------------PLS-2 Algorithm-------------------------------

#install.packages("psych")
library("psych")
source("new_XY.R")
source("itr_latentScores.R")
source("predict.R")

# --------------------------prep of feature matrix X and response matrix Y-------
# ------------------------------------------------------------------------------
# preration data Wine tasting

source("data_prep.R")
#source("data_prep_wine.R")




#--------------------------calculating the 1st latent variables 'u', 't' and loading-----------------
#-------------------------call user defined function to find first latent scores

l<-itr_latentScores(mX, mY)
 
# initialising with the first score vector for storing subsequent scores-------- 
# ------------------------------------------------------------------------------

t_scores <- as.matrix(l[[1]][,1])
u_scores <- as.matrix(l[[1]][,2])
x_weights <- as.matrix(l[[2]])
y_weights <- as.matrix(l[[3]])

# LV : I change size of p_loads to generalize
# and add a B matrix to store the b scalars
p_loads <- matrix(,nrow = ncol(mX),1)
B <- matrix(NA,1)

#--------------------------subtracting identified components from X and Y matrices
# ------------------------------------------------------------------------------
 
  # initilising matrix 
  E_<- mX
  F_<- mY
  n_itr<-0
  tshold<--1
  
  # using function new_XY(l, E_, F_) subtracting energy spanned by principal components, and updating matrices E_, F_ 
   
  sn<- TRUE # initialising loop condition 
  n <-0
  while (sn == TRUE) {
      
      ml<-new_XY(l,E_,F_)  # call function new_XY
      
      E_ <- ml[[1]] # retrieving updated matrix E_
      F_ <- ml[[2]] # retrieving updated matrix F_
  
      l<-itr_latentScores(E_, F_) # call function itr_latentScores to find the next component
      
      # storing scores and weights 
      t_scores<-cbind(t_scores,as.matrix(l[[1]][,1]))
      u_scores<-cbind(u_scores,as.matrix(l[[1]][,2]))
      x_weights<-cbind(x_weights,as.matrix(l[[2]]))
      y_weights<-cbind(y_weights,as.matrix(l[[3]]))
      
      #print(dim(as.matrix(ml[[3]])))
    #  p_loads <- cbind(as.matrix(p_loads),as.matrix(ml[[3]]))
      p_loads <- cbind(p_loads,as.matrix(ml[[3]]))
    
     # n<-n+1
     # print(cat('load p', n))
      
      #LV : here i store the b scalars into B matrix (as diagonal terms)
      B <- diag(c(diag(B),ml[4]))
      
      #print(cat('load B', n))
      #n_itr<-n_itr+1
      #print(cat("number of iteration before convergence: ", n_itr))
      #print(log10(sum(abs(E_))))
      
      sn<-log10(sum(abs(E_))) > tshold # monitoring sn while condition as E approaching NULL
  }
    
  b <- t(t_scores[,ncol(t_scores)])%*%u_scores[,ncol(u_scores)]
  B <- diag(c(diag(B),b))
  p <- t(E_)%*%t_scores[,ncol(t_scores)]
  p_loads<- cbind(p_loads,p)
  
  # call "predict" function to fit the new feature matrix
 # X_tst <- mX
  
  
#  Y_pd<-predict(X_tst,y_weights,p_loads,B)
  

 #sapply(1:length(Y_pd), function(x) {which.max(Y_pd[x,])})
  

  