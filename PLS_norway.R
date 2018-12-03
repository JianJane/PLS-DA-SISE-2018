# preration data Wine tasting

source("data_prep.R")
#source("data_prep_wine.R")

# ---------------------------------PLS-2 Algorithm using SVD-------------------------------

#install.packages("psych")

library("psych")
library('matlib')

#library(pls)
#source("new_XY.R")
#source("itr_latentScores.R")
#source("predict.R")

# ---------------------------------PLS-2 Algorithm using SVD-------------------------------
# -----------------------------------------------------------------------------------------

E0 <- mX
F0 <- mY

# corss product matrix S 
M0 <- t(E0)%*%F0
#UDV <- svd(M0)

# 
T.scores = matrix(NA,nrow(E0),1)
U.scores = matrix(NA,nrow(F0),1)

X.weights = matrix(NA,ncol(E0),1)
Y.weights = matrix(NA,ncol(F0),1)

p.loadings = matrix(NA,nrow = ncol(E0), 1)
q.loadings = matrix(NA,nrow = ncol(F0), 1)


cndit <- TRUE
n<-0

while (cndit ==1){
  
  UDV <- svd(M0)
  
  X.weights<-cbind(X.weights,UDV$u[,1])
  Y.weights<-cbind(Y.weights,UDV$v[,1])
  
  tau<- E0 %*% UDV$u[,1]
  tau_norm <- tau/sqrt(sum(tau^2))
  
  T.scores<-cbind(T.scores,tau_norm)
  
  u <- F0 %*% UDV$v[,1]
  u_norm <- u/sqrt(sum(u^2))
  
  U.scores<-cbind(U.scores,u_norm)
  
  # X and Y loadings are obtained by regressing against the same vector t 
  p <- t(E0)%*%tau_norm
  q <- t(F0)%*%tau_norm
  
  p.loadings<- cbind(p.loadings,p)
  q.loadings<- cbind(q.loadings,q)
  
  E_ <- E0 - tau_norm%*%t(p)
  F_ <- F0 - tau_norm%*%t(q)
  
  cndit <- log10(sum(abs(E0-E_))) > 0
  
  M_ <- t(E_)%*%F_
  M0<-M_
  E0<-E_
  F0<-F_
  
  #print("n")
  #print(n<-n+1)
}

#inv_pw <- inv(t(p.loadings[,-1])%*%X.weights[,-1])
inv_pw<-inv(t(p.loadings[,c(-1,-6)])%*%X.weights[,c(-1,-6)])

R.loadings <- X.weights[,c(-1,-6)]%*%inv_pw

B <- R.loadings %*% t(q.loadings[,c(-1,-6)])

Y_new <- mX%*%B

# use max to find the selected class
res <- sapply(1:nrow(Y_new), function(x) {which.max(Y_new[x,])})
res <- as.matrix(res, nrow=length(res), 1)