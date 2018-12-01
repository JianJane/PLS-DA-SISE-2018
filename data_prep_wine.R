# preration data Wine tasting

X_wine<-matrix(c(7, 7, 13, 7, 4, 3, 14, 7, 10, 5, 12, 5, 16, 7, 11, 3, 13, 3, 10, 3 ),5,4, byrow=TRUE)
Y_wine<-matrix(c(14, 7, 8, 10, 7, 6, 8, 5, 5, 2, 4, 7, 6, 2, 4 ),5,3, byrow=TRUE)

mX<-as.matrix(X_wine)
mY<-as.matrix(Y_wine)

mX<-scale(mX, center=TRUE, scale=TRUE)
mY<-scale(mY, center=TRUE, scale=TRUE)
