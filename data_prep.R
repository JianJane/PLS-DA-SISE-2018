# loading and preparing data

#library("pschy")
data("iris")

iris <- data.frame(iris)
#iris <- iris[sample(1:nrow(iris)),]

# train-test split
t_perct<-1
n_train<- ceiling(t_perct*nrow(iris))
n_test<- nrow(iris)- n_train

# Y matri
tIris <- iris[1:n_train,]
#response dummy matrix 
Y <- tIris$Species 
mY<- dummy.code(Y) 
# feature matrix
tIris$Species<-NULL 
mX <- tIris

# test Ymatrix for validation
X_tst <- iris[(n_train+1):nrow(iris),]
Y_tst <- X_tst$Species
Y_tst <- dummy.code(Y_tst)
X_tst$Species<-NULL

# column centering and divided by column std, or Z-scoring matrix

#mX<-scale(mX, center=TRUE, scale=TRUE)
mX<-scale(mX, center=TRUE)

#mY<-scale(mY, center=TRUE, scale=TRUE)
 mY<-scale(mY, center=TRUE)
