# function for finding u and Tau

itr_latentScores <- function(mX, mY)
{
  
# function returns a data frame whoes columns are 't', 'u', 'w'
# the latent scores of X and Y respectively
  
  U <-as.matrix(rnorm(length(mY[,1])),ncol=1)
  
  #print(U)
  t_ <-matrix(0, length(U),1) # initialising t_ for t convergence
  
  t_conv<-c() # storing norms of difference between two consecutive 'u' 
  n_itr <- 20 # number of iteration for convergence (converged after 4 iterations)
  
for ( i in 1:n_itr) {
    
    # estimate X weights, and normalise
      w <- t(mX)%*%U   
      norm_w <- sqrt(sum(w^2))  # the euclidian norm
      #print(w)
      #print(norm_w)
      w<- w/norm_w # vector w
      #print(w)
      #print(i)
    # estimate X scores, and normalise
      t <- as.matrix(mX)%*%as.numeric(w)
      norm_Tau <- sqrt(sum(t^2))
      t <- t/norm_Tau
    # estimate Y weights, and normalise
      c<-t(mY)%*%t
      c_norm <- sqrt(sum(c^2))
      c<-c/c_norm
      #print(c)
    # estimate Y scores
      U<-mY%*%c
      print(U)
      print(i)
    # storing norm of difference 
    t_conv<-c(t_conv,sqrt(sum((t-t_)^2)))
    #print(t_conv)
    t_ <- t
}       
  
  latents <- cbind(t,U)
  
  b <- t(t)%*%U
  p <- t(mX)%*%t
 
  colnames(latents) <- c("t","u")
  l<-list(data.frame(latents),as.numeric(w),as.numeric(c))
  
  return(l)
  
}