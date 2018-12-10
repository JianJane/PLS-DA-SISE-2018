function (X, Y, ncomp = 2, scale = TRUE, mode = c("regression", 
                                                  "canonical", "invariant", "classic"), tol = 1e-06, max.iter = 100, 
          near.zero.var = FALSE, logratio = "none", multilevel = NULL, 
          all.outputs = TRUE) 
{
  if (is.null(multilevel)) {
    if (is.null(Y)) 
      stop("'Y' has to be something else than NULL.")
    if (is.null(dim(Y))) {
      Y = factor(Y)
    }
    else {
      stop("'Y' should be a factor or a class vector.")
    }
    if (nlevels(Y) == 1) 
      stop("'Y' should be a factor with more than one level")
    Y.mat = unmap(Y)
    colnames(Y.mat) = levels(Y)
  }
  
  else {
    multilevel = data.frame(multilevel)
    if ((nrow(X) != nrow(multilevel))) 
      stop("unequal number of rows in 'X' and 'multilevel'.")
    if (ncol(multilevel) != 1) 
      stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
    if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0, 1, 2)) 
      stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
    multilevel = data.frame(multilevel, Y)
    multilevel[, 1] = as.numeric(factor(multilevel[, 1]))
    Y.mat = NULL
  }
  
  result = internal_wrapper.mint(X = X, Y = Y.mat, ncomp = ncomp, 
                                 scale = scale, near.zero.var = near.zero.var, mode = mode, 
                                 max.iter = max.iter, tol = tol, logratio = logratio, 
                                 multilevel = multilevel, DA = TRUE, all.outputs = all.outputs)
  
  out = list(call = match.call(), X = result$A[-result$indY][[1]], 
             Y = if (is.null(multilevel)) {
               Y
             } else {
               result$Y.factor
             }, ind.mat = result$A[result$indY][[1]], ncomp = result$ncomp, 
             mode = result$mode, variates = result$variates, loadings = result$loadings, 
             loadings.star = result$loadings.star, names = result$names, 
             tol = result$tol, iter = result$iter, max.iter = result$max.iter, 
             nzv = result$nzv, scale = scale, logratio = logratio, 
             explained_variance = result$explained_variance, input.X = result$input.X, 
             mat.c = result$mat.c)
  
  class(out) = c("plsda", "pls", "DA")
  
  if (!is.null(multilevel)) {
    out$multilevel = multilevel
    class(out) = c("mlplsda", class(out))
  }
  
  return(invisible(out))
}