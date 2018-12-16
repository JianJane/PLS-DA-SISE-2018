#' Summary function for PLSDA objects
#'
#' @param objectPLSDA
#'
#' @return Nothing
#' @export
#'
#'
summary.PLSDA <- function(objectPLSDA){
  cat("------PLSDA model------\n",sep="")
  cat("Classes : ")
  cat(objectPLSDA$levels,sep=", ")
  cat(".\n",sep="")
  cat("Number of explanatory variables : ",ncol(objectPLSDA$X),".\n",sep="")
  cat("Selected variables : ",length(objectPLSDA$selected.var),".\n\n",sep="")
  cat("Number of components : ",objectPLSDA$comp.selected,".\n",sep="")
  cat("Explained variance :\n\n",sep="")
  print(objectPLSDA$explained.var)
  cat("VIP :\n")
  print(objectPLSDA$VIP)
}

#' Summary function for PLSDAprediction objects
#'
#' @param object
#'
#' @return Nothing
#' @export
#'
#'
summary.PLSDAprediction <- function(object){
  cat("------PLSDA prediction object------\n",sep="")
  cat("Number of explanatory variables : ",ncol(object$new.data),".\n",sep="")
  cat("Number of new observations : ",nrow(object$new.data),".\n",sep="")


}
