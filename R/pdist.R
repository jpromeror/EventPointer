#' Distance between matrices
#'
#' Calculates the distance between two matrices given by the user
#'
#'
#' @keywords internal
#'
#' @param X matrix
#' @param Y matrix
#'
#' @return Returns a matrix of distances
#'

pdist2<-function(X,Y)
{
  X1<-rowSums(X*X)
  Y1<-rowSums(Y*Y)
  Z<-outer(X1,Y1,"+")-2*X%*%t(Y)
}
