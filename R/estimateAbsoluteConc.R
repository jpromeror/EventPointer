# Estimation of the concentration of isoforms
#
#' Internal function called by EvenPointer
#'
#' Estimate the concentration of the isoforms mapped to each of the three paths given the corresponding signal.
#'
#' @keywords internal
#'
#' @param Signal1 Signal for Path 1
#' @param Signal2 Signal for Path 2
#' @param SignalR Signal for Path R
#' @param lambda Penalty factor
#'
#'
#' @return list with concentrations of both isoforms and relative error of the estimation
#'

##########################################################
# Given the signals of the three paths,
# estimate the concentrations of each of the isoforms
##########################################################

# lambda parameter included to regularize the affinities

estimateAbsoluteConc <- function(Signal1, Signal2, SignalR, lambda ) {
  # require(nnls)
  Signal1 <- as.numeric(Signal1)
  Signal2 <- as.numeric(Signal2)
  SignalR <- as.numeric(SignalR)
  
  A <- cbind(Signal1, Signal2) 
  b <- SignalR
  Salida <- nnls(A,b) # non negative least squares
  if (lambda == 0) {
    Salida <- Salida$x
    u <- Salida[1]
    v <- Salida[2]
    w <- 0
    offset <- w / (1-u-v) # some times the offset is way too large (1-u-v = 0)
    T1est <- Signal1 * u
    T2est <- Signal2 * v
    Relerror <- as.numeric(crossprod((A[,1:2])%*%c(u,v)-b)/crossprod(b))
    return(list(T1est = T1est, T2est=T2est, offset = offset, Relerror = Relerror))
  }
  # Add a new equation to make the values of u and v close to each other
  if (is.null(lambda)) lambda <- 1
  penalty <- sum(Salida$residuals)/abs(diff(Salida$x))*lambda*abs(log(Salida$x[1]+1e-3)-log(Salida$x[2]+1e-3))/sqrt(3)
  A <- rbind(A,c(penalty,-penalty),c(penalty,0),c(0,penalty))
  b <- c(SignalR,0,penalty,penalty)
  Salida <- nnls(A,b)$x # non negative least squares
  u <- Salida[1]
  v <- Salida[2]
  w <- 0 # No offset used in this model
  offset <- w / (1-u-v) # some times the offset is way too large (1-u-v = 0)
  T1est <- Signal1 * u
  T2est <- Signal2 * v
  Relerror <- as.numeric(crossprod(cbind(Signal1, Signal2)%*%c(u,v)-SignalR)/crossprod(SignalR))
  # if(Relerror==0){browser()}
  return(list(T1est = T1est, T2est=T2est, offset = offset, Relerror = Relerror))
}

