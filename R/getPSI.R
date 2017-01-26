# Percent Spliced In estimation
#
#' Internal function called by EvenPointer
#'
#' Function to calculate the PSI for a particular alternative splicing event
#'
#' @keywords internal
#'
#' @param ExFit data.frame with expression of alternative splicing events
#'
#' @return matrix with PSI values along events and samples
#'

getPSI <- function(ExFit) {
  # Create matrix to fill with PSI values (1 per event and sample)
  PSI <- matrix(0, nrow = nrow(ExFit)/3, ncol = ncol(ExFit)-5)
  colnames(PSI)  <- colnames(ExFit[6:ncol(ExFit)])
  rownames(PSI)  <- ExFit[seq(1,nrow(ExFit),by = 3),1]

  NCols<-ncol(ExFit)
  # Perform the operations for every detectable alternative splicing event
  for (n in 1:(nrow(ExFit)/3)) {

    # Get expression signal from path 1
    Signal1 <- ExFit[1+3*(n-1),6:NCols]

    # Get expression signal from path 2
    Signal2 <- ExFit[2+3*(n-1),6:NCols]

    # Get expression signal from Reference
    SignalR <- ExFit[3+3*(n-1),6:NCols]

    # Function to estimate concentrations from the interrogated isoforms
    Output <- estimateAbsoluteConc(Signal1, Signal2, SignalR, lambda = 1)

    # Compute the actual PSI value (T1/T1+T2)
    psi <- Output$T1est / (Output$T1est + Output$T2est)
    PSI[n,] <- psi
  }
  return(PSI)
}
