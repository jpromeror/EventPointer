# Flux traversing the splicing graph
#
#' Internal function called by EvenPointer
#'
#' Function to determine the flux that is traversing each of the edges in the splicing graph
#'
#' @keywords internal
#'
#' @param Incidence Incidence matrix of the splicing graph
#' @param ncol number of columns for the final matrix (>1 will return a matrix of fluxes rather than a vector)
#'
#' @return matrix/vector of fluxes traversing the splicing graph
#'

getRandomFlow <- function(Incidence, ncol = 1)
{
  # With the incidence matrix, it is possible to get its null-space and generate an
  # arbitrary flow on it. Using the flow it is possible to get the triplets of events.

  # The seed is set to ensure the order of events remains the same
    set.seed("0xABBA")

  # Solve the Null Space for the Incidence Matrix
    solh <- Null(t(Incidence))

  # COndition to ensure that everything that exits the Start node (-1),
  # exits at the End Node (1)
    solp <- ginv(Incidence) %*% c(-1,rep(0,nrow(Incidence)-2),1)

  # Matrix of fluxes, with as many columns as specified by the user
    v <- matrix(runif(ncol(solh)*ncol),ncol=ncol)
    randSol <- as.vector(solp) + solh %*% v

  return(randSol)

}
