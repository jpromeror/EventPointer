#' Bootstrap PSI values from BAM files
#'
#' @description
#' Example dataset containing PSI (Percent Spliced In) values with bootstrap
#' replicates obtained from RNA-Seq BAM files analysis.
#'
#' @format A 3-dimensional array with dimensions (events, bootstrap samples, samples).
#' The first layer contains original PSI estimates, followed by bootstrap replicates.
#'
#' @source Generated from EventsDetection_BAM function using example BAM files.
#'
#' @examples
#' data(PSI_boots)
#' dim(PSI_boots)
#' 
"PSI_boots"