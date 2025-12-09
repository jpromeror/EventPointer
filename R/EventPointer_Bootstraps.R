#' EventPointer_Bootstraps
#' 
#' @description Statistical analysis of alternative splicing events with bootstrap technique.
#' 
#' @param PSI Array or matrix that contains the values of PSI calculated in the function GetPSIFromTranRef.
#'   If bootstrap option was selected in GetPSIFromTranRef, input must be an array. If not, input must be a matrix.
#' @param Design A matrix defining the experimental design. Rows represent samples and columns conditions.
#' @param Contrast A numeric matrix with contrasts to be tested. Rows correspond to coefficients in the design matrix, and columns correspond to contrasts.
#' @param nbootstraps How many layers, Bootstraps or samplings are going to be used. Caution, high numbers increase computational time.
#' @param UsePseudoAligBootstrap TRUE (default) if bootstrap data from pseudoalignment want to be used or FALSe if not.
#' @param Threshold it assigns a threshold to compute the pvalues. Default value is 0.
#' @param cores The number of cores desired to use.
#' @param ram How many ram memory is used,in Gb.
#' 
#' @examples
#'        data(PSIss)
#'        PSI <- PSIss$PSI
#'        
#'        Dmatrix <- cbind(1,rep(c(0,1),each=2))
#'        Cmatrix <- matrix(c(0,1),nrow=2)
#'        
#'        Fit <- EventPointer_Bootstraps(PSI = PSI,
#'                                       Design = Dmatrix,
#'                                       Contrast = Cmatrix,
#'                                       nbootstraps = 10,
#'                                       UsePseudoAligBootstrap = TRUE,
#'                                       Threshold = 0,
#'                                       cores = 1,
#'                                       ram = 1)
#' 
#' @return A list containing the summary of the Bootstrap analysis: DeltaPSI, Pvalues, FDR. This info can be
#' obtained in a simple table with the function ResulTable.
#' 
#' 
#' @export
#' @import Matrix
#' @import doParallel
#' @import foreach
#' @import iterators
#' @importFrom parallel makeCluster stopCluster
#' @importFrom lpSolve lp
#' @importFrom matrixStats iqr
#' @importFrom abind abind
#' @importFrom qvalue qvalue
#' @importFrom cobs cobs
#' @importFrom stats pnorm ppoints nlminb
#' @importFrom IRanges median IQR
#' @importFrom S4Vectors na.omit
#' 

EventPointer_Bootstraps <- function(PSI, Design, Contrast, nbootstraps=10000,
                                    UsePseudoAligBootstrap=TRUE, Threshold = 0,
                                    cores=1,ram=0.1){
  
  
  if(is.null(PSI)){
    stop("PSI field is empty")
  }
  if(is.null(Design)){
    stop("Design field is empty")
  }
  if(is.null(Contrast)){
    stop("Contrast field is empty")
  }
  
  
  # set.seed(1)
  result <- checkContrastDesignMatrices(Contrast, Design)
  
  if (result == TRUE){
    # Bootstrap Test: ----
    
    table <- mclapplyPSI_Bootstrap(PSI_boots = PSI,
                                   Design = Design,
                                   Contrast = Contrast,
                                   nbootstraps = nbootstraps,
                                   KallistoBootstrap = UsePseudoAligBootstrap,
                                   th = Threshold, 
                                   cores = cores,
                                   ram = ram)

    cat("\n The program has succesfully ended. \n", sep ="\n")
  }
  return(table)
}

   #Histograms ----
    # if (dim(Cmatrix)[1]==1){
    #   hist(table[[2]],1000, main = paste("Histogram of contrast 1"), xlab = "p-values")
    #   }else{
    #     for (i in 1:dim(Cmatrix)[2]) {
    #       hist(table[[2]][,i],1000, main = paste("Histogram of contrast " , i), xlab = "p-values")
    #     }
    #   }