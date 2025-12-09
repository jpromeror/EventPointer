#' Statistical analysis of alternative splicing events from transcript quantification
#'
#' @description
#' Performs statistical testing on PSI values from transcript-level quantification
#' to identify differential splicing events between conditions. Uses bootstrap
#' resampling for robust statistical inference.
#'
#' @param PSI PSI_ST.RData variable resulting from \code{\link{EventsDetection_ST}}.
#'   Contains PSI values and bootstrap information for each event.
#' @param Design A matrix defining the linear model. Each row corresponds to a
#'   sample, and each column corresponds to a coefficient (such as the baseline
#'   and treatment effects).
#' @param Contrast A numeric matrix with contrasts to be tested. Rows correspond
#'   to coefficients in the design matrix, and columns correspond to contrasts.
#' @param BootstrapStats Logical. If TRUE, bootstrap-based statistics are computed.
#'   Default is TRUE.
#' @param nbootstraps Integer. Number of bootstrap iterations to use for statistical
#'   testing. Higher numbers increase computational time but improve statistical
#'   power. Default is 10000.
#' @param UsePseudoAligBootstrap Logical. If TRUE (default), uses bootstrap data
#'   from pseudoalignment (Kallisto/Salmon) for statistical inference.
#' @param Threshold Numeric. Threshold value for computing p-values. Default is 0.
#' @param cores Integer. Number of cores to use for parallel processing. Default is 1.
#' @param ram Numeric. Amount of RAM (in GB) to use for computations. Default is 4.
#' @param pathResult Path where results will be saved. A subdirectory
#'   "EventPointerStatsSTResult" will be created. Default is current directory ("./").
#'
#' @return Creates result files in the specified path containing statistical
#'   analysis results, including p-values and effect sizes for each contrast.
#'   Files are saved in CSV format in the "bootstrapResult" subdirectory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(PSI_boots)
#'
#' Design <- cbind(1, rep(c(0, 1), each = 2))
#' Contrast <- matrix(c(0, 1), nrow = 2)
#'
#' pathResult <- tempdir()
#'
#' EventPointerStats_ST(
#'   PSI,
#'   Design,
#'   Contrast,
#'   BootstrapStats = TRUE,
#'   nbootstraps = 10000,
#'   UsePseudoAligBootstrap = TRUE,
#'   Threshold = 0,
#'   cores = 1,
#'   ram = 4,
#'   pathResult = pathResult
#' )
#' }
EventPointerStats_ST <- function(PSI,
                                 Design,
                                 Contrast,
                                 BootstrapStats = TRUE,
                                 nbootstraps = 10000,
                                 UsePseudoAligBootstrap = TRUE,
                                 Threshold = 0,
                                 cores = 1,
                                 ram = 4,
                                 pathResult = "./"){
  
  if(is.null(Design)){
    stop("Design field is empty")
  }
  if(is.null(Contrast)){
    stop("Contrast field is empty")
  }
  pathResult <- paste0(pathResult, "/EventPointerStatsSTResult/")
  checkMatrices <- checkContrastDesignMatrices(Contrast, Design)
  if(checkMatrices){
    if (!file.exists(pathResult)) {
      dir.create(pathResult)
    } 
    
    if(BootstrapStats){
      resBootstrap <- EventPointer_Bootstraps(PSI$PSI,
                                              Design,
                                              Contrast,
                                              cores = cores,
                                              ram = ram,
                                              nbootstraps = nbootstraps,
                                              UsePseudoAligBootstrap = UsePseudoAligBootstrap,
                                              Threshold = Threshold)
      pathResultBootstrap <- paste0(pathResult, "bootstrapResult/")
      dir.create(pathResultBootstrap)
      for (coef in c(1:dim(resBootstrap$Pvalues)[2])){
        tableRes <- ResulTable(resBootstrap, coef = coef)
        write.csv(tableRes,file = paste0(pathResultBootstrap,"ResBootstrapContrast",coef,".csv"))
      }
    }
    
  }
  
}
