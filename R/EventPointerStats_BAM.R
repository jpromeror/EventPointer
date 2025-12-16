#' Statistical analysis of alternative splicing events from BAM files
#'
#' @description
#' Performs statistical testing on PSI values derived from BAM files to identify
#' differential splicing events between conditions. Uses bootstrap resampling
#' for robust statistical inference.
#'
#' @param PSI_boots PSI_boots.RData obtained from \code{\link{EventsDetection_BAM}}
#'   function. Contains PSI values and bootstrap information for each event.
#' @param Design A matrix defining the linear model. Each row corresponds to a
#'   sample, and each column corresponds to a coefficient (such as the baseline
#'   and treatment effects).
#' @param Contrast A numeric matrix with contrasts to be tested. Rows correspond
#'   to coefficients in the design matrix, and columns correspond to contrasts.
#' @param Threshold Numeric. Threshold value for computing p-values. Default is 0.
#' @param nbootstraps Integer. Number of bootstrap iterations to use for statistical
#'   testing. Higher numbers increase computational time but improve statistical
#'   power. Default is 1000.
#' @param cores Integer. Number of cores to use for parallel processing. Default is 1.
#' @param ram Numeric. Amount of RAM (in GB) to use for computations. Default is 0.1.
#' @param pathResult Path where results will be saved. A subdirectory
#'   "EventPointerStatsResult" will be created containing a table with the results
#'   of the differential $\Psi$ analysis for each contrast. The table presents
#'   the $\Delta \Psi$ associated with each event and its corresponding significance
#'   parameters. Note that a separate table will be generated for each contrast
#'   specified in the contrast matrix. Default is current directory ("./").
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
#' Design <- cbind(rep(1, 9), rep(c(1, 0, 0), 3), rep(c(0, 1, 0), 3))
#' Contrast <- cbind(c(0, 1, 0), c(0, 0, 1))
#'
#' PathSGResult <- tempdir()
#'
#' EventPointerStats_BAM(
#'   PSI_boots,
#'   Design,
#'   Contrast,
#'   Threshold = 0,
#'   nbootstraps = 1000,
#'   cores = 1,
#'   ram = 0.1,
#'   pathResult = PathSGResult
#' )
#' }
EventPointerStats_BAM <- function(PSI_boots,
                                  Design,
                                  Contrast,
                                  Threshold = 0,
                                  nbootstraps = 1000,
                                  cores = 1,
                                  ram = 0.1,
                                  pathResult = "./"){

  if(is.null(Design)){
    stop("Design field is empty")
  }
  if(is.null(Contrast)){
    stop("Contrast field is empty")
  }
  pathResult <- paste0(pathResult, "/EventPointerStatsResult/")
  
  if (!file.exists(pathResult)) {
    dir.create(pathResult)
  } 
  
  result <- checkContrastDesignMatrices(Contrast, Design)
  if (result == TRUE){
    UseBootstrap <- T
    resBootstrap <- EventPointer_Bootstraps(PSI=PSI_boots, Design=Design,
                                            Contrast=Contrast,nbootstraps=nbootstraps,
                                            UsePseudoAligBootstrap =T,
                                            Threshold =Threshold,
                                            cores=cores, ram=ram)
    pathResultBootstrap <- paste0(pathResult, "bootstrapResult/")
    dir.create(pathResultBootstrap,showWarnings = FALSE)
    for (coef in c(1:dim(resBootstrap$Pvalues)[2])){
      tableRes <- ResulTable(resBootstrap, coef = coef)
      write.csv(tableRes,file = paste0(pathResultBootstrap,"ResBootstrapContrast",coef,".csv"))
    }
    return(resBootstrap)
  }
  
  
}
