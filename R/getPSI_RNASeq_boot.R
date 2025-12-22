#' getPSI_RNASeq_boot
#' 
#' @description Calculates PSI (Percent Spliced In) values with bootstrap resampling for RNA-Seq data.
#' This function estimates PSI values using a median polish approach and generates bootstrap samples
#' to assess the variability of PSI estimates.
#' 
#' @param Result An EventPointer result object from RNA-Seq analysis (e.g., from EventsDetection_BAM).
#'   Must contain count matrices that can be extracted with getCountMatrix.
#' @param lambda Regularization parameter for PSI estimation. Default is 0.1 if NULL.
#'   Controls the smoothness of PSI estimates.
#' @param cores The number of cores to use for parallel computation. Default is 1.
#' @param nboot Number of bootstrap iterations to perform. Default is 20.
#'   Higher values provide more stable estimates but increase computation time.
#' 
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts count matrices from the Result object
#'   \item Applies median polish to log-transformed count ratios to estimate normalization factors
#'   \item Performs parallel bootstrap resampling to generate multiple PSI estimates
#'   \item Combines original PSI with bootstrap samples into a single array
#' }
#' 
#' The bootstrap procedure uses parallel processing to distribute computations across multiple cores.
#' 
#' @return A 3-dimensional array with dimensions (events, bootstrap samples + 1, samples).
#'   The first layer (bootstrap sample 1) contains the original PSI estimates,
#'   followed by nboot bootstrap replicates.
#' 
#' @examples
#' \dontrun{
#' # Assuming 'eventsResult' is an EventPointer RNA-Seq result object
#' PSI_boot <- getPSI_RNASeq_boot(Result = eventsResult,
#'                                lambda = 0.1,
#'                                cores = 4,
#'                                nboot = 20)
#' 
#' # Access original PSI estimates
#' original_PSI <- PSI_boot[, 1, ]
#' 
#' # Access bootstrap samples
#' bootstrap_PSI <- PSI_boot[, 2:(nboot+1), ]
#' }
#' 
#' @export
#' @importFrom parallel makePSOCKcluster stopCluster clusterApplyLB
#' @importFrom aroma.light medianPolish
#' @importFrom abind abind
#' 

getPSI_RNASeq_boot<- function(Result, lambda = NULL, cores=1, nboot=20){
  if (is.null(lambda)) {
    lambda <- 0.1
  }
  # print(lambda)
  FCountMatrix <- getCountMatrix(Result)
  CountMatrix <- getCountMatrix(Result, modeFill = "Counts")
  
  # library(aroma.light) # Now imported via @importFrom in EventsDetection_BAM.R
  leq <- CountMatrix/FCountMatrix 
  Prueba <- medianPolish(log(leq), na.rm = TRUE)
  l1eq <- exp(Prueba$overall + Prueba$row[c(TRUE,FALSE,FALSE)])
  l2eq <- exp(Prueba$overall + Prueba$row[c(FALSE,TRUE,FALSE)])
  lReq <- exp(Prueba$overall + Prueba$row[c(FALSE,FALSE,TRUE)])
  
  l1eq[is.na(l1eq)] <- 1
  l2eq[is.na(l2eq)] <- 1
  lReq[is.na(lReq)] <- 1
  
  # nboot <- 20
  
  # PSI <- estimatePSI(CountMatrix, l1eq, l2eq, lReq, lambda = NULL, refine = F)
  cl <- makePSOCKcluster(cores) #not to overload your computer
  # registerDoParallel(cl)
  rowCounts <- seq_len(nrow(CountMatrix)/3)
  rowProcess <- split(rowCounts, ceiling(rowCounts/as.integer(length(rowCounts)/cores)))
  PSI_boot <- array(NA, c(ncol(CountMatrix), nrow(CountMatrix)/3, nboot))
  # PSI_boot <- calcBootstrapPSI(rowProcess = rowCounts,
  #                  PSI_boot=PSI_boot, 
  #                   CountMatrix = CountMatrix,
  #                   l1eq = l1eq,
  #                   l2eq = l2eq,
  #                   lReq = lReq,
  #                   lambda = lambda,
  #                   refine = F,
  #                   nboot = nboot)
  # resCalcPSI <- foreach(i = 1:length(rowProcess)) %dopar% {
  #   calcBootstrapPSI(rowProcess[[i]], PSI_boot, CountMatrix, l1eq, l2eq, lReq, lambda, refine, nboot)
  # }
  resCalcPSI <- clusterApplyLB(cl = cl, x = rowProcess, fun = calcBootstrapPSI,
                               PSI_boot=PSI_boot,
                               CountMatrix = CountMatrix,
                               l1eq = l1eq,
                               l2eq = l2eq,
                               lReq = lReq,
                               lambda = lambda,
                               nboot = nboot)
  stopCluster(cl)
  closeAllConnections()
  gc()
  for (i in c(1:length(resCalcPSI))) {
    minRow <- min(rowProcess[[i]])
    maxRow <- max(rowProcess[[i]])
    PSI_boot[, c(minRow:maxRow),] <-  resCalcPSI[[i]][, c(minRow:maxRow),]
  }
  PSI_boot <- aperm(PSI_boot,c(2,3,1))
  PSI <- estimatePSI(CountMatrix, l1eq, l2eq, lReq, lambda = lambda)
  PSI_boot <- abind(PSI, PSI_boot, along = 2)
  return(PSI_boot)
}