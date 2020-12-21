#' ResulTable
#' 
#' @description Extract a table of the top-ranked events from the output of EventPointer_Bootstraps.
#' 
#' @usage ResulTable(EP_Result,coef = 1,number = Inf)
#' 
#' @param EP_Result The output of the function EventPointer_Bootstraps
#' @param coef Number specifying which coefficient or contrast of the model is of interest.
#' @param number Maximum number of events to list
#' 
#' @return A dataframe with a row for the number of top events and the following columns:
#' 
#' deltaPSI: the difference of PSI between conditions
#' 
#' pvalue: raw p-value
#' 
#' lfdr: local false discovery rate
#' 
#' qvalue: adjusted p-value or q-value
#' 
#' 
#' @examples
#'
#'      data(PSIss)
#'      PSI <- PSIss$PSI
#'      
#'      Dmatrix <- cbind(1,rep(c(0,1),each=2))
#'      Cmatrix <- matrix(c(0,1),nrow=2)
#'      
#'      Fit <- EventPointer_Bootstraps(PSI = PSI,
#'                                     Design = Dmatrix,
#'                                     Contrast = Cmatrix,
#'                                     cores = 1,
#'                                     ram = 1,
#'                                     nBootstraps = 10,
#'                                     UsePseudoAligBootstrap = TRUE)
#'      
#'      ResulTable(EP_Result = Fit,coef = 1,number = 5)
#'
#' @export

ResulTable <- function(EP_Result,coef=1,number=Inf){
  
  numofcoeffs <- dim(EP_Result$Pvalues)[2]
  
  if(is.null(EP_Result)){
    stop("EP_Result field is empty")
  }
  
  if(coef > numofcoeffs){
    stop("coef out of bound")
  }
  
  if(numofcoeffs == 1){
    deltaPSI <- EP_Result$deltaPSI
    lfdr <- EP_Result$LocalFDR$lfdr2
    qvalues <- EP_Result$LocalFDR$qvalues
  }else{
    deltaPSI <- EP_Result$deltaPSI[,coef]
    lfdr <- EP_Result$LocalFDR[[coef]]$lfdr2
    qvalues <- EP_Result$LocalFDR[[coef]]$qvalues
  }
  # deltaPSI
  pvalues <- EP_Result$Pvalues[,coef]
  
  
  
  # identical(rownames(deltaPSI),names(pvalues))
  
  table <- data.frame(deltaPSI = deltaPSI,
                      pvalue = pvalues,
                      lfdr = lfdr,
                      qvalues = qvalues)
  
  oo <- order(table[,2],decreasing = FALSE)
  table <- table[oo,]
  if(number == Inf){
    return(table)
  }else{
    return(table[seq_len(number),])
  }
  
}