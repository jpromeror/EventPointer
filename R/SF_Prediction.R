#' Splicing Factor Prediction
#' 
#' Methodology to predict context-specific splicing factors
#' 
#' @param P_value_PSI A data.frame with the p.values of the experiment.
#' @param ExS The ExS matrix biuldt in CreateExSmatrix function.
#' @param nSel Top ranked events to be considered as spliced events.
#' @param significance Threshold of P.value to consider which events are 
#' deferentially spliced. A vector of length equal to the number of contrasts.
#' If null it will consider the nSel top ranked events.
#' @param method methodology to apply: "Fisher" for Fisher's exact test 
#' (default), "PoiBin" for Poisson Binomial test, "Wilcoxon" for a wilcoxon 
#' test or  "Gsea" for a test of kolmogorov smirnov
#' 
#' @return The function returs a list. This list has for each contrast a data.frame containing
#' the results of the prediction.
#' 
#' 
#' 
#' @import glmnet
#' @import poibin
#' @import Matrix
#' @importFrom stats binomial phyper qhyper pnorm dnorm gaussian coef
#' @importFrom IRanges IRanges
#' @importFrom fgsea fgsea
#' @importFrom matrixStats rowRanks



SF_Prediction <- function(P_value_PSI,ExS,nSel=1000,significance=NULL,method="Fisher",valueRanking="Pvalue"){
  
  
  resPred <- vector(mode="list",length = ncol(P_value_PSI))
  ExS <- ExS[rownames(P_value_PSI), ]
  N <- nrow(ExS) #the same for each Fisher's test
  
  switch(method,
         Fisher={
           resPred <- hyperGeometricApproach(ExS, nSel, P_value_PSI, significance, resPred,N)
         },
         PoiBin={
           resPred <- poissonBinomialApproach(ExS, nSel, P_value_PSI, significance, resPred,N)
         },
         Wilcoxon={
           # nmTopEv <- significanceFunction (P_value_PSI, nSel=NULL, significance)
           # ExS <- ExS[nmTopEv, ]
           # resPred <- Wilcoxon.z.matrix(ExprT = t(abs(P_value_PSI[,4])),GeneGO = ExS)
           if (is.null(significance)) {
             significance = c(0.05,0.05,0.05,0.05)
           }
           resPred <- WilcoxonApproach(P_value_PSI, ExS, significance=significance, resPred, nSel = nSel, N = N)
          
           
         },
         Gsea ={
           resPred <- GseaApproach(P_value_PSI,ExS, significance=c(0.05,0.05,0.05,0.05), resPred)
           
         } )
  
  return(resPred)
}