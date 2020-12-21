#' Protein_Domain_Enrichment
#'
#' @description Analyze whether the presence of a protein domain increases or decreases in the condition under study.
#'  
#' @param PathsxTranscript the output of EventDetection_transcriptome.
#' @param TxD matrix that relates transcripts with Protein domain. Users can get it from BioMart
#' @param Diff_PSI matrix with the difference of psi of the condition under study. Can get it from the output of EventPointer_Bootstraps 
#' @param method a character string indicating which correlation coeffcient is to be calculated. "spearman" (default) or "pearson"
#'  can be selected.
#' 
#' 
#' @return A list containing the results of the protein domain enrichment anaylisis. This list contains 3 matrices in which the rows indicate
#'  the protein domains and the columns the number of contrasts. The 3 matrices are the following:
#'  
#'  -mycor: correlation value between the deltaPSI and the DifProtDomain matrix (see more details in vignette)
#'  
#'  -STATISTIC: the values of the test statistic
#'  
#'  -PVAL: the pvalues of the test statistic
#'  
#' 
#' @examples 
#' 
#'    data("EventXtrans")
#'    data("TxD")
#'    data("Fit")
#'    
#'    #same annotation in TxD and EventXtrans
#'    transcriptnames <- EventXtrans$transcritnames
#'    transcriptnames <- gsub("\\..*","",transcriptnames) 
#'    EventXtrans$transcritnames <- transcriptnames
#'    
#'    Result_PDEA <- Protein_Domain_Enrichment(PathsxTranscript = EventXtrans,
#'                                             TxD = TxD,
#'                                             Diff_PSI = Fit$deltaPSI)
#' 
#' @export
#' @import Matrix
#' @importFrom matrixStats rowRanks rowSums2
#' @importFrom stats model.matrix

Protein_Domain_Enrichment <- function(PathsxTranscript,TxD,Diff_PSI,method="spearman"){
  
  if(is.null(PathsxTranscript)){
    stop("PathsxTranscript field is empty")
  }
  if(is.null(TxD)){
    stop("TxD field is empty")
  }
  if(is.null(Diff_PSI)){
    stop("Diff_PSI field is empty")
  }
  
  
  
  transcriptnames <- PathsxTranscript$transcritnames
  
  jjx <- match(rownames(TxD),transcriptnames)
  
  ExTP1 <- PathsxTranscript$ExTP1[,jjx]
  ExTP2 <- PathsxTranscript$ExTP2[,jjx]
  
  #we put the same order:
  #events(path1) x transcripts
  #events(path2) x transcripts
  ## colnames same order with the rownames of TxD
  
  ## we obtain matrix of event(path_x) x Domain
  eventsxpfam_path1 <- ExTP1 %*% TxD
  eventsxpfam_path2 <- ExTP2 %*% TxD
  
  #every transcript of the event must have the protein domain
  numbtrans_path1 <- rowSums(ExTP1)
  numbtrans_path2 <- rowSums(ExTP2)
  
  eventsxpfam_path1 <- Diagonal(n = length(numbtrans_path1),x = (1/numbtrans_path1) ) %*% eventsxpfam_path1
  eventsxpfam_path1@x <- (eventsxpfam_path1@x==1)*1
  eventsxpfam_path1 <- drop0(eventsxpfam_path1)
  
  eventsxpfam_path2 <- Diagonal(n = length(numbtrans_path2),x = (1/numbtrans_path2) ) %*% eventsxpfam_path2
  eventsxpfam_path2@x <- (eventsxpfam_path2@x==1)*1
  eventsxpfam_path2 <- drop0(eventsxpfam_path2)
  
  
  FinalExPF <- eventsxpfam_path1 - eventsxpfam_path2
  #FinalExPF == 1: Domain that win the path 1
  #FinalExPF == -1: Domain that win the path 2
  #FinalExPF == 0: Both paths have the Domain
  #FinalExPF == . (equal to 0): Domain that neither path1 nor path2 have
  
  # table(FinalExPF@x)
  
  
  ## statistical analysis
  # rownames(Results_Table$deltaPSI)
  FinalExPF <- FinalExPF[match(rownames(Diff_PSI),rownames(FinalExPF)),]
  FinalExPF <- drop0(FinalExPF)
  
  ## calculate fasta correlation spearman
  ## we calculate correlation only for those protein domains that have -1,0,1 values
  dim(FinalExPF)
  FinalExPF_2 <- FinalExPF
  miend <- FinalExPF_2@p[-1]
  mistart <- c(1,miend+1)
  mistart <- mistart[-length(mistart)]
  jjx <- which((miend - mistart) <0) 
  ### these are those with all 0 in the column (we have to remove them)
  
  # iix <- which((miend - mistart + 1) == nrow(FinalExPF_2) ) 
  ### don't have any 0 (we want to keep them)
  
  if(length(jjx) > 0){
    FinalExPF_2 <- FinalExPF_2[,-jjx]
  }
  miend <- FinalExPF_2@p[-1]
  mistart <- c(1,miend+1)
  mistart <- mistart[-length(mistart)]
  # jjx <- which((miend - mistart) <0) #jjx must be integer(0)
  # midata <- data.frame(value=FinalExPF_2@x,lacolumn=NA)
  # midata$lacolumn <- rep(seq_len(ncol(FinalExPF_2)),diff(FinalExPF_2@p))
  lacolumn <- rep(seq_len(ncol(FinalExPF_2)),diff(FinalExPF_2@p))
  lacolumn <- as.factor(lacolumn)
  Dmatrix <- model.matrix(~0+lacolumn)
  dim(Dmatrix)
  length(FinalExPF_2@x)
  ppp <- crossprod(Dmatrix,FinalExPF_2@x)
  qqq <- colSums(Dmatrix)
  jjx_2 <- which((qqq-ppp)==0)
  if(length(jjx_2) > 0){
    FinalExPF_2 <- FinalExPF_2[,-jjx_2]
  }
  
  
  rm(FinalExPF)
  rm(eventsxpfam_path1)
  rm(eventsxpfam_path2)
  rm(ExTP1)
  rm(ExTP2)
  rm(transcriptnames)
  rm(Dmatrix)
  gc()
  
  results <- calculateCorrelationTest(A = t(FinalExPF_2),B = Diff_PSI,method = method)
  # results_pear <- calculateCorrelationTest(A = t(FinalExPF_2),B = Results_Table$deltaPSI,method = "pearson")
  
  #check results
  # m <- sample(1:ncol(FinalExPF_2),size=1)
  # mm <- match(colnames(FinalExPF_2)[m],rownames(results$mycor))
  # # cor.test(FinalExPF_2[,m],Results_Table$deltaPSI[,1],method = "pearson")
  # cor.test(FinalExPF_2[,m],Results_Table$deltaPSI,method = "spearman")
  # results$mycor[mm,1]
  # results$PVAL[mm,1]
  # results$STATISTIC[mm,1]
  
  return(results)
  
  
  
}

































