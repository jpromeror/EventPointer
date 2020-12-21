#' GetPSI_FromTranRef
#'
#' @description Get the values of PSI. A filer expression is applied if
#' the user select the option of filter.
#'
#' @param PathsxTranscript the output of EventDetection_transcriptome.
#' @param Samples  matrix or list containing the expression of the samples. 
#' @param Bootstrap Boolean variable to indicate if bootstrap data from pseudo-alignment is used.
#' @param Filter Boolean variable to indicate if an expression filter is applied. Default TRUE.
#' @param Qn Quartile used to filter the events (Bounded between 0-1, Qn would be 0.25 by default).
#'
#'
#' @return The output is a list containing two elements: a matrix with the
#' values of PSI  and a list containing as many matrices as number of events.
#' In each matrix is stored the expression of the different paths of an event along the samples.
#'
#' @examples
#'      data(EventXtrans)
#'      
#'      PathSamples <- system.file("extdata",package="EventPointer")
#'      PathSamples <- paste0(PathSamples,"/output")
#'      PathSamples <- dir(PathSamples,full.names = TRUE)
#'      
#'      data_exp <- getbootstrapdata(PathSamples = PathSamples,type = "kallisto")
#'      
#'      #same annotation
#'      rownames(data_exp[[1]]) <- gsub("\\|.*","",rownames(data_exp[[1]]))
#'      
#'      #Obtain values of PSI
#'      PSI_List <- GetPSI_FromTranRef(PathsxTranscript = EventXtrans,Samples = data_exp,Bootstrap = TRUE, Filter = FALSE)
#'      PSI <- PSI_List$PSI
#'      Expression_List <- PSI_List$ExpEvs
#'
#' @export
#' @importFrom matrixStats rowMins rowQuantiles



GetPSI_FromTranRef <- function(Samples,
                               PathsxTranscript,
                               Bootstrap = FALSE,
                               Filter = TRUE,
                               Qn = 0.25){
  
  
  
  
  
  if (is.null(PathsxTranscript)) {
    stop("PathsxTranscript field is empty")
  }
  if (is.null(Samples)) {
    stop("Samples field is empty")
  }
  
  ## check annotation of expression data and eventsxtranscripts can merge
  
  
  
  if(Bootstrap==FALSE){
    trannames_gtf <- PathsxTranscript$transcritnames
    trannames_Samples <- rownames(Samples)
    
    
    if (length(trannames_Samples) != length(trannames_gtf)) {
      warning("\nNot the same number of Transcripts in GTF and Fasta used...\n")
    }
    
    if (any(trannames_Samples %in% trannames_gtf == 
            FALSE)) {
      cat("\nRemoving transcripts from Samples data...")
      trannames_Samples <- trannames_Samples[trannames_Samples %in% trannames_gtf]
      if(length(trannames_Samples) == 0){
        stop("\nNo match between Samples and Events annotation\n")
      }
    }
    
    
    if (any(trannames_gtf %in% trannames_Samples == 
            FALSE)) {
      cat("\nRemoving Events that include transcripts not annotated in Samples data...\n")
      trannames_gtf <- trannames_gtf[trannames_gtf %in% trannames_Samples]
      if(length(trannames_gtf)==0){
        stop("\nNo match between Samples and Events annotation\n")
      }
      
      
    }
    
    
    Samples <- Samples[trannames_Samples,]
    jjx <- which(PathsxTranscript$transcritnames %in% trannames_gtf)
    Path1 <- PathsxTranscript$ExTP1[,jjx]
    Path2 <- PathsxTranscript$ExTP2[,jjx]
    PathRef <- PathsxTranscript$ExTPRef[,jjx]
    
    
    
    if (!identical(trannames_gtf, trannames_Samples)) {
      iix <- match(trannames_gtf, trannames_Samples)
      Samples <- Samples[iix, ]
    }
    # identical(rownames(Samples), trannames_gtf)
    
    # length(which(Matrix::rowSums(PathsxTranscript$ExTPRef)==0))
    rrx <- which(Matrix::rowSums(PathRef)==0)
    if(length(rrx)>0){
      Path1 <- Path1[-rrx,]
      Path2 <- Path2[-rrx,]
      PathRef <- PathRef[-rrx,]
    }
    
    ConcentrationPath1 <- as.matrix(Path1 %*% 
                                      Samples)
    ConcentrationPath2 <- as.matrix(Path2 %*% 
                                      Samples)
    ConcentrationPathRef <- as.matrix(PathRef %*% 
                                        Samples)
    
    
    if (Filter) {
      Min_p1 <- rowMins(ConcentrationPath1)
      Min_p2 <- rowMins(ConcentrationPath2)
      Min_pRef <- rowMins(ConcentrationPathRef)
      
      Max_p1 <- rowMaxs(ConcentrationPath1)
      Max_p2 <- rowMaxs(ConcentrationPath2)
      Max_pRef <- rowMaxs(ConcentrationPathRef)
      
      names(Min_p1) <- names(Min_p2) <- names(Max_pRef) <- names(Max_p1) <- names(Max_p2) <- names(Min_pRef) <- rownames(ConcentrationPath1)
      
      Qt_p1 <- rowQuantiles(ConcentrationPath1, 
                            probs = 0.8)
      Qt_p2 <- rowQuantiles(ConcentrationPath2, 
                            probs = 0.8)
      
      th_p <- quantile(c(Max_p1, Max_p2), 
                       Qn)
      th_pRef <- quantile(Max_pRef, Qn)
      
      Filt <- which((Min_pRef > th_pRef) & 
                      (Qt_p1 > th_p) & (Qt_p2 > th_p))
      
      
      ConcentrationPath1 <- ConcentrationPath1[Filt, 
                                               ]
      ConcentrationPath2 <- ConcentrationPath2[Filt, 
                                               ]
      ConcentrationPathRef <- ConcentrationPathRef[Filt, 
                                                   ]
      
      PSI <- ConcentrationPath1/ConcentrationPathRef
      
      ExpEvs <- vector(mode = "list", length = dim(ConcentrationPath1)[1])
      names(ExpEvs) <- rownames(ConcentrationPath1)
      
      matsamples <- matrix(0, nrow = dim(ConcentrationPath1)[2], 
                           ncol = 3)
      colnames(matsamples) <- c("P1", "P2", 
                                "PRef")
      rownames(matsamples) <- colnames(ConcentrationPath1)
      
      for (i in seq_len(dim(ConcentrationPath1)[1])) {
        matsamples[, 1] <- ConcentrationPath1[i, 
                                              ]
        matsamples[, 2] <- ConcentrationPath2[i, 
                                              ]
        matsamples[, 3] <- ConcentrationPathRef[i, 
                                                ]
        ExpEvs[[i]] <- matsamples
      }
      
      Result <- list(PSI = PSI, ExpEvs = ExpEvs)
      return(Result)
    } else {
      PSI <- ConcentrationPath1/ConcentrationPathRef
      
      ExpEvs <- vector(mode = "list", length = dim(ConcentrationPath1)[1])
      names(ExpEvs) <- rownames(ConcentrationPath1)
      
      matsamples <- matrix(0, nrow = dim(ConcentrationPath1)[2], 
                           ncol = 3)
      colnames(matsamples) <- c("P1", "P2", 
                                "PRef")
      rownames(matsamples) <- colnames(ConcentrationPath1)
      
      for (i in seq_len(dim(ConcentrationPath1)[1])) {
        matsamples[, 1] <- ConcentrationPath1[i, 
                                              ]
        matsamples[, 2] <- ConcentrationPath2[i, 
                                              ]
        matsamples[, 3] <- ConcentrationPathRef[i, 
                                                ]
        ExpEvs[[i]] <- matsamples
      }
      
      Result <- list(PSI = PSI, ExpEvs = ExpEvs)
      return(Result)
    }
    
  }
  else{
    
    trannames_gtf <- PathsxTranscript$transcritnames
    trannames_Samples <- rownames(Samples[[1]])
    
    
    if (length(trannames_Samples) != length(trannames_gtf)) {
      warning("\nNot the same number of Transcripts in GTF and Fasta references used...\n")
    }
    
    if (any(trannames_Samples %in% trannames_gtf == 
            FALSE)) {
      cat("\nRemoving transcripts from Samples data...")
      trannames_Samples <- trannames_Samples[trannames_Samples %in% trannames_gtf]
      if(length(trannames_Samples) == 0){
        stop("\nNo match between Samples and Events annotation\n")
      }
    }
    
    
    if (any(trannames_gtf %in% trannames_Samples == 
            FALSE)) {
      cat("\nRemoving Events that include transcripts not annotated in Samples data...\n")
      trannames_gtf <- trannames_gtf[trannames_gtf %in% trannames_Samples]
      if(length(trannames_gtf)==0){
        stop("\nNo match between Samples and Events annotation\n")
      }
      
      
    }
    
    if(!length(trannames_Samples) == nrow(Samples[[1]])){
      # Samples <- Samples[trannames_Samples,]
      aax <- which(rownames(Samples[[1]]) %in% trannames_Samples)
      Samples <- lapply(Samples,function(X){
        X[aax,]
      })
    }
    
    jjx <- which(PathsxTranscript$transcritnames %in% trannames_gtf)
    Path1 <- PathsxTranscript$ExTP1[,jjx]
    Path2 <- PathsxTranscript$ExTP2[,jjx]
    PathRef <- PathsxTranscript$ExTPRef[,jjx]
    
    
    
    if (!identical(trannames_Samples,trannames_gtf)) {
      iix <- match(trannames_Samples,trannames_gtf)
      Path1 <- Path1[,iix]
      Path2 <- Path2[,iix]
      PathRef <- PathRef[,iix]
      # identical(rownames(Samples[[1]]), trannames_gtf[iix])
    }
    
    
    # length(which(Matrix::rowSums(PathsxTranscript$ExTPRef)==0))
    rrx <- which(Matrix::rowSums(PathRef)==0)
    if(length(rrx)>0){
      Path1 <- Path1[-rrx,]
      Path2 <- Path2[-rrx,]
      PathRef <- PathRef[-rrx,]
    }
    
    ConcentrationPath1_b <- lapply(Samples,function(X){
      return(as.matrix(Path1%*%X))
    })
    ConcentrationPath2_b <- lapply(Samples,function(X){
      return(as.matrix(Path2%*%X))
    })
    ConcentrationPathRef_b <- lapply(Samples,function(X){
      return(as.matrix(PathRef%*%X))
    })
    
    #expression max. likelihood
    ConcentrationPath1 <- sapply(ConcentrationPath1_b, function(X){
      X[,1]
    })
    ConcentrationPath2 <- sapply(ConcentrationPath2_b, function(X){
      X[,1]
    })
    ConcentrationPathRef <- sapply(ConcentrationPathRef_b, function(X){
      X[,1]
    })
    
    
    
    if (Filter) {
      Min_p1 <- rowMins(ConcentrationPath1)
      Min_p2 <- rowMins(ConcentrationPath2)
      Min_pRef <- rowMins(ConcentrationPathRef)
      
      Max_p1 <- rowMaxs(ConcentrationPath1)
      Max_p2 <- rowMaxs(ConcentrationPath2)
      Max_pRef <- rowMaxs(ConcentrationPathRef)
      
      names(Min_p1) <- names(Min_p2) <- names(Max_pRef) <- names(Max_p1) <- names(Max_p2) <- names(Min_pRef) <- rownames(ConcentrationPath1)
      
      Qt_p1 <- rowQuantiles(ConcentrationPath1, 
                            probs = 0.8)
      Qt_p2 <- rowQuantiles(ConcentrationPath2, 
                            probs = 0.8)
      
      th_p <- quantile(c(Max_p1, Max_p2), 
                       Qn)
      th_pRef <- quantile(Max_pRef, Qn)
      
      Filt <- which((Min_pRef > th_pRef) & 
                      (Qt_p1 > th_p) & (Qt_p2 > th_p))
      
      
      PSI <- vector(mode = "list",length = length(Samples))
      names(PSI)<-names(Samples)
      for (kk in 1:length(Samples)){
        PSI[[kk]] <- ConcentrationPath1_b[[kk]][Filt,]/ConcentrationPathRef_b[[kk]][Filt,]
      }
      
      PSI <- array(unlist(PSI,use.names = FALSE), c(dim(PSI[[1]]),length(PSI)),
                   dimnames = list(rownames(PSI[[1]]),
                                   c(),
                                   names(PSI)))
      
      ConcentrationPath1 <- ConcentrationPath1[Filt, 
                                               ]
      ConcentrationPath2 <- ConcentrationPath2[Filt, 
                                               ]
      ConcentrationPathRef <- ConcentrationPathRef[Filt, 
                                                   ]
      
      ExpEvs <- vector(mode = "list", length = dim(ConcentrationPath1)[1])
      names(ExpEvs) <- rownames(ConcentrationPath1)
      
      matsamples <- matrix(0, nrow = dim(ConcentrationPath1)[2], 
                           ncol = 3)
      colnames(matsamples) <- c("P1", "P2", 
                                "PRef")
      rownames(matsamples) <- colnames(ConcentrationPath1)
      
      for (i in seq_len(dim(ConcentrationPath1)[1])) {
        matsamples[, 1] <- ConcentrationPath1[i, 
                                              ]
        matsamples[, 2] <- ConcentrationPath2[i, 
                                              ]
        matsamples[, 3] <- ConcentrationPathRef[i, 
                                                ]
        ExpEvs[[i]] <- matsamples
      }
      
      Result <- list(PSI = PSI, ExpEvs = ExpEvs)
      return(Result)
    } else {
      PSI <- vector(mode = "list",length = length(Samples))
      names(PSI)<-names(Samples)
      for (kk in 1:length(Samples)){
        PSI[[kk]] <- ConcentrationPath1_b[[kk]]/ConcentrationPathRef_b[[kk]]
      }
      
      PSI <- array(unlist(PSI,use.names = FALSE), c(dim(PSI[[1]]),length(PSI)),
                   dimnames = list(rownames(PSI[[1]]),
                                   c(),
                                   names(PSI)))
      
      ExpEvs <- vector(mode = "list", length = dim(ConcentrationPath1)[1])
      names(ExpEvs) <- rownames(ConcentrationPath1)
      
      matsamples <- matrix(0, nrow = dim(ConcentrationPath1)[2], 
                           ncol = 3)
      colnames(matsamples) <- c("P1", "P2", 
                                "PRef")
      rownames(matsamples) <- colnames(ConcentrationPath1)
      
      for (i in seq_len(dim(ConcentrationPath1)[1])) {
        matsamples[, 1] <- ConcentrationPath1[i, 
                                              ]
        matsamples[, 2] <- ConcentrationPath2[i, 
                                              ]
        matsamples[, 3] <- ConcentrationPathRef[i, 
                                                ]
        ExpEvs[[i]] <- matsamples
      }
      
      Result <- list(PSI = PSI, ExpEvs = ExpEvs)
      return(Result)
    }
    
  }
  
  
  
  
}