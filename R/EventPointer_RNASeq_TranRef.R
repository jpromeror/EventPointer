#' EventPointer_RNASeq_TranRef
#' 
#' Statistical analysis of alternative splicing events with the output of GetPSI_FromTranRef
#'
#' @param Count_Matrix The list containing the expression data taken from
#' the ouput of GetPSI_FromTranRef
#' @param Statistic The type of statistic to apply.
#' Default = 'LogFC' (can be 'logFC, 'Dif_LogFC','DRS')
#' @param Design The design matrix of the experiment.
#' @param Contrast The Contrast matrix of the experiment.
#'
#'
#' @return a data.frame with the information of the names of the event,
#' its p.values and the corresponding z.value. If there is more than
#' one contrast, the function returns as many data.frames as number
#' of contrast and all these data.frame are sotred in an unique list.
#'
#'
#' @examples
#'    data(EventXtrans)
#'    data(PSIss)
#'    # Design and contrast matrix:
#'
#'    Design <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
#'    Contrast <- matrix(c(0,1),nrow=2)
#'
#'    # Statistical analysis:
#'
#'
#'    Fit <- EventPointer_RNASeq_TranRef(Count_Matrix =  PSIss$ExpEvs,
#'                                       Statistic = 'LogFC',Design = Design,
#'                                       Contrast = Contrast)
#'
#'
#' @export




EventPointer_RNASeq_TranRef <- function(Count_Matrix, 
    Statistic = "LogFC", Design, Contrast) {
    # if (classsss(Count_Matrix) != "list") {
    if (!is(Count_Matrix,"list")) {
        stop("Wrong Count_Matrix input: must be a list")
    }
    
    if (Statistic == "LogFC") {
        stt <- "Logarithm of the fold change of both isoforms"
    } else if (Statistic == "Dif_LogFC") {
        stt <- "Relative concentrations of both isoforms"
    } else if (Statistic == "DRS") {
        stt <- "Difference of the logarithm of the fold change of both isoforms"
    } else {
        stop("Wrong statistical test given")
    }
    
    # if (classsss(Design) != "matrix" 
    # | classsss(Contrast) != "matrix") {
    if ( !is(Design,"matrix") | 
         !is(Contrast,"matrix")) {
        stop("Wrong Design and/or Contrast matrices")
    }
    
    
    
    if (Statistic == "LogFC" | Statistic == 
        "Dif_LogFC" | Statistic == "DRS") {
        AuxM <- matrix(c(1, 0, 0, 1, 1, 0, 
            1, 1, 1), nrow = 3, byrow = TRUE)
        
        D <- kronecker(Design, AuxM)
        
        Count_Matrix <- sapply(Count_Matrix, 
            function(X) return(t(X[, c(3, 
                1, 2)])))
        
        
        # Limma Pipeline
        
        NormCounts <- voom(t(Count_Matrix), 
            D)
        
        fit <- lmFit(object = NormCounts, 
            design = D)
        
        FinalResult <- vector("list", length = ncol(Contrast))
        names(FinalResult) <- colnames(Contrast)
        
        for (mm in seq_len(ncol(Contrast))) {
            Cused <- Contrast[, mm, drop = FALSE]
            
            # The contrasts we are interested in are
            # the ones related with each Path, and we
            # apply a kronecker product of the
            # contrast matrix with the corresponding
            # vector for each Path (P1 = 1 1 0 ; P2 =
            # 1 1 1)
            
            if (Statistic == "LogFC" | Statistic == 
                "Dif_LogFC") {
                if (Statistic == "LogFC") {
                  P1 <- kronecker(Cused, 
                    matrix(c(1, 1, 0), nrow = 3))
                  P2 <- kronecker(Cused, 
                    matrix(c(1, 1, 1), nrow = 3))
                } else if (Statistic == "Dif_LogFC") {
                  P1 <- kronecker(Cused, 
                    matrix(c(0, 1, 0), nrow = 3))
                  P2 <- kronecker(Cused, 
                    matrix(c(0, 1, 1), nrow = 3))
                }
                
                
                C <- cbind(P1, P2)
                
                fit2 <- contrasts.fit(fit, 
                  C)
                
                fit2 <- eBayes(fit2)
                
                # Merge the results from both contrasts
                # in one table
                
                T2 <- topTable(fit2, coef = 1, 
                  number = Inf)
                T3 <- topTable(fit2, coef = 2, 
                  number = Inf)
                
                EvsIds <- rownames(T2)
                ii3 <- match(EvsIds, rownames(T3))
                T3 <- T3[ii3, ]
                
                colnames(T3) <- letters[seq_len(ncol(T3))]
                T34_345 <- cbind(T2, T3)
                
                # Irwin Hall Pvalue Summarization
                Values1 <- IHsummarization(T34_345[, 
                  4], T34_345[, 3], T34_345[, 
                  10], T34_345[, 9])
                
                Final <- data.frame(Event_ID = rownames(T34_345), 
                  Pvalue = Values1$Pvalues, 
                  Zvalue = Values1$Tstats, 
                  stringsAsFactors = FALSE)
                
                rownames(Final) <- Final$Event_ID
                
                # EventsN <- PrepareOutput(Events, Final)
            } else if (Statistic == "DRS") {
                DRS <- kronecker(Cused, matrix(c(0, 
                  0, 1), nrow = 3))
                
                # Compute estimated coefficients and
                # standard errors for the given contrasts
                fit2 <- contrasts.fit(fit, 
                  DRS)
                
                # Empirical Bayesian Statistics
                fit2 <- eBayes(fit2)
                
                # Obtain the ranking of events for each
                # of the contrasts
                T2 <- topTable(fit2, number = Inf)
                
                Final <- data.frame(rownames(T2), 
                  T2[, 4], T2[, 3], stringsAsFactors = FALSE)
                
                colnames(Final) <- c("Event_ID", 
                  "Pvalue", "Zvalue")
                
                rownames(Final) <- Final$Event_ID
                
                # EventsN <- PrepareOutput(Events, Final)
            }
            
            # Add extra information (Gene Name and
            # Event Classification) and Sort
            # data.frame by pvalue
            
            
            FinalResult[[mm]] <- Final
        }
        
        if (ncol(Contrast) == 1) {
            FinalResult <- FinalResult[[1]]
        }
        
        cat("Done")
        
        cat("\n Analysis Finished")
        
        cat(paste("\n Done \n", sep = ""))
        
        # Return the Result to the user
        cat("\n", format(Sys.time(), "%X"), 
            " Analysis Completed \n")
        return(FinalResult)
    } else {
        stop("Wrong Statistical Analysis Given")
    }
}
