#' Statistical analysis of alternative splcing events for RNASeq data
#'
#' Statistical analysis of all the alternative splicing events found in the given bam files.
#'
#' @param Events Output from EventDetection function
#' @param Design The design matrix for the experiment.
#' @param Contrast The contrast matrix for the experiment.
#' @param Statistic Statistical test to identify differential splicing events, must be one of : LogFC, Dif_LogFC and DRS.
#' @param PSI Boolean variable to indicate if PSI should be calculated for every splicing event.
#'
#'
#' @return Data.frame ordered by the splicing p.value . The object contains the different information for each splicing event
#' such as Gene name, event type, genomic position, p.value, z.value and delta PSI.
#'
#' @examples
#'
#'    data(AllEvents_RNASeq)
#'    Dmatrix<-matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1),ncol=2,byrow=FALSE)
#'    Cmatrix<-t(t(c(0,1)))
#'    Events <- EventPointer_RNASeq(AllEvents_RNASeq,Dmatrix,Cmatrix,Statistic='LogFC',PSI=TRUE)
#' @export

EventPointer_RNASeq <- function(Events, Design, 
    Contrast, Statistic = "LogFC", PSI = FALSE) {
    ### Screen output for users new output
    
    # stopifnot(Statistic == 'LogFC' |
    # Statistic == 'Dif_LogFC' | Statistic ==
    # 'DRS')
    
    if (is.null(Events)) {
        stop("Missing alternative splicing events")
    }
    options(warn = -1)
    
    if (Statistic == "LogFC") {
        stt <- "Logarithm of the fold change of both isoforms"
    } else if (Statistic == "Dif_LogFC") {
        stt <- "Relative concentrations of both isoforms"
    } else if (Statistic == "DRS") {
        stt <- "Difference of the logarithm of the fold change of both isoforms"
    } else {
        stop("Wrong statistical test provided")
    }
    
    if (PSI) {
        psi_m <- " Delta PSI will be calculated"
    }
    
    # if (classss(Design) != "matrix" | 
    #   classss(Contrast) != "matrix") {
    if (!is(Design,"matrix") | 
        !is(Contrast,"matrix")) {
        stop("Wrong Design and/or Contrast matrices")
    }
    
    TimeS <- paste(format(Sys.time(), "%X"), 
        sep = "")
    cat(paste(TimeS, " Running EventPointer: ", 
        sep = ""), "\n")
    cat(paste("\t** Statistical Analysis: ", 
        stt, sep = ""), "\n")
    
    if (PSI) {
        MPSI <- paste("\t**", psi_m, sep = "")
        cat(paste(MPSI, sep = ""), "\n")
    }
    
    cat(paste(paste(rep(" ", length(unlist(strsplit(TimeS, 
        "")))), sep = "", collapse = ""), 
        " ----------------------------------------------------------------", 
        sep = ""), "\n")
    
    ########################################################################################## 
    
    if (PSI == TRUE) {
        Msg <- paste("\t** Calculating PSI", 
            sep = "")
        cat(paste(Msg, "...", sep = ""))
        
        PSIss <- getPSI_RNASeq(Events)
        PSIs <- PSIss$PSI
        
        DPSIs <- vector("list", length = ncol(Contrast))
        
        fit <- lmFit(PSIs, design = Design)
        fit2 <- contrasts.fit(fit, Contrast)
        fit2 <- eBayes(fit2)
        
        for (jj in seq_len(ncol(Contrast))) {
            TopPSI <- topTable(fit2, coef = jj, 
                number = Inf)[, 1, drop = FALSE]
            DPSIs[[jj]] <- TopPSI
        }
    }
    
    
    
    Count_Matrix <- PrepareCountData(Events)
    
    # Auxiliary matrix for Kronecker Product
    
    if (Statistic == "LogFC" | Statistic == 
        "Dif_LogFC" | Statistic == "DRS") {
        AuxM <- matrix(c(1, 0, 0, 1, 1, 0, 
            1, 1, 1), nrow = 3, byrow = TRUE)
        
        D <- kronecker(Design, AuxM)
        
        # Limma Pipeline
        
        NormCounts <- voom(t(Count_Matrix), 
            D)
        
        fit <- lmFit(object = NormCounts, 
            design = D)
        
        FinalResult <- vector("list", length = ncol(Contrast))
        
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
                
                Final <- data.frame(Gen = rownames(T34_345), 
                  Pvalue = Values1$Pvalues, 
                  ZValue = Values1$Tstats, 
                  stringsAsFactors = FALSE)
                
                EventsN <- PrepareOutput(Events, 
                  Final)
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
                
                colnames(Final) <- c("Gene", 
                  "Pvalue", "Zvalue")
                
                EventsN <- PrepareOutput(Events, 
                  Final)
            }
            
            # Add extra information (Gene Name and
            # Event Classification) and Sort
            # data.frame by pvalue
            
            
            if (PSI) {
                IIx <- match(rownames(EventsN), 
                  rownames(DPSIs[[mm]]))
                EventsN <- cbind(EventsN, 
                  DPSIs[[mm]][IIx, ])
                colnames(EventsN)[6] <- "Delta PSI"
            }
            
            FinalResult[[mm]] <- EventsN
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
