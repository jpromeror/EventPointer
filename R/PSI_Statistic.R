#' PSI_Statistic
#'
#' Statistical analysis of the alternative splicing events. This function takes as input the values of PSI.
#' Perform a statistical analysis based on permutation test
#'
#' @param PSI A matrix with the values of the PSI.
#' @param Design The design matrix for the experiment.
#' @param Contrast The contrast matrix for the experiment.
#' @param nboot The number of random analysis.
#'
#' @return The output of these functions is a list containing: two data.frame (deltaPSI and Pvalues) with the values of the deltaPSI
#'  and the p.values for each contrast, and a third element (LocalFDR) with the information of the local false discovery rate.
#'
#' @examples
#'       data(ArraysData)
#'       PSI_Arrays_list<-EventPointer:::getPSI(ArraysData)
#'       PSI_Arrays <- PSI_Arrays_list$PSI
#'       Design <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
#'       Contrast <- matrix(c(0,1),nrow=1)
#'       
#'       # Statistical analysis:
#'       
#'       table <- PSI_Statistic(PSI_Arrays,Design = Design, Contrast = Contrast, nboot = 50)
#'       
#'
#' @export
#' @importFrom IRanges relist disjoin %over% reduce
#' @importFrom stats pbeta predict
#' @importFrom qvalue qvalue
#' @importFrom cobs cobs
#' @importFrom matrixStats rowVars






PSI_Statistic <- function(PSI, Design, Contrast, 
    nboot) {
    if (is.null(PSI)) {
        stop("PSI field is empty")
    }
    
    # if (classsss(Design) != "matrix" | 
    #   classsss(Contrast) != "matrix") {
    if (!is(Design,"matrix") | 
        !is(Contrast,"matrix")) {
        stop("Wrong Design and/or Contrast matrices")
    }
    
    if (is.null(nboot)) {
        stop("nboot field is empty")
    }
    
    eps <- 10 * .Machine$double.eps
    V <- Contrast %*% (solve(crossprod(Design)) %*% 
        t(Design))
    incrPSI_original <- (tcrossprod(V, PSI) + 
        1)/2
    ncontrastes <- dim(Contrast)[1]
    l <- dim(PSI)[1]
    combboots <- rep(list(matrix(NA, l, nboot)), 
        ncontrastes)  # Intialize matrix for the increase in PSI
    
    for (boot2 in seq_len(nboot)) {
        PSI1 <- PSI[, sample(ncol(PSI), replace = TRUE)]  # Samples the Yb (mixes data)
        output <- tcrossprod(V, PSI1)  # Obtain the increase in PSI
        for (boot3 in seq_len(ncontrastes)) {
            combboots[[boot3]][, boot2] <- output[boot3, 
                ]  # Fills matrix of increase in PSI
        }
    }
    
    table <- get_beta(combboots, incrPSI_original, 
        ncontrastes)
    Pvalues <- table$pvalues
    
    # return(table) estimation of the lfdr.
    # We need to remove the NaN
    if (nrow(Contrast) > 1) {
        LocalFDR <- apply(Pvalues, 2, function(X) {
            result <- qvalue(X, trunc = FALSE, 
                monotone = FALSE)
            salida <- cobs(x = log(result$pvalues + 
                eps), y = result$lfdr, constraint = "incr", 
                pointwise = matrix(c(-1, 
                  0, 1, 1, min(log(result$pvalues + 
                    eps), na.rm = TRUE), 
                  0), byrow = TRUE, ncol = 3), 
                lambda = 1, print.mesg = FALSE)
            iix <- which(is.na(X))
            if (length(iix) > 0) {
                jjx <- which(!is.na(X))
                if (length(jjx) > 0) {
                  result$lfdr2 <- result$lfdr
                  result$lfdr2[jjx] <- predict(salida, 
                    log(X[jjx] + eps))[, 
                    2]
                } else {
                  Prueba$lfdr2 <- Prueba$lfdr
                }
            } else {
                result$lfdr2 <- predict(salida, 
                  log(X + eps))[, 2]
            }
            return(result)
        })
    } else {
        LocalFDR <- qvalue(Pvalues, trunc = FALSE, 
            monotone = FALSE)
        salida <- cobs(x = log(LocalFDR$pvalues + 
            eps), y = LocalFDR$lfdr, constraint = "incr", 
            pointwise = matrix(c(-1, 0, 1, 
                1, min(log(LocalFDR$pvalues + 
                  eps), na.rm = TRUE), 0), 
                byrow = TRUE, ncol = 3), 
            lambda = 1, print.mesg = FALSE)
        iix <- which(is.na(Pvalues))
        if (length(iix) > 0) {
            jjx <- which(!is.na(Pvalues))
            if (length(jjx) > 0) {
                LocalFDR$lfdr2 <- LocalFDR$lfdr
                LocalFDR$lfdr2[jjx] <- predict(salida, 
                  log(Pvalues[jjx] + eps))[, 
                  2]
            } else {
                LocalFDR$lfdr2 <- LocalFDR$lfdr
            }
        } else {
            LocalFDR$lfdr2 <- predict(salida, 
                log(Pvalues + eps))[, 2]
        }
    }
    
    tablefinal <- list(deltaPSI = table$deltaPSI, 
        Pvalues = Pvalues, LocalFDR = LocalFDR)
    
    return(tablefinal)
}
