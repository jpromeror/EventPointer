# FPKMs/Counts associated with an alternative splicing event
#
#' Internal function called by EvenPointer
#'
#' Function to obtain either counts or fpkms for each of the three paths of an alternative splicing event
#'
#' @keywords internal
#'
#' @param Events Detected alternative splicing events for a particular gene
#' @param sg_txiki Splicing graph of the gene
#' @param type Either counts of fpkms
#'
#'
#' @return matrix with Counts/FPKMs for alternative splicing events
#'

GetCounts <- function(Events,sg_txiki, type = "counts") {
  readsC <- counts(sg_txiki)
  readsF <- FPKM(sg_txiki)
  countsEvents <- lapply(Events, getPathCounts,readsC)
  countsEvents<-lapply(countsEvents, getPathFPKMs,readsF)
  return(countsEvents)
}

getPathCounts <- function(x, readsC, widthinit) {
  reads <- rbind(colSums(readsC[x$P1$featureID,,drop = FALSE]),
                 colSums(readsC[x$P2$featureID,,drop = FALSE]),
                 colSums(readsC[x$Ref$featureID,,drop = FALSE]))
  
  rownames(reads) <- c("P1","P2","Ref")
  x$Counts<-reads
  return(x)
}

getPathFPKMs <- function(x, readsC, widthinit) {
  reads <- rbind(colSums(readsC[x$P1$featureID,,drop = FALSE]),
                 colSums(readsC[x$P2$featureID,,drop = FALSE]),
                 colSums(readsC[x$Ref$featureID,,drop = FALSE]))
  
  rownames(reads) <- c("P1","P2","Ref")
  x$FPKM<-reads
  return(x)
}
