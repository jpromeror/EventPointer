#' Events X RBPS matrix creation
#'
#' Generates the Events x RBP matrix for the splicing factor enrichment analysis.
#'
#' @param pathtoeventstable Path to eventsFound.txt with the information of all the events
#' @param SG_List List with the information of the splicing graph of the genes. Returned by the funciotn EventDetectio_transcriptome
#' @param nt Number of nt up and down for the splicing regions of each event 
#' @param Peaks Table with the peaks 
#' @param POSTAR Table with peaks of POSTAR
#' @param EventsRegions Events regions if calculated prevously. Not need to calculated again.
#' @param cores Number of cores if user want to run in parallel.
#'
#' @return The function returns a list with the ExS matrix and with the splicing regions of the events. If the Splicign regions
#'     is an input of the function then only the ExS matrix will be returned. The ExS matrix is the input for the Splicing Factor
#'     enrichment analysis. 
#'
#'
#'
#'
#' @import GenomicRanges
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom GenomeInfoDb seqlevels `seqlevels<-`
#' @importFrom stats pt


CreateExSmatrix <- function(pathtoeventstable,
                            SG_List,
                            nt = 400,
                            Peaks,
                            POSTAR,
                            EventsRegions=NULL,
                            cores=1){
  
  
  #######################################################################
  ### get seq of each event
  #######################################################################
  # Each event has key regions where we would like to match the RBPs
  # The correct selection of this regions is cruzial to get good results. If we choose a little region PWMs won't
  # hit anyway, and on the contrary, if we select a huge region, PWMs will present too much hits.
  
  # for select this regions we first select the key points which will be the center of 
  # Intervals. Then we extend the interval from that points adding +- XXX nt
  
  # A cassette Event, for example, has de following structure.
  # The key points in this case are the points 1,2,3,4
  # REF5' -------------||||||||------------ REF3'
  # REF5' --------------------------------- REF3'
  #       1            2       3          4
  
  # For a complex event and in a general way
  # The key points in this case are the points 1 to 8
  # REF5' ---||--------||||||||-------------REF3'
  # REF5' ---||||---------------------|||    REF3'
  #       1  2  3      4      5       6 7  8
  
  #load data
  EventsFound <- read.delim(file=pathtoeventstable,stringsAsFactors = FALSE)
  EventsFound$EventID <- paste0(EventsFound$GeneName,"_",EventsFound$EventNumber)
  
  #type of events
  typeA <- c("Retained Intron","Cassette Exon","Alternative 3' Splice Site","Alternative 5' Splice Site")
  
  
  if(is.null(EventsRegions)){
    
    GRseq3 <- callGRseq_parallel(EventsFound,SG_List,cores,typeA,nt)
    
  } else{
    GRseq3 <- EventsRegions
  }
  
  
  
  #######################################################################
  ### Load clip peaks 
  #######################################################################
  
  
  POSTAR_L <- split.data.frame(x = POSTAR, f = POSTAR$RBP)
  Peaks_L <- split.data.frame(x = Peaks, f = Peaks$RBP)
  
  # table(names(POSTAR_L) %in% names(Peaks_L))
  # table(names(Peaks_L) %in% names(POSTAR_L))
  # we have repetead  
  
  mySF <- unique(c(names(POSTAR_L),names(Peaks_L)))
  mySF <- mySF[!mySF %in% c("AGO2MNASE", "HURMNASE", "PAN-ELAVL")]
  nSF <- length(mySF)
  ExS <- matrix(0, nrow = nrow(EventsFound), ncol = nSF)
  rownames(ExS) <- EventsFound$EventID
  colnames(ExS) <- mySF
  for(i in seq_len(nSF)){
    # i <- 1
    # i <- 5
    SF<- mySF[i]
    # cat(i, "/", nSF, " ",as.character(SF), "\n")
    jjx <- match(SF, names(POSTAR_L))
    iix <- match(SF, names(Peaks_L))
    if(!is.na(jjx)){
      peaks <- POSTAR_L[[jjx]]
      iD <- which(!(peaks$seqname %in% c(paste0("chr", c(1:22)), "chrX", "chrY")))
      if(length(iD) > 0){
        peaks <- peaks[-iD, ]
      }
      peaks_GR <- GRanges(peaks)
      seqlevels(peaks_GR) <- sub("chr", "", seqlevels(peaks_GR))
      peaks_GR<-reduce(peaks_GR)
      Overlaps <- findOverlaps(peaks_GR, GRseq3)
      
      EvMatch <- as.character(elementMetadata(GRseq3)$EventID[subjectHits(Overlaps)])
      if(length(EvMatch)>0){
        ExS[EvMatch,i] <- 1
      }
      
    }else if(!is.na(iix)){
      peaks <- Peaks_L[[iix]]
      iD <- which(!(peaks$seqname %in% c(paste0("chr", c(1:22)), "chrX", "chrY")))
      if(length(iD) > 0){
        peaks <- peaks[-iD, ]
      }
      peaks_GR <- GRanges(peaks)
      seqlevels(peaks_GR) <- sub("chr", "", seqlevels(peaks_GR))
      peaks_GR<-reduce(peaks_GR)
      Overlaps <- findOverlaps(peaks_GR, GRseq3)
      
      EvMatch <- as.character(elementMetadata(GRseq3)$EventID[subjectHits(Overlaps)])
      if(length(EvMatch)>0){
        ExS[EvMatch,i] <- 1
      }
    }
  }
  
  if(is.null(EventsRegions)){
    return(list(ExS=ExS,EventsRegions = GRseq3))    
  }else{
    return(ExS)    
  }
  
  
  
  
  
}














