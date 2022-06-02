#' Events_ReClassification
#' 
#' EventPointer can detect splicing events that cannot
#'  be cataloged in any of the canonical types 
#'  (Cassette Exon, Alternative 3' or 5' splice site, 
#'  retained intron and mutually exclusive exon). 
#'  These events are classified as "Complex Events". 
#'  With this function, EventPointer reclassifies these
#'  complex events according to how similar the event 
#'  is to the canonical events. 
#'  The same complex event can have several types.
#' Further, EP adds a new type of event: "multiple skipping exon". 
#' These events are characterized by presenting several exons in a 
#' row as alternative exons. If there is only one alternative exon
#'  we would be talking about a "Casstte Exon".
#' 
#' 
#'
#'
#' @param EventTable Table returned by EventDetection_transcriptome. 
#' Can be easily loaded using the function read.delim as data.frame.
#' @param SplicingGraph A list with the splicing graph of all the genes of
#'  a reference transcriptome. This data is returned
#'   by the function EventDetection_transcriptome.
#'
#' @return A data.frame containing a new column with the new classification ('EventType_new'):
#' 
#'
#' @examples
#'
#'    #load splicing graph
#' data("SG_reclassify")
#' 
#' #load table with info of the events
#' PathFiles<-system.file("extdata",package="EventPointer")
#' inputFile <- paste(PathFiles,"/Events_found_class.txt",sep="")
#' EventTable <- read.delim(file=inputFile)
#' #this table has the information of 5 complex events.
#' 
#' EventTable_new <- Events_ReClassification(EventTable = EventTable,
#'                                           SplicingGraph = SG_reclassify)
#'          
#'     
#'
#' @export
#' @import Matrix
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar combn write.table read.table
#' @importFrom stringr str_count
#' @importFrom igraph graph_from_data_frame as_adj clusters graph_from_adjacency_matrix graph.data.frame
#' @importFrom MASS Null ginv
#' @importFrom nnls nnls
#' @importFrom methods as new slot 'slot<-'
#' @importFrom graph ftM2graphNEL
#' @importFrom prodlim row.match
#' @importFrom 'methods' is
#' @importFrom S4Vectors queryHits subjectHits split
#' @importFrom IRanges relist disjoin %over% reduce

Events_ReClassification <- function(EventTable,SplicingGraph){
  
  if (is.null(EventTable)) {
    stop("not EventTable")
  }
  if (is.null(SplicingGraph)) {
    stop("not SplicingGraph")
  }
  
  EventTable$ID <- paste0(EventTable$GeneName,"_",EventTable$EventNumber)
  EventTable$EventType_new <- EventTable$EventType
  
  uux <- which(EventTable$EventType=="Complex Event")
  if(length(uux)==0){
    return(EventTable)
  }
  cat("\nStarting reclassification\n")
  pb <- txtProgressBar(min = 0, max = length(uux), 
                       style = 3)
  misnewtype <- c()
    for(jj in seq_len(length(uux))){
      setTxtProgressBar(pb, jj)
      # cat(jj,"\n")
      mievento <- EventTable$ID[uux][jj]
      # mievento
      ggx <- match(EventTable$GeneName[uux][jj],names(SplicingGraph))
      # ggx
      SG <- SplicingGraph[[ggx]]
      
      pp1 <- EventTable$Path.1[uux][jj]
      pp2 <- EventTable$Path.2[uux][jj]
      ppref <- EventTable$Path.Reference[uux][jj]
      misnewtype <- c(misnewtype,reclasify_intern(SG = SG,mievento = mievento,
                                                  pp1=pp1,pp2=pp2,ppref=ppref))
      
    }
  
  close(pb)
  cat("\nReclassification Finished\n")
  # Result <- unlist(Result)
  # EventTable$EventType_new[uux] <- Result
  
  EventTable$EventType_new[uux] <- misnewtype
  
  return(EventTable)
  
}