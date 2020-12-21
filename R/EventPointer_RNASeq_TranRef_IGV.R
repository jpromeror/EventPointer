#' EventPointer RNASeq from reference transcriptome IGV Visualization
#'
#' @description Generates of files to be loaded in IGV for visualization and interpretation of events detected
#' from a reference transcriptome (see EventDetection_transcriptome).
#'
#' @param SG_List List with the Splicing Graph information of the events. This list is created by EventDetection_transcriptome function.
#' @param pathtoeventstable Complete path to the table returned by EventDetection_transcriptome that contains the information of each event,
#'  or table with specific events that the user want to load into IGV to visualize.
#' @param PathGTF Directory where to write the GTF files.
#'
#'
#' @return The function displays a progress bar to show the user the progress of the function. Once the progress bar reaches 100%, one .gtf
#' file is written to the specified directory in PathGTF. The created file is named 'paths_RNASeq.gtf'. 
#'
#' @examples
#'
#'   ###### example using all the events found in a reference transcriptome
#'   data("EventXtrans")
#'   SG_List <- EventXtrans$SG_List
#'   PathEventsTxt<-system.file('extdata',package='EventPointer')
#'   PathEventsTxt <- paste0(PathEventsTxt,"/EventsFound_Gencode24_2genes.txt")
#'   PathGTF <- tempdir()
#'   
#'   EventPointer_RNASeq_TranRef_IGV(SG_List = SG_List,pathtoeventstable = PathEventsTxt,PathGTF = PathGTF)
#'   
#'   
#'      
#' @export
#' @import Matrix
#' @importFrom RBGL connectedComp
#' @importFrom  graph ftM2graphNEL
#' @importFrom  prodlim row.match



EventPointer_RNASeq_TranRef_IGV <- function(SG_List,pathtoeventstable,PathGTF){
  
  if(is.null(SG_List)){
    stop("SG_List field is empty")
  }
  if(is.null(pathtoeventstable)){
    stop("pathtoeventstable field is empty")
  }
  if(is.null(PathGTF)){
    stop("PathGTF field is empty")
  }
  
  Result2 <- read.delim(file=pathtoeventstable,stringsAsFactors = FALSE)
  
  Result2$GeneName <- paste0(Result2$GeneName,"_",Result2$EventNumber)
  
  cat("\nCreating .GTF ...")
  # Create file to store gtf for patths
  # (events)
  FILE.paths <- paste(PathGTF, "/paths_RNASeq.gtf", 
                      sep = "")
  cat(file = FILE.paths, paste("#track name=", 
                               shQuote("paths", type = "cmd"), " gffTags=", 
                               shQuote("on", type = "cmd"), sep = ""), 
      "\n")
  
  pb <- txtProgressBar(min = 0, max = nrow(Result2), 
                       style = 3)
  
  # 1:nrow(Result2)
  for (jj in seq_len(nrow(Result2))) {
    setTxtProgressBar(pb, jj)
    
    Gene <- Result2$GeneName[jj]
    Gene <- gsub("_.*","",Gene)
    EvtEdg <- SG_List[[Gene]]$Edges
    
    if (unique(EvtEdg[, "Strand"]) == 
        "") {
      EvtEdg[, "Strand"] <- "-"
    }
    
    # iixx <- which(Result2$GeneName==Gene)
    
    EventPaths <- GetIGVPaths(Result2[jj, 
                                      ], EvtEdg)
    class(EventPaths[, 2]) <- "integer"
    class(EventPaths[, 3]) <- "integer"
    WriteGTF_RNASeq(PathGTF, Result2[jj, 
                                     ], EventPaths)
  }
  
  close(pb)
  
  cat("\n.GTF created")
  
  
  
}


