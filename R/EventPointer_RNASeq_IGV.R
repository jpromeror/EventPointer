#' EventPointer RNASeq IGV Visualization
#'
#' Generates files to be loaded in IGV for visualization and interpretation of events
#'
#' @param SG_RNASeq Output from EventDetection_BAM function. Contains splicing graphs components.
#' @param EventsCSV Path to EventsFound.txt file generated with EventDetection_BAM function
#' @param PathGTF Directory where to write the GTF files.
#'
#'
#' @return Creates paths_RNASeq.gtf file that represents the alternative splicing events.
#'
#' @examples
#' \dontrun{
#'   data(AllEvents_RNASeq)
#'   data(SG_RNASeq)
#'
#'    # Run EventPointer
#'
#'    Dmatrix<-matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1),ncol=2,byrow=FALSE)
#'    Cmatrix<-t(t(c(0,1)))
#'    Events <- EventPointer_RNASeq(AllEvents_RNASeq,Dmatrix,Cmatrix,Statistic='LogFC',PSI=TRUE)
#'
#'    # IGV Visualization
#'
#'    EventsCSV<-paste(system.file('extdata',package='EventPointer'),'/EventsFound_RNASeq.txt',sep='')
#'    PathGTF<-tempdir()
#'    EventPointerBAM_IGV(Events,SG_RNASeq,EventsCSV,PathGTF)
#'    }
#' @export

EventPointerBAM_IGV <- function(SG_RNASeq, EventsCSV, PathGTF) {
    
    if (is.null(SG_RNASeq)) {
        stop("Missing splicing graphs information")
    }
  if (is.null(EventsCSV) | 
      !is(EventsCSV,"character")) {
    stop("Wrong or missing EventsCSV field")
    }
    
    # if (is.null(PathGTF) | 
    #  classsss(PathGTF) != "character") {
    if (is.null(PathGTF) | 
        !is(PathGTF,"character")) {
        stop("Wrong or missing PathGTF field")
    }
    
    
    
    SgF <- rowRanges(SG_RNASeq)
    
    ######################## iiP<-which(strand(SgF)@values=='+')
    ######################## iiN<-which(strand(SgF)@values=='-')
    ######################## strand(SgF)@values[iiP]<-'-'
    ######################## strand(SgF)@values[iiN]<-'+'
    
    # Create file to store gtf for patths
    # (events)
    FILE.paths <- paste(PathGTF, "/paths_RNASeq.gtf", 
        sep = "")
    cat(file = FILE.paths, paste("#track name=", 
        shQuote("paths", type = "cmd"), " gffTags=", 
        shQuote("on", type = "cmd"), sep = ""), 
        "\n")
    
    # Read EventsFound txt
  EventsInfo <- read.delim(file = EventsCSV, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    
    cat("\n Generating GTF Files...")
    
  for (jj in seq_len(nrow(EventsInfo))) {
    
    EventXX <- as.numeric(unlist(strsplit(EventsInfo$EventID[jj], "_")))
        Gene <- EventXX[1]
        EvNumb <- EventXX[2]
        # Gene<-EventsInfo[jj,1]
        
        SG_Gene <- SgF[unlist(geneID(SgF)) == 
            Gene]
        
        SG_Edges <- SG_Info(SG_Gene)$Edges
        
        if (unique(SG_Edges[, "Strand"]) == 
            "") {
            SG_Edges[, "Strand"] <- "-"
        }
        
    iixx <- jj
        
        EventPaths <- GetIGVPaths(EventsInfo[iixx, 
            ], SG_Edges)
        class(EventPaths[, 2]) <- "integer"
        class(EventPaths[, 3]) <- "integer"
        WriteGTF_RNASeq(PathGTF, EventsInfo[iixx, 
            ], EventPaths)
    }
    
    cat("\n")
}
