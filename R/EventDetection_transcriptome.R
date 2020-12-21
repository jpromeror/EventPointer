#' EventDetection_transcriptome
#' 
#' Finds all the possible alternative splicing (AS) events given a reference transcriptome. This function use parallel foreach. User must
#' set the value of cores (by default equal to one). Moreover, it will create a .txt file with the relative information of all the
#' AS events found. Besides, it will return a list with main information of the splicing graph of each event. This list will be used as an
#' input in downstream functions (Get_PSI_FromTranRef, FindPrimers, and EventPointer_RNASeq_TranRef_IGV)
#'
#'
#' @param inputFile Path to the GTF file of the reference transcriptome.
#' @param Transcriptome Name of the transcriptome
#' @param Pathtxt Directory to save the .txt of the events found
#' @param cores Number of cores using in the parallel processing (by default = 1)
#'
#' @return a list is returned with the following information:
#' 
#' ExTP1  a sparce matrix of Events x Transcripts that relates which isoform build up the path1 of each event.
#' 
#' ExTP2  a sparce matrix of Events x Transcripts that relates which isoform build up the path2 of each event.
#' 
#' ExTPRef  a sparce matrix of Events x Transcripts that relates which isoform build up the pathRef of each event.
#' 
#' transcritnames  a vector with the annotation names of the isoforms.
#' 
#' SG_List  A list containing the information of the splicing graph of each gene.
#'
#' @examples
#'
#'          PathFiles<-system.file("extdata",package="EventPointer")
#'          inputFile <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep="")
#'          Transcriptome <- "Gencode24_2genes"
#'          Pathtxt <- tempdir()
#'          
#'          
#'          # Run the function 
#'          
#'          EventXtrans <- EventDetection_transcriptome(inputFile = inputFile,
#'                                                       Transcriptome = Transcriptome,
#'                                                       Pathtxt=Pathtxt,
#'                                                       cores=1)
#'
#' @export
#' @import Matrix
#' @import SGSeq
#' @import SummarizedExperiment
#' @importFrom affxparser writeCdf
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom GenomicFeatures makeTxDbFromBiomart makeTxDbFromUCSC makeTxDbFromGFF
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar combn write.table read.table
#' @importFrom stringr str_count
#' @importFrom GenomeInfoDb 'seqlevelsStyle<-' seqlevelsStyle seqnames
#' @importFrom igraph graph_from_data_frame as_adj clusters graph_from_adjacency_matrix graph.data.frame
#' @importFrom MASS Null ginv
#' @importFrom stats dist qnorm quantile runif setNames
#' @importFrom nnls nnls
#' @importFrom limma lmFit contrasts.fit eBayes topTable voom
#' @importFrom matrixStats rowMaxs
#' @importFrom RBGL connectedComp
#' @importFrom methods as new slot 'slot<-'
#' @importFrom graph ftM2graphNEL
#' @importFrom prodlim row.match
#' @importFrom 'methods' is
#' @importFrom GenomicRanges makeGRangesFromDataFrame granges
#' @importFrom S4Vectors queryHits subjectHits split
#' @importFrom IRanges relist disjoin %over% reduce




EventDetection_transcriptome <- function(inputFile = NULL, 
                                         Transcriptome = NULL,
                                         Pathtxt = NULL,
                                         cores = 1){
  if (is.null(inputFile)) {
    stop("not inputFile")
  }
  if (is.null(Transcriptome)) {
    stop("not Transcriptome")
  }
  if (is.null(Pathtxt)) {
    stop("not Pathtxt")
  }
  
  
  cat("Creating SG Information...")
  
  
  TxDb <- makeTxDbFromGFF(file = inputFile, 
                          format = "gtf", dataSource = "Custom GTF")
  TranscriptFeatures <- convertToTxFeatures(TxDb)
  
  transcritnames <- txName(TranscriptFeatures)
  transcritnames <- unique(unlist(transcritnames))
  
  genes22 <- unique(unlist(geneName(TranscriptFeatures)))
  numberofgenes <- length(genes22)
  
  Result2 <- vector("list", length = numberofgenes)
  
  SG_List <- vector("list", length = numberofgenes)
  names(SG_List) <- genes22
  
  pb <- txtProgressBar(min = 0, max = numberofgenes, 
                       style = 3)
  
  registerDoParallel(cores = cores)
  Result <- foreach(jj = seq_len(numberofgenes)) %dopar%
    {
      # jj <- 1
      setTxtProgressBar(pb, jj)
      Gene <- genes22[jj]
      SG_Gene <- TranscriptFeatures[which(any(geneName(TranscriptFeatures) == Gene)), ]
      
      if(length(unique(type(SG_Gene)))==1){
        # next()  #if there is only unique exons in transcripts
        return(list(Result2 = NULL,SG_List=NULL))
      } else{
        
        if(length(unique(unlist(txName(SG_Gene)))) == 1){
          # next() #if there is only one transcript
          return(list(Result2 = NULL,SG_List=NULL))
        }else{
          
          SG <- try(SG_creation_fast(SG_Gene),silent=TRUE)
          # if(classsss(SG)=="try-error"){
          if(is(SG,"try-error")){
            return(list(Result2 = NULL,SG_List=NULL))
          }else{
            # SG_List[[jj]] <- SG
            # plot(ftM2graphNEL(as.matrix(SG$Edges[,1:2]),edgemode='directed'))
            
            randSol <- getRandomFlow(SG$Incidence, 
                                     ncol = 10)
            
            if (any(round(randSol) != 1)) {
              
              Events <- findTriplets(randSol)
              
              if (nrow(Events$triplets) > 
                  0) {
                twopaths <- which(rowSums(Events$triplets != 
                                            0) == 3)
                
                Events <- getEventPaths(Events, 
                                        SG)
                
                if (length(Events) > 0) {
                  
                  Events <- ClassifyEvents(SG, 
                                           Events, twopaths)
                  
                  GenI <- jj
                  Events <- AnnotateEvents_KLL(Events, 
                                               Gene, GenI)
                  
                  if (is.null(Events)) {
                    # Result2[[jj]] <- Events
                    return(list(Result2 = NULL,SG_List = SG))
                  } else {
                    
                    transcrits <- sacartranscritos(edgetr = SG$Edges, 
                                                   Events)
                    
                    Events$tran_P1 <- ""
                    Events$tran_P2 <- ""
                    Events$tran_Ref <- ""
                    
                    Events$tran_P1 <- as.vector(transcrits$p1)
                    Events$tran_P2 <- as.vector(transcrits$p2)
                    Events$tran_Ref <- as.vector(transcrits$ref)
                    
                    # Result2[[jj]] <- Events
                    return(list(Result2 = Events,SG_List = SG))
                  }
                  
                }
                return(list(Result2 = NULL,SG_List=NULL))
              }
              return(list(Result2 = NULL,SG_List=NULL))
            }
            return(list(Result2 = NULL,SG_List=NULL))
            
            
            
          }
        }  
      }
      
      
    }
  
  close(pb)
  
  cat("\nCreating .txt ...")
  
  Result2 <- lapply(Result,function(X){
    X$Result2
  })
  SG_List <- lapply(Result,function(X){
    X$SG_List
  })
  names(SG_List) <- genes22
  
  
  Result2 <- Filter(Negate(is.null), Result2)
  Result2 <- do.call(rbind, Result2)
  
  goodone2 <- comprobaciontranscritos2(Result2)
  
  Result2 <- Result2[goodone2, ]
  
  
  colnames(Result2) <- c("GeneName", "GeneID", 
                         "EventNumber", "EventType", "GPos", 
                         "Path.1", "Path.2", "Path.Reference", 
                         "tran_P1", "tran_P2", "tran_Ref")
  
  write.table(Result2, file = paste(Pathtxt, 
                                    "/EventsFound_", Transcriptome, ".txt", 
                                    sep = ""), sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
  
  cat("\n.txt created")
  
  
  # paths x transcipts matrices
  
  cat("\nCreating the sparseMatrix of paths x transcripts...")
  
  rnames <- paste(Result2$GeneName, Result2$EventNumber, 
                  sep = "_")
  
  
  ListaNombres <- strsplit(Result2$tran_P1, 
                           "\\|")
  nTperPath <- sapply(ListaNombres, length)
  # TrNames <- unique(unlist(ListaNombres))
  # # TrNames includes all the possible
  # transcripts.  In transcritnames is
  # sotored the names of all the Isoforms
  TrInPaths <- unlist(ListaNombres)
  i <- rep(seq_len(length(ListaNombres)), 
           nTperPath)
  j <- match(TrInPaths, transcritnames)
  ExTP1 <- sparseMatrix(i, j, x = 1, dims = c(max(i), 
                                              length(transcritnames)))
  rownames(ExTP1) <- rnames
  # image(ExTP1)
  
  ListaNombres <- strsplit(Result2$tran_P2, 
                           "\\|")
  nTperPath <- sapply(ListaNombres, length)
  # TrNames <- unique(unlist(ListaNombres))
  # # TrNames includes all the possible
  # transcripts.  In transcritnames is
  # sotored the names of all the Isoforms
  TrInPaths <- unlist(ListaNombres)
  i <- rep(seq_len(length(ListaNombres)), 
           nTperPath)
  j <- match(TrInPaths, transcritnames)
  ExTP2 <- sparseMatrix(i, j, x = 1, dims = c(max(i), 
                                              length(transcritnames)))
  rownames(ExTP2) <- rnames
  # image(ExTP2)
  
  ListaNombres <- strsplit(Result2$tran_Ref, 
                           "\\|")
  nTperPath <- sapply(ListaNombres, length)
  # TrNames <- unique(unlist(ListaNombres))
  # # TrNames includes all the possible
  # transcripts.  In transcritnames is
  # sotored the names of all the Isoforms
  TrInPaths <- unlist(ListaNombres)
  i <- rep(seq_len(length(ListaNombres)), 
           nTperPath)
  j <- match(TrInPaths, transcritnames)
  ExTPRef <- sparseMatrix(i, j, x = 1, 
                          dims = c(max(i), length(transcritnames)))
  rownames(ExTPRef) <- rnames
  # image(ExTPRef)
  
  # identical(ExTP1+ExTP2,ExTPRef)
  
  PathsxTranscript <- list(ExTP1 = ExTP1, 
                           ExTP2 = ExTP2, ExTPRef = ExTPRef, 
                           transcritnames = transcritnames,
                           SG_List=SG_List)
  
  return(PathsxTranscript)
  
  cat("\n\n\t******FINISHED******\n")
  
}






