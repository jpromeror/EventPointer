#' Events .gtf from trancriptome .gtf
#'
#'
#' @param inputFile If input is GTF, inputFile should point to the GTF file to be used.
#' @param PathGTF Directory where the output will be saved
#' @param Transcriptome the name of the transcriptome
#' @param Pathtxt Directory to save the .txt of the events founded
#'
#' @return a list containing four elements: three sparce matrices that relate which isoforms build up the paths 
#' (path1,path2 and pathRef) of each event. The fourth element contains the name of the reference annotation: 
#' only appear the name of the transcript.  
#'
#' @examples 
#'    PathFiles<-system.file("extdata",package="EventPointer")
#'    inputFile <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep="")
#'    Transcriptome <- "Gencode24_2genes"
#'    Pathtxt <- tempdir()
#'    PathGTF <- tempdir()
#'    
#'    # Run the function 
#'    
#'    EventXtrans <- EventsGTFfromTrancriptomeGTF(inputFile = inputFile,
#'                                                Transcriptome = Transcriptome,
#'                                                Pathtxt=Pathtxt,PathGTF=PathGTF)
#'
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


EventsGTFfromTrancriptomeGTF <- function(inputFile = NULL,Transcriptome = NULL, Pathtxt=NULL, PathGTF=NULL){
  
  if(is.null(inputFile)){
    stop('not inputFile')
  }
  if(is.null(Transcriptome)){
    stop('not Transcriptome')
  }
  if(is.null(Pathtxt)){
    stop('not Pathtxt')
  }
  if(is.null(PathGTF)){
    stop('not PathGTF')
  }
  
  cat("Creating SG Information...")
  
  TxDb <- makeTxDbFromGFF(file = inputFile, format = "gtf", dataSource = "Custom GTF")
  TranscriptFeatures <- convertToTxFeatures(TxDb)
  #SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
  
  #Nombre de todos los transcritos posibles
  transcritnames <- txName(TranscriptFeatures)
  transcritnames <- unique(unlist(transcritnames))
  # length(transcritnames) must be 199169 for gencode.v24 (the same than the fasta)
  

  
  
  #Metodo de SG_Seq corregido
  
  
  #SplicingGraphFeatures2 <- convertToSGFeatures2(TranscriptFeatures)
  
  genes22 <- unique(unlist(geneName(TranscriptFeatures)))
  numberofgenes<-length(genes22)
  # genes2 <- unique(unlist(geneID(SplicingGraphFeatures2)))
  # 
  # GeneIDs2 <- as.data.frame(SplicingGraphFeatures2)
  # GeneIDs2 <- GeneIDs2[, c("geneID", "geneName")]
  # IdName2 <- apply(GeneIDs2, 1, function(X) {
  #   A <- unlist(X[2])[1]
  #   return(A)
  # })
  # GeneIDs2[, 2] <- IdName2
  # GeneIDs2 <- unique(GeneIDs2)

  
  # Result <- vector("list", length = length(genes22))
  Result2 <- vector("list", length = numberofgenes)
  
  SG_Edges <- vector("list",length = numberofgenes)
  names(SG_Edges) <- genes22
  
  pb <- txtProgressBar(min = 0, max = numberofgenes, style = 3)
  
  # 1:numberofgenes
  
  
  for (jj in 1:numberofgenes) {
    
    #cat(jj,'\n')
    setTxtProgressBar(pb, jj)
    
    Gene <- genes22[jj]
    
    SG_Gene <- TranscriptFeatures[which(any(geneName(TranscriptFeatures)==Gene)),]
    
    
    SG_Gene <- convertToSGFeatures2(SG_Gene)
    
    
    # Get the number of transcripts that appear for the gene that is being analyzed
    
    Txx <- unique(unlist(as.data.frame(SG_Gene)[, "txName"]))
    
    # In order to have alternative splicing a gene must have at least more than one
    # exon and more than one transcript
    
    SG_Gene2 <- SG_Gene[which((start(SG_Gene)-end(SG_Gene))!=0)]
    SG_Gene_SoloE <- SG_Gene2[type(SG_Gene2)=="E"]
    
    
    
    if (nrow(as.data.frame(SG_Gene)) != 1 & length(Txx) != 1 & nrow(as.data.frame(SG_Gene_SoloE))>1) {
      # Function to obtain the information of the SG, the information returned
      # Corresponds to the Edges, Adjacency matrix and Incidence Matrix of the SG that
      # is being analyzed
      
      #SG <- SG_Info(SG_Gene)
      SG <- SG_creation_RNASeq(SG_Gene)
      
      SG_Edges[[jj]]<-SG$Edges
      
      #plot(ftM2graphNEL(as.matrix(SG$Edges[,1:2]), edgemode="directed"))
    
      # Using Ax=b with A the incidence matrix and b=0, we solve to obtain the null
      # space and obtain the flux through each edge if we relate the SG with an
      # hydraulic installation
      
      # The new parameter (ncol) is used to generate 10 fluxes to ensure that the
      # eventdetection is exaclty the same every step.
      
      randSol <- getRandomFlow(SG$Incidence, ncol = 10)
      
      # To find events, we need to have fluxes different from 1 in different edges
      
      if (any(round(randSol) != 1)) {
        # Identify the alternative splicing events using the fluxes obtained previously.
        # The number of fluxes correspond to the number of edges so we can relate events
        # with edges
        
        Events <- findTriplets(randSol)
        
        # There should be at least one event to continue the algorithm
        
        if (nrow(Events$triplets) > 0) {
          
          twopaths <- which(rowSums(Events$triplets != 0)==3)
          
          # Transform events into triplets, related with edges to the corresponding paths.
          # Also the edges are divided into Path 1, Path 2 and Reference
          
          Events <- getEventPaths(Events, SG)
          
          # Condition to ensure that there is at least one event
          
          if (length(Events) > 0) {
            # Get PSRs that map to the gene that is being used
            
            # Classify the events into canonical categories using the adjacency matrix
            
            
            
            Events <- ClassifyEvents(SG, Events,twopaths)
            
            # Obtain the corresponding PSRs/JunctionProbes for each of the paths of every
            # event
            
            GenI <- jj
            Events <- AnnotateEvents_KLL(Events, Gene,GenI)
            
            if(is.null(Events)){
              Result2[[jj]] <- Events
            }else{
              edgetr <- transfromedge(SG,SG_Gene)
              
              transcrits <- sacartranscritos(edgetr,Events)
              
              Events$tran_P1 <- ""
              Events$tran_P2 <- ""
              Events$tran_Ref <- ""
              
              Events$tran_P1 <- as.vector(transcrits$p1)
              Events$tran_P2 <- as.vector(transcrits$p2)
              Events$tran_Ref <- as.vector(transcrits$ref)
              
              Result2[[jj]] <- Events
            }
            
            
          }
          
        }
        
        
      }
      
    }
    
  }
  
  close(pb)
  
  
  
  ############################
  
  #Result metodo de siempre
  #Result2 metodo corregido

  ############################
  # writhe the .txt
  ############################
  cat("\nCreating .txt ...")
  Result2 <- Filter(Negate(is.null), Result2)
  Result2 <- do.call(rbind, Result2)
  
  goodone2 <- comprobaciontranscritos2(Result2)
  
  Result2<-Result2[goodone2,]
  

  colnames(Result2)<-c("GeneName","GeneID","EventNumber","EventType","GPos",       
                        "Path.1","Path.2","Path.Reference","tran_P1","tran_P2","tran_Ref")
  
  write.table(Result2, file = paste(Pathtxt, "/EventsFound_",Transcriptome,".txt", sep = ""), sep = "\t", 
              quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  cat("\n.txt created")
  
  ############################
  # writhe the GTF
  ############################
  cat("\nCreating .GTF ...")
  # Create file to store gtf for patths (events)
  FILE.paths <- paste(PathGTF, "/paths_RNASeq.gtf", sep = "")
  cat(file = FILE.paths, paste("#track name=", shQuote("paths", type = "cmd"), 
                               " gffTags=", shQuote("on", type = "cmd"), sep = ""), "\n")
  
  pb <- txtProgressBar(min = 0, max = nrow(Result2), style = 3)
  
  # 1:nrow(Result2)
  for (jj in 1:nrow(Result2)){
    
    setTxtProgressBar(pb, jj)
    
    Gene <- Result2$GeneName[jj]
    EvtEdg <- SG_Edges[[Gene]]
    
    if (unique(EvtEdg[, "Strand"]) == "") {
      EvtEdg[, "Strand"] <- "-"
    }
    
    #iixx <- which(Result2$GeneName==Gene)
    
    EventPaths <- GetIGVPaths(Result2[jj, ], EvtEdg)
    class(EventPaths[, 2]) <- "integer"
    class(EventPaths[, 3]) <- "integer"
    WriteGTF_RNASeq(PathGTF, Result2[jj, ], EventPaths)
    
  }
  
  close(pb)
  
  cat("\n.txt created")
  
  
  
  ############################
  # sacar matrix path x transcrito
  ############################
  
  cat("\nCreating the sparseMatrix of paths x transcripts...")
  
  rnames <- paste(Result2$GeneName,Result2$EventNumber,sep="_")
  
  
  ListaNombres <- strsplit(Result2$tran_P1, "\\|")
  nTperPath <- sapply(ListaNombres, length)
  #TrNames <- unique(unlist(ListaNombres)) # TrNames includes all the possible transcripts.
  #In transcritnames is sotored the names of all the Isoforms 
  TrInPaths <- unlist(ListaNombres)
  i <- rep(1:length(ListaNombres), nTperPath)
  j <- match(TrInPaths, transcritnames)
  ExTP1 <- sparseMatrix(i, j, x = 1, dims = c(max(i),length(transcritnames)))
  rownames(ExTP1)<-rnames
  #image(ExTP1)  
  
  ListaNombres <- strsplit(Result2$tran_P2, "\\|")
  nTperPath <- sapply(ListaNombres, length)
  #TrNames <- unique(unlist(ListaNombres)) # TrNames includes all the possible transcripts.
  #In transcritnames is sotored the names of all the Isoforms 
  TrInPaths <- unlist(ListaNombres)
  i <- rep(1:length(ListaNombres), nTperPath)
  j <- match(TrInPaths, transcritnames)
  ExTP2 <- sparseMatrix(i, j, x = 1, dims = c(max(i),length(transcritnames)))
  rownames(ExTP2)<-rnames
  #image(ExTP2)  
  
  ListaNombres <- strsplit(Result2$tran_Ref, "\\|")
  nTperPath <- sapply(ListaNombres, length)
  #TrNames <- unique(unlist(ListaNombres)) # TrNames includes all the possible transcripts.
  #In transcritnames is sotored the names of all the Isoforms 
  TrInPaths <- unlist(ListaNombres)
  i <- rep(1:length(ListaNombres), nTperPath)
  j <- match(TrInPaths, transcritnames)
  ExTPRef <- sparseMatrix(i, j, x = 1, dims = c(max(i),length(transcritnames)))
  rownames(ExTPRef)<-rnames
  #image(ExTPRef)  
  
  #identical(ExTP1+ExTP2,ExTPRef)
  
  PathsxTranscript<-list(ExTP1=ExTP1,ExTP2=ExTP2,ExTPRef=ExTPRef,transcritnames=transcritnames)
  
  cat("\n\n\t******FINISHED******\n")
  
  
  return(PathsxTranscript)
  
  
}




