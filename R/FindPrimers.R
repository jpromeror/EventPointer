#' FindPrimers
#' 
#' @description  FindPrimers is the main function of the primers design option.
#' The aim of this function is the design of PCR primers and 
#' TaqMan probes for detection and quantification of alternative splicing.
#' 
#' @description  Depending on the assay we want to carry out the the algorithm will design 
#' the primers for a conventional PCR or the primers and TaqMan 
#' probes if we are performing a TaqMan assay.
#' 
#' @description In the case of a conventional PCR we will be able to detect the alternative
#'  splicing event. Besides, the algorithm gives as an output the length of the PCR
#'  bands that are going to appear. In the case of a TaqMan assay, we will not only
#'  detect but also quantify alternative splicing.
#' 
#' 
#' @param SG Information of the graph of the gene where the selected event belongs.
#'  This information is avaible in the output of EventsGTFfromTranscriptomeGTF function.
#' @param EventNum The "EventNum" variable can be found in the returned .txt file from
#'  the EventsGTFfromTranscriptomeGTF function in the column "EventNumber" or in the output of  EventPointer_RNASeq_TranRef,
#'  the number after the "_" character of the 'Event_ID'.  
#' @param Primer3Path Complete path where primer3_core.exe is placed.
#' @param Dir Complete path where primer3web_v4_0_0_default_settings.txt file and primer3_config directory are stored.
#' @param mygenomesequence genome sequence of reference
#' @param taqman TRUE if you want to get probes and primers for taqman. FALSE if you want to get primers for conventional PCR.
#' @param nProbes Number of probes for Taqman experiments. By default 1.
#' @param nPrimerstwo Number of potential exon locations for primers using
#'  two primers (one forward and one reverse). By default 3.
#' @param ncommonForward Number of potential exon locations for primers using
#'  one primer in forward and two in reverse. By default 3.
#' @param ncommonReverse Number of potential exon locations for primers using 
#'  two primer in forward and one in reverse. By default 3.
#' @param nExons Number of combinations of ways to place primers in exons to 
#'  interrogate an event after sorting. By default 5.
#' @param nPrimers Once the exons are selected, number of primers combination sequences to search within the whole set of potential sequences.
#'  By default 5.
#' @param shortdistpenalty Penalty for short exons following an exponential funciton(A * exp(-dist * shortdistpenalty)).
#'  By defautl 2000.
#' @param maxLength Max length of exons that are between primers and for paths once we have calculated the sequence.
#'  By default 1000.
#' @param minsep Distance from which it is penalized primers for being too close
#'  By default 100.
#' @param wminsep Weigh of the penalization to primers for being too close
#'  By default 200.
#' @param valuethreePenalty penalization for cases that need three primers instead of 2.
#'  By default 1000.
#' @param minexonlength Minimum length that a exon has to have to be able to contain a primer.
#'  By default 25.
#' @param wnpaths Penalty for each existing path
#'  By default 200.
#' @param qualityfilter Results will show as maximum 3 combinations with a punctuation higher than qualityfilter
#'  By default 5000.
#'  
#' @return 
#' The output of the function is a `data.frame` whose columns are:
#' 
#' For1Seq: Sequence of the first forward primer.
#' 
#' For2Seq: Sequence of the second forward primer in case it is needed.
#' 
#' Rev1Seq: Sequence of the first reverse primer.
#' 
#' Rev2Seq: Sequence of the second reverse primer in case it is needed.
#' 
#' For1Exon: Name of the exon of the first forward primer.
#' 
#' For2Exon: Name of the exon of the second forward primer in case it is needed.
#' 
#' Rev1Exon: Name of the exon of the first reverse primer.
#' 
#' Rev2Exon: Name of the exon of the second reverse primer in case it is needed.
#' 
#' FINALvalue: Final punctuation for that combination of exons and sequences. The lower it is this score, the better it is the combination. 
#' 
#' DistPath1: Distances of the bands, in base pairs, that interrogate Path1  when we perform the conventional PCR experiment.
#' 
#' DistPath2: Distances of the bands, in base pairs, that interrogate Path2 `when we perform the conventional PCR experiment.
#' 
#' DistNoPath: Distances of the bands, in base pairs, that they do not interrogate any of the two paths when we perform the conventional PCR experiment.
#' 
#' SeqProbeRef: Sequence of the TaqMan probe placed in the Reference.
#' 
#' SeqProbeP1: Sequence of the TaqMan probe placed in the Path1.
#' 
#' SeqProbeP2: Sequence of the TaqMan probe placed in the Path2.
#' 
#' 
#' @examples
#' \dontrun{
#' 
#' data("EventXtrans")
#' #From the output of EventsGTFfromTranscriptomeGTF we take the splicing graph information
#' SG_list <- EventXtrans$SG_List
#' #SG_list contains the information of the splicing graphs for each gene
#'
#' #Let's supone we want to design primers for the event 1 of the gene ENSG00000254709.7 
#' 
#' #We take the splicing graph information of the required gene
#' SG <- SG_list$ENSG00000254709.7 
#' 
#' #We point the event number
#' EventNum <- 1
#' 
#' #Define rest of variables:
#' Primer3Path <- Sys.which("primer3_core")
#' Dir <- "C:\\PROGRA~2\\primer3\\"
#' 
#' MyPrimers <- FindPrimers(SG = SG,
#'                          EventNum = EventNum,
#'                          Primer3Path = Primer3Path,
#'                          Dir = Dir,
#'                          mygenomesequence = BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
#'                          taqman = 1,
#'                          nProbes=1,
#'                          nPrimerstwo=4,
#'                          ncommonForward=4,
#'                          ncommonReverse=4,
#'                          nExons=10, 
#'                          nPrimers =5,
#'                          maxLength = 1200)  
#' 
#' 
#'}
#' 
#' @export
#' @importFrom igraph graph_from_adjacency_matrix distances as_adjacency_matrix all_simple_paths as_ids
#' @importFrom nnls nnls
#' @importFrom Matrix Diagonal
#' @importFrom BSgenome getSeq 
#' @importFrom stringr str_length str_replace
#' @importFrom S4Vectors isEmpty 
#' @importFrom Biostrings reverseComplement DNAString 

FindPrimers<-function(SG,EventNum,Primer3Path,Dir,mygenomesequence,taqman = NA, nProbes=1, nPrimerstwo=3,ncommonForward=3,
                           ncommonReverse=3,nExons = 5, nPrimers =15, shortdistpenalty=2000,
                           maxLength=1000 , minsep = 100,wminsep = 200,valuethreePenalty= 1000,
                           minexonlength=25,wnpaths = 200,qualityfilter = 5000)
{
  
  if(is.na(taqman)){
    stop("taqman field is empty")
  }
  
  if(!(taqman == TRUE | taqman == FALSE)){
    stop("taqman variable should be equal to TRUE or FALSE")
  }
  
  
  
  if (is.null(SG)){
    stop("SG field is empty")
  }
  
  if (is.null(EventNum)){
    stop("EventNum field is empty or is not > 0")
  }
  
  if (EventNum <=0){
    stop("EventNum field is empty or is not > 0")
  }
  
  
  
  if (is.null(Primer3Path)){
    stop("Primer3Path field is empty")
  }
  
  if (is.null(Dir)){
    stop("Dir field is empty")
  }
  
  if(length(dir(Dir))==0){
    stop("Dir = Complete path where primer3web_v4_0_0_default_settings.txt file and primer3_config directory are stored")
  }
  
  
  if(!file.exists(Primer3Path)){
    stop("Primer3Path field should point to the .exe of Primers3.
         Check in the vignette for more information")
  }
  
  
  
  
  
  if(length(grep("chr",SG$Edges$Chr))==0){
    SG$Edges$Chr <- paste0("chr",SG$Edges$Chr)
  }
  
  
  
  
  randSol <- getRandomFlow(SG$Incidence, ncol = 10)
  Events <- findTriplets(randSol)
  Events <- getEventPaths(Events, SG)
  
  # Misprimers <- vector(mode = "list",length = length(EventNum))
  # names(Misprimers) <- EventNum
  
  # for(ttx in 1:length(EventNum)){
  
  Event<- Events[[EventNum]]
  
  
  # We check if the Event is placed in a reverse strand:
  if (Event$P1$Strand[1]=="-"){
    genreverse <- 1
    print("The event has been designed with reverse strand.")
  }else{
    genreverse <- 0
  }
  # Get general data for all possible combinations of exons.
  generaldata <- getgeneraldata(SG,Event,shortdistpenalty)
  
  # if(any(sapply(generaldata$exonsPathsandRef,length)==0)){
  #   Misprimers[[ttx]] <- "Not possible to place primers due to the structure of the Event."
  #   next()
  # }
  
  # Search and ranking of suitable combinations of exons.
  FinalExons <- try(getFinalExons(generaldata, 
                                  maxLength,
                                  nPrimerstwo,
                                  ncommonForward,
                                  ncommonReverse,
                                  nExons, 
                                  minsep,wminsep,
                                  valuethreePenalty,
                                  minexonlength),silent = TRUE)
  # if(classsss(FinalExons)=="try-error"){
  if(is(FinalExons,"try-error")){
    # Misprimers[[ttx]] <- "Not possible to place primers due to the structure of the Event."
    # next()
    FinalInfo <- "Not possible to place primers due to the structure of the Event."
    return(FinalInfo)
  }
  # Get sequences for best suitables combinations.
  FinalSeq <- try(PrimerSequenceGeneral(taqman, FinalExons,
                                        generaldata,SG,Dir, 
                                        nPrimers,
                                        Primer3Path=Primer3Path,maxLength,minsep,wminsep,
                                        valuethreePenalty,wnpaths,qualityfilter,mygenomesequence),silent = TRUE)
  
  # if(classsss(FinalSeq)=="try-error"){
  if(is(FinalSeq,"try-error")){
    # Misprimers[[ttx]] <- "Not possible to place primers due to the structure of the Event."
    # next()
    FinalInfo <- "Not possible to place primers due to the structure of the Event."
    return(FinalInfo)
  }
  
  # if (classsss(FinalSeq)!="character"){
  if (!is(FinalSeq,"character")){
    PrimerProbes <- data.frame(cbind(rep(NA,nrow(FinalSeq)),
                                     rep(NA,nrow(FinalSeq)),
                                     rep(NA,nrow(FinalSeq))))[,FALSE]
    # Get TaqMan probes for  if they are needed
    if (taqman==1){
      PrimerProbes <- ProbesSequence(SG,FinalSeq,generaldata,Dir,
                                     Primer3Path=Primer3Path,
                                     nProbes,mygenomesequence)
    }
    # BindInfo
    FinalInfo <- cbind(FinalSeq[,-c(5,6,7,8,13)],PrimerProbes)
    
    # Arrange the output if the gene is situated on the negative strand.
    if (genreverse == 1){
      FinalInfo <- genreverse(FinalInfo,taqman)
    }
    
  }else{
    FinalInfo <-FinalSeq
  }
  
  return(FinalInfo)
  # Misprimers[[ttx]] <- FinalInfo
  # }
  
  # return(Misprimers)
  
}
