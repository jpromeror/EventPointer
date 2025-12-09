#' Detection of alternative splicing events from transcriptome quantification
#'
#' @description
#' Identifies alternative splicing events from transcript-level abundance data
#' (Kallisto, Salmon, RSEM) using a reference transcriptome GTF file.
#'
#' @param PathSamplesAbundance Path to quantification data of gene isoforms
#'  obtained from the corresponding pseudoaligner.
#' @param typeAbundance Type of abundance files provided for the PSI analysis.
#' The options are "salmon" or "kallisto".
#' @param PathTranscriptomeGTF Path to the reference transcriptome GTF file.
#' @param EventsTranscriptome Optional. Pre-computed events from transcriptome.
#'   If NULL, events will be detected from the GTF file.
#' @param Bootstrap Logical. If TRUE, bootstrap replicates are used for PSI
#'   calculation. Default is TRUE.
#' @param Filter Logical. If TRUE, filters events based on expression levels.
#'   Default is FALSE.
#' @param Qn Numeric. Quantile threshold for filtering (between 0 and 1).
#'   Default is 0.25.
#' @param cores Integer. Number of cores to use for parallel processing.
#'   Default is 1.
#' @param PathEventsGTFResults Path where results will be saved. Default is
#'   current directory ("."). The following files are obtained:
#'      EventsTranscriptome.RData: Detected events data in RData format.
#'      PSI_ST.RData: $\Psi$ per event and sample in RData format.
#'
#' @return A list containing PSI (Percent Spliced In) values and associated
#'   statistics for detected splicing events across samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' PathFiles <- system.file("extdata", package = "EventPointer")
#' PathTranscriptomeGTF <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep = "")
#' PathSamplesAbundance <- paste0(PathFiles, "/output")
#' PathSamplesAbundance <- dir(PathSamplesAbundance, full.names = TRUE)
#' Pathtxt <- tempdir()
#'
#' EventsPSI <- EventsDetection_ST(
#'   PathSamplesAbundance,
#'   typeAbundance = "kallisto",
#'   PathTranscriptomeGTF = PathTranscriptomeGTF,
#'   EventsTranscriptome = NULL,
#'   Bootstrap = TRUE,
#'   Filter = FALSE,
#'   Qn = 0.25,
#'   cores = 1,
#'   PathEventsGTFResults = Pathtxt
#' )
#' }

EventsDetection_ST <- function(PathSamplesAbundance = NULL,
                               typeAbundance = "kallisto",
                               PathTranscriptomeGTF = NULL,
                               EventsTranscriptome = NULL,
                               Bootstrap = TRUE,
                               Filter = FALSE,
                               Qn = 0.25,
                               cores = 1,
                               PathEventsGTFResults = "."){

  #Event detection from transcriptome gtf
  if(!is.null(PathTranscriptomeGTF) & is.null(EventsTranscriptome)){ 
    EventsTranscriptome <- EventDetection_transcriptome(inputFile = PathTranscriptomeGTF,
                                                        Pathtxt=PathEventsGTFResults,
                                                        cores=cores)
    save(EventsTranscriptome,file=paste0(PathEventsGTFResults,"/EventsTranscriptome.RData"))
  }
  
  
  
  if(!is.null(EventsTranscriptome) & !is.null(PathSamplesAbundance)){ 
    data_exp <- getbootstrapdata(PathSamples = PathSamplesAbundance,type = typeAbundance)
    
    rownames(data_exp[[1]]) <- sapply(strsplit(rownames(data_exp[[1]]),"\\|"),function(X) return(X[1]))
    
    PSI <- GetPSI_FromTranRef(PathsxTranscript = EventsTranscriptome,
                              Samples = data_exp,
                              Bootstrap = Bootstrap,
                              Filter = Filter,
                              Qn = Qn)
  }else{
    stop("Events not provided")
  }
  save(PSI,file=paste0(PathEventsGTFResults,"/PSI_ST.RData"))
  return(PSI)
  
  
}

