#' Bam files preparation for EventPointer
#'
#' Prepares the information contained in .bam files to be analyzed by EventPointer
#'
#' @param PathSamplesAbundance Path to BAM and BAI files or path to folder
#'   containing BAM and BAI files.
#' @param PathTranscriptomeGTF Path to file containing the regions to be analysed
#'   from the BAM files in GTF format.
#' @param region Numerical vector indicating the index of positions (at chromosomal
#'   level) to be analysed from the GTF. Default is NULL so that all regions are
#'   analysed.
#' @param min_junction_count Minimum number of junctions detected in the alignment
#'   to be considered in the splicing graph. Default is 2.
#' @param max_complexity Maximum allowed complexity. If a locus exceeds this
#'   threshold, it is skipped with a warning. Complexity is defined as the maximum
#'   number of unique predicted splice junctions overlapping a given position.
#'   High complexity regions are often due to spurious read alignments and can
#'   slow down processing. Default is 30.
#' @param min_n_sample Minimum number of samples that a junction must have to be
#'   considered. Default is NULL (automatically set to minimum of 3 or total
#'   number of samples).
#' @param min_anchor Minimum number of aligned bases at one end of an exon to
#'   consider a junction. Default is 6.
#' @param cores Number of cores to use for parallel processing. Default is 1.
#' @param PathSGResult Path where results will be saved. The following 4 files
#'   are generated:
#'   \itemize{
#'     \item TotalEventsFound.csv: General data for total events detected in CSV format.
#'     \item EventsDetection_EPBAM.RData: Raw data per event, paths of splicing graph and counts.
#'     \item SgFC.RData: Contains the splicing graph in RData format.
#'     \item PSI_boots.RData: \eqn{\Psi} per event and sample in RData format.
#'   }
#'   Default is current directory (".").
#' @param verbose Logical indicating whether to show warnings and messages.
#'   If FALSE, warnings from internal functions (e.g., makeTxDbFromGRanges)
#'   will be suppressed. Default is TRUE.
#'
#' @return SGFeaturesCounts object. It contains a GRanges object with the corresponding elements to build
#' the different splicing graphs found and the counts related to each of the elements.
#'
#' @export
#' @importFrom txdbmaker makeTxDbFromBiomart makeTxDbFromUCSC makeTxDbFromGFF


PrepareBam_EP <- function(PathSamplesAbundance,
                                PathTranscriptomeGTF = NULL,
                                region = NULL,
                                min_junction_count = 2,
                                max_complexity = 30,
                                min_n_sample = NULL,
                                min_anchor = 6,
                                cores = 1,
                                PathSGResult = ".",
                                verbose = TRUE) {
    #first step: creating splicing graph ##########
  cat("Getting BAM general information.../n")
  
  Bam_Info <- getBamInfo(PathSamplesAbundance,region = region, cores = cores)
  
  cat("Obtaining Reference Transcriptome...")

  
  if (verbose) {
    TxDb <- makeTxDbFromGFF(file = PathTranscriptomeGTF,
                            format = "gtf", dataSource = "External Transcriptome")
  } else {
    TxDb <- suppressWarnings(makeTxDbFromGFF(file = PathTranscriptomeGTF,
                                             format = "gtf", dataSource = "External Transcriptome"))
  }
  
  TxF_Ref <- convertToTxFeatures(TxDb)
  
  cat("Done")
  
  cat("\n Predicting Features from BAMs...")
  if (is.null(region)) {
    si <- seqinfo(BamFile(Bam_Info$file_bam[1]))
    sl <- rep(seqlevels(si),2)
    st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
    which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
    cat("\n Using this regions:\n")
    print(data.frame(which@seqnames,which@ranges))
  }else{
    si <- seqinfo(BamFile(Bam_Info$file_bam[1]))
    sl <- rep(seqlevels(si),2)
    st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
    which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
    which <- which[region]
    cat("\n Using this regions:\n")
    print(data.frame(which@seqnames,which@ranges))
    
  }
  
  if(is.null(min_n_sample)){
    min_n_sample<-min(c(length(Bam_Info$file_bam),3))
  }
  
  if (verbose) {
    seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(which)
  } else {
    suppressWarnings(seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(which))
  }
  
  TxF_RefLevels <- TxF_Ref[which(as.vector(seqnames(TxF_Ref)) %in% as.vector(seqnames(which)))]

  cat("\n Creating the splicing graph from the alignment files. This will take some time...")
  if (!is.null(which)) {
    TxF_mod <- predictTxFeatures(Bam_Info, cores = cores, which = which, 
                                 min_junction_count = min_junction_count,
                                 min_n_sample = min_n_sample, 
                                 max_complexity = max_complexity,
                                 min_anchor = min_anchor)  
  }else{
    TxF_mod <- predictTxFeatures(Bam_Info, cores = cores, 
                                 min_junction_count = min_junction_count,
                                 min_n_sample = min_n_sample,
                                 max_complexity = max_complexity,
                                 min_anchor = min_anchor)  
  }
  closeAllConnections()
  if (verbose) {
    seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(TxF_mod)
  } else {
    suppressWarnings(seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(TxF_mod))
  }
  
  features <- convertToSGFeatures(TxF_mod)
  features <- annotate(features, TxF_RefLevels)
  valid_rows <- which(lengths(features@geneName)>0)
  features <- features[valid_rows]
  cat("\n Assigning read counts to the paths in the splicing graph...")
  SgFC <- getSGFeatureCounts(Bam_Info, features,
                             min_anchor = min_anchor,
                             cores = cores)
  closeAllConnections()
  #cat("\n Detect event from splicing graph...")
  SgFC <- annotate(SgFC, TxF_RefLevels)
  
  save(SgFC,file=paste0(PathSGResult,"/SgFC.RData"))
  #end of first step.
  return(SgFC)
}
