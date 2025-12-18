#' Detection of alternative splicing events from BAM alignment files
#'
#' @description
#' Identifies alternative splicing events directly from BAM alignment files using
#' a splicing graph approach. This function builds a splicing graph from read
#' alignments, detects events, and calculates PSI values with bootstrap resampling.
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
#' @param nboot Number of resamples of the quantification of the samples to
#'   perform bootstrap. Default is 20.
#' @param lambda Lambda parameter for PSI calculation. Default is NULL.
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
#'   will be suppressed. Default is FALSE
#'
#' @return Invisibly returns NULL. Results are saved to files in PathSGResult.
#'
#' @export
#' @importFrom SGSeq convertToTxFeatures convertToSGFeatures annotate TxFeatures mergeTxFeatures processTerminalExons makeSGFeatureCounts geneID type featureID
#' @importFrom parallel makePSOCKcluster stopCluster makeCluster parLapplyLB clusterApplyLB
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom Rsamtools BamFile scanBamFlag ScanBamParam scanBam bamFlagTest idxstatsBam
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments grglist njunc junctions coverage
#' @importFrom GenomicRanges seqinfo GRanges strand findOverlaps
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges IRanges IRangesList flank width start end reduce ranges runLength RleList IntegerList slice
#' @importFrom GenomeInfoDb seqlevelsStyle seqnames 'seqlevelsStyle<-' seqlevels seqlengths
#' @importFrom S4Vectors Rle queryHits subjectHits Hits DataFrame mcols 'mcols<-' findMatches elementNROWS
#' @importFrom matrixStats rowMeans2
#' @importFrom abind abind
#' @importFrom methods as is
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils read.table
#' @importFrom aroma.light medianPolish
#' @importFrom stringr str_sub
#' @importFrom nnls nnls
#'
#' @examples
#' \dontrun{
#' PathSamplesAbundance <- system.file("extdata/bams", package = "EventPointer")
#' PathTranscriptomeGTF <- list.files(PathSamplesAbundance, "*.gtf", full.names = TRUE)
#'
#' PathSGResult <- tempdir()
#'
#' EventsDetection_BAM(
#'   PathSamplesAbundance,
#'   PathTranscriptomeGTF,
#'   region = 16,
#'   min_junction_count = 2,
#'   max_complexity = 30,
#'   min_n_sample = NULL,
#'   min_anchor = 6,
#'   nboot = 20,
#'   lambda = NULL,
#'   cores = 1,
#'   PathSGResult = PathSGResult
#' )
#' }
EventsDetection_BAM <- function(PathSamplesAbundance,
                                PathTranscriptomeGTF = NULL,
                                region = NULL,
                                min_junction_count = 2,
                                max_complexity = 30,
                                min_n_sample = NULL,
                                min_anchor = 6,
                                nboot = 20,
                                lambda = NULL,
                                cores = 1,
                                PathSGResult = ".",
                                verbose = FALSE) {
  
  
  #first step: creating splicing graph ##########
  # cat("Getting BAM general information.../n")
  # 
  # Bam_Info <- getBamInfo(PathSamplesAbundance,region = region, cores = cores)
  # 
  # cat("Obtaining Reference Transcriptome...")
  # 
  # 
  # if (verbose) {
  #   TxDb <- makeTxDbFromGFF(file = PathTranscriptomeGTF,
  #                           format = "gtf", dataSource = "External Transcriptome")
  # } else {
  #   TxDb <- suppressWarnings(makeTxDbFromGFF(file = PathTranscriptomeGTF,
  #                                            format = "gtf", dataSource = "External Transcriptome"))
  # }
  # 
  # TxF_Ref <- convertToTxFeatures(TxDb)
  # 
  # cat("Done")
  # 
  # cat("\n Predicting Features from BAMs...")
  # if (is.null(region)) {
  #   si <- seqinfo(BamFile(Bam_Info$file_bam[1]))
  #   sl <- rep(seqlevels(si),2)
  #   st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
  #   which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
  #   cat("\n Using this regions:\n")
  #   print(data.frame(which@seqnames,which@ranges))
  # }else{
  #   si <- seqinfo(BamFile(Bam_Info$file_bam[1]))
  #   sl <- rep(seqlevels(si),2)
  #   st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
  #   which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
  #   which <- which[region]
  #   cat("\n Using this regions:\n")
  #   print(data.frame(which@seqnames,which@ranges))
  #   
  # }
  # 
  # if(is.null(min_n_sample)){
  #   min_n_sample<-min(c(length(Bam_Info$file_bam),3))
  # }
  # 
  # if (verbose) {
  #   seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(which)
  # } else {
  #   suppressWarnings(seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(which))
  # }
  # 
  # TxF_RefLevels <- TxF_Ref[which(as.vector(seqnames(TxF_Ref)) %in% as.vector(seqnames(which)))]
  # 
  # cat("\n Creating the splicing graph from the alignment files. This will take some time...")
  # if (!is.null(which)) {
  #   TxF_mod <- predictTxFeatures(Bam_Info, cores = cores, which = which, 
  #                                min_junction_count = min_junction_count,
  #                                min_n_sample = min_n_sample, 
  #                                max_complexity = max_complexity,
  #                                min_anchor = min_anchor)  
  # }else{
  #   TxF_mod <- predictTxFeatures(Bam_Info, cores = cores, 
  #                                min_junction_count = min_junction_count,
  #                                min_n_sample = min_n_sample,
  #                                max_complexity = max_complexity,
  #                                min_anchor = min_anchor)  
  # }
  # closeAllConnections()
  # if (verbose) {
  #   seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(TxF_mod)
  # } else {
  #   suppressWarnings(seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(TxF_mod))
  # }
  # 
  # features <- convertToSGFeatures(TxF_mod)
  # features <- annotate(features, TxF_RefLevels)
  # valid_rows <- which(lengths(features@geneName)>0)
  # features <- features[valid_rows]
  # cat("\n Assigning read counts to the paths in the splicing graph...")
  # SgFC <- getSGFeatureCounts(Bam_Info, features,
  #                            min_anchor = min_anchor,
  #                            cores = cores)
  # closeAllConnections()
  # cat("\n Detect event from splicing graph...")
  # SgFC <- annotate(SgFC, TxF_RefLevels)
  # 
  # save(SgFC,file=paste0(PathSGResult,"/SgFC.RData"))
  
  SgFC <- PrepareBam_EP(PathSamplesAbundance = PathSamplesAbundance,
                          PathTranscriptomeGTF = PathTranscriptomeGTF,
                          region = region,
                          min_junction_count = min_junction_count,
                          max_complexity = max_complexity,
                          min_n_sample = min_n_sample,
                          min_anchor = min_anchor,
                          cores = cores,
                          PathSGResult = PathSGResult,
                          verbose = verbose)
  
  
  #end of first step.
  #second step: Detecting splicing events ##########
  EventsDetection_EPBAM<-EventDetection(SgFC, cores=cores)
  
  cat("\n Checking that detected events are previously annotated...")

  save(EventsDetection_EPBAM, file=paste0(PathSGResult,"/EventsDetection_EPBAM.RData"))
  
  event_list <- list()
  index <- 1
  
  for (geneEvent in EventsDetection_EPBAM) {
    for (event in geneEvent) {
      event_list[[index]] <- event$Info
      index <- index + 1
    }
  }
  
  totalEventTable <- do.call(rbind, event_list)
  
  write.csv(totalEventTable, file = paste0(PathSGResult, "/TotalEventsFound.csv"), row.names = FALSE)
  #end of step 2.
  
  
  #third step: compute PSI values ##########
  PSI_boots <- getPSI_RNASeq_boot(EventsDetection_EPBAM, lambda = lambda, cores = cores, nboot = nboot)
  save(PSI_boots, file=paste0(PathSGResult,"/PSI_boots.RData"))
  
  
  cat("\n DONE!")
  closeAllConnections()

}

