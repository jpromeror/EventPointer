#' Bam files preparation for EventPointer
#'
#' Prepares the information contained in .bam files to be analyzed by EventPointer
#'
#' @param Samples Name of the .bam files to be analyzed (Sample1.bam,Sample2.bam,...,etc).
#' @param SamplePath Path where the bam files are stored.
#' @param Ref_Transc Reference transcriptome used to name the genes found in bam files. Options are: Ensembl, UCSC or GTF.
#' @param fileTransc Path to the GTF reference transcriptome ff Ref_Transc is GTF.
#' @param cores Number of cores used for parallel processing.
#' @param Alpha Internal SGSeq parameter to include or exclude regions
#'
#' @return SGFeaturesCounts object. It contains a GRanges object with the corresponding elements to build
#' the different splicing graphs found and the counts related to each of the elements.
#'
#' @examples
#' \dontrun{
#'  # Obtain the samples and directory for .bam files
#'
#'    BamInfo<-si
#'    Samples<-BamInfo[,2]
#'    PathToSamples <- system.file("extdata/bams", package = "SGSeq")
#'    PathToGTF<-paste(system.file("extdata",package="EventPointer"),"/FBXO31.gtf",sep="")
#'
#'   # Run PrepareBam function
#'    SG_RNASeq<-PrepareBam_EP(Samples=Samples,
#'                             SamplePath=PathToSamples,
#'                             Ref_Transc="GTF",
#'                             fileTransc=PathToGTF,
#'                             cores=1)
#'}
#' @export


PrepareBam_EP<-function(Samples,SamplePath,Ref_Transc="Ensembl",fileTransc=NULL,cores=1,Alpha=2)
{
  #Event Pointer for RNASeq Data
  cat("Preparing BAM files for EventPointer...")

  # Create DataFrame (required by SGSeq) with two columns:
  # 1) Sample Name 2) Path to the .bam file

  Location <- paste(SamplePath,"/",Samples,sep="")
  Bams <- data.frame(sample_name=Samples,file_bam=Location,stringsAsFactors = FALSE)

  cat("\n Obtaining Bam Information")
  cat("\n")

  # Get additional information from each bam with getBamInfo
  Bam_Info <- cbind(Bams,getBamInfo(Bams,yieldSize=NULL,cores=cores))

  # A "Reference" transcriptome is needed to, later, annotate
  # the splicing events with corresponding genes.
  # The user should provide either ensembl,UCSC,or GTF File

  cat("Done")

  cat("\n Obtaining Reference Transcriptome...")
  
  stopifnot(input=="Ensembl" | input=="UCSC" | input=="GTF")

  if(Ref_Transc =="Ensembl")
  {
    TxDb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="grch37.ensembl.org")

  }else if(Ref_Transc=="UCSC")
  {
    TxDb <- makeTxDbFromUCSC(genome="hg19", tablename="knownGene")

  }else if(Ref_Transc=="GTF")
  {
    
    stopifnot(!is.null(inputFile))
    
    TxDb <- makeTxDbFromGFF(file = fileTransc, format = "gtf", dataSource = "External Transcriptome")

  }else{

    stop("Unknown Reference Transcriptome")
  }


  # Steps for the Reference Transcriptome

  # Convert the TxDb to Features (GenomicFeatures)
  TxF_Ref <- convertToTxFeatures(TxDb)

  # Convert from TxFeatures to Splicing Graph (Genomic Features)
  SgF_Ref <- convertToSGFeatures(TxF_Ref)

  # Steps for bam files

  cat("Done")

  cat("\n Predicting Features from BAMs...")

  # Predict TxFeatures from the input bam files
  TxF <- predictTxFeatures(Bam_Info,cores=cores,alpha=Alpha)

  # Convert predicted Features to Splicing Graph
  SgF <- convertToSGFeatures(TxF)

  # Get the reads for each subexon and junction
  SgFC <- getSGFeatureCounts(Bam_Info, SgF,cores=cores)

  # Relate the SG Features with the Reference Transcriptome
  seqlevelsStyle(TxF_Ref)<-seqlevelsStyle(SgFC)
  SgFC <- annotate(SgFC, TxF_Ref)
  # SgF<-GenomicRanges:::rowRanges(SgFC)

  # Result<-list(SgF=SgF,SgFC=SgFC)

  return(SgFC)

}

