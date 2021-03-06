% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrepareBam_EP.R
\name{PrepareBam_EP}
\alias{PrepareBam_EP}
\title{Bam files preparation for EventPointer}
\usage{
PrepareBam_EP(
  Samples,
  SamplePath,
  Ref_Transc = "Ensembl",
  fileTransc = NULL,
  cores = 1,
  Alpha = 2
)
}
\arguments{
\item{Samples}{Name of the .bam files to be analyzed (Sample1.bam,Sample2.bam,...,etc).}

\item{SamplePath}{Path where the bam files are stored.}

\item{Ref_Transc}{Reference transcriptome used to name the genes found in bam files. Options are: Ensembl, UCSC or GTF.}

\item{fileTransc}{Path to the GTF reference transcriptome ff Ref_Transc is GTF.}

\item{cores}{Number of cores used for parallel processing.}

\item{Alpha}{Internal SGSeq parameter to include or exclude regions}
}
\value{
SGFeaturesCounts object. It contains a GRanges object with the corresponding elements to build
the different splicing graphs found and the counts related to each of the elements.
}
\description{
Prepares the information contained in .bam files to be analyzed by EventPointer
}
\examples{
\dontrun{
 # Obtain the samples and directory for .bam files

   BamInfo<-si
   Samples<-BamInfo[,2]
   PathToSamples <- system.file('extdata/bams', package = 'SGSeq')
   PathToGTF<-paste(system.file('extdata',package='EventPointer'),'/FBXO31.gtf',sep='')

  # Run PrepareBam function
   SG_RNASeq<-PrepareBam_EP(Samples=Samples,
                            SamplePath=PathToSamples,
                            Ref_Transc='GTF',
                            fileTransc=PathToGTF,
                            cores=1)
}
}
