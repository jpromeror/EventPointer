% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EventPointer_RNASeq.R
\name{EventPointer_RNASeq}
\alias{EventPointer_RNASeq}
\title{Statistical analysis of alternative splcing events for RNASeq data}
\usage{
EventPointer_RNASeq(Events, Design, Contrast, Statistic = "LogFC", PSI = FALSE)
}
\arguments{
\item{Events}{Output from EventDetection function}

\item{Design}{The design matrix for the experiment.}

\item{Contrast}{The contrast matrix for the experiment.}

\item{Statistic}{Statistical test to identify differential splicing events, must be one of : LogFC, Dif_LogFC and DRS.}

\item{PSI}{Boolean variable to indicate if PSI should be calculated for every splicing event.}
}
\value{
Data.frame ordered by the splicing p.value . The object contains the different information for each splicing event
such as Gene name, event type, genomic position, p.value, z.value and delta PSI.
}
\description{
Statistical analysis of all the alternative splicing events found in the given bam files.
}
\examples{

   data(AllEvents_RNASeq)
   Dmatrix<-matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1),ncol=2,byrow=FALSE)
   Cmatrix<-t(t(c(0,1)))
   Events <- EventPointer_RNASeq(AllEvents_RNASeq,Dmatrix,Cmatrix,Statistic='LogFC',PSI=TRUE)
}
