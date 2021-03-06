% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateExSmatrix.R
\name{CreateExSmatrix}
\alias{CreateExSmatrix}
\title{Events X RBPS matrix creation}
\usage{
CreateExSmatrix(
  pathtoeventstable,
  SG_List,
  nt = 400,
  Peaks,
  POSTAR,
  EventsRegions = NULL,
  cores = 1
)
}
\arguments{
\item{pathtoeventstable}{Path to eventsFound.txt with the information of all the events}

\item{SG_List}{List with the information of the splicing graph of the genes. Returned by the funciotn EventDetectio_transcriptome}

\item{nt}{Number of nt up and down for the splicing regions of each event}

\item{Peaks}{Table with the peaks}

\item{POSTAR}{Table with peaks of POSTAR}

\item{EventsRegions}{Events regions if calculated prevously. Not need to calculated again.}

\item{cores}{Number of cores if user want to run in parallel.}
}
\value{
The function returns a list with the ExS matrix and with the splicing regions of the events. If the Splicign regions
    is an input of the function then only the ExS matrix will be returned. The ExS matrix is the input for the Splicing Factor
    enrichment analysis.
}
\description{
Generates the Events x RBP matrix for the splicing factor enrichment analysis.
}
