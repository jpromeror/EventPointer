# Arrange junction array probes
#
#' Internal function called by EvenPointer
#'
#' Function to arrange the .txt files of exon and junction probes in the format required
#' by EventPointer.
#'
#' @keywords internal
#'
#' @param Probes Probes.txt file passed to EventPointer
#' @param Class Specify either PSR for Exon probe or Junction for junction probes
#'
#'
#' @return Probes arranged in required format
#'

PrepareProbes<-function(Probes,Class)
{

if(Class =="PSR")
{

  Probes<-Probes[,c(1:4,8:11,7)]
  colnames(Probes)<-c("Probe ID","X Coord","Y Coord","Gene","Chr","Start","Stop","Strand","Probe Sequence")


}else if(Class=="Junction")
  {

  # There are some probes that have more than 2 alignments (3,4,5),
  # for now we will discard those probes. We should ask Affy.

  Probes<-Probes[,c(1:4,9:11,8)]
  Probes[,5]<-paste("chr",Probes[,5],sep="")
  ix<-str_count(Probes[,6],",")
  ix<-which(ix==1)
  Probes<-Probes[ix,]
  ProbSS<-matrix(as.numeric(unlist(strsplit(Probes[,6],"[,-]"))),ncol=4,byrow=2)[,2:3]
  Probes<-cbind(Probes[,1:5],ProbSS,Probes[,7:8])
  colnames(Probes)<-c("Probe ID","X Coord","Y Coord","Gene","Chr","Start","Stop","Strand","Probe Sequence")
  # Probes<-Probes[ix,]

  }


return(Probes)

}
