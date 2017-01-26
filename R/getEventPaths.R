# Event Paths
#
#' Internal function called by EvenPointer
#'
#' Function to get the corresponding exons and junctions for each of the paths (P1,P2,Ref) for an
#' alternative splicing event detected by EventPointer
#'
#' @keywords internal
#'
#' @param Events Detected alternative splicing events for a particular gene
#' @param SG Splicing graph representation for a particular gene
#'
#'
#' @return Events list with Type slot
#'

getEventPaths<-function(Events,SG)
{
  Exx<-vector("list",length=nrow(Events$triplets))
  
  for(ii in 1:(nrow(Events$triplets)))
  {
    P1<-SG$Edges[which(Events$groups==Events$triplets[ii,1]),]
    P2<-SG$Edges[which(Events$groups==Events$triplets[ii,2]),]
    Ref<-SG$Edges[which(Events$groups==Events$triplets[ii,3]),]
    
    if(nrow(P1)>nrow(P2))
    {
      
      Exx[[ii]]$P1<-P1
      Exx[[ii]]$P2<-P2
      Exx[[ii]]$Ref<-Ref
      
    }else{
      
      Exx[[ii]]$P1<-P2
      Exx[[ii]]$P2<-P1
      Exx[[ii]]$Ref<-Ref
      
    }
    
    
  }
  
  return(Exx)
  
}