# Percent Spliced In estimation (RNASeq)
#
#' Internal function called by EvenPointer
#'
#' Function to calculate the PSI for a particular alternative splicing event
#'
#' @keywords internal
#'
#' @param ExFit data.frame with fpkms of alternative splicing events
#'
#' @return matrix with PSI values along events and samples
#'

getPSI_RNASeq<-function(Result)
{
  CountMatrix<-vector("list",length=length(Result))
  Vec<-c()
  
  for(jj in 1:length(Result))
  {
    # print(jj)
    A<-Result[[jj]]
    
    if(!is.null(A))
    {
      Evs_Counts<-lapply(A,function(X){Res<-X$FPKM;return(Res)})
      names(Evs_Counts)<-1:length(Evs_Counts)
      Ids<-paste(A[[1]]$Gene,"_",names(Evs_Counts),sep="")
      Ids<-rep(Ids,each=3)
      Vec<-c(Vec,Ids)
      Evs_Counts<-do.call(rbind,Evs_Counts)
      
      if(!any(is.na(Evs_Counts)))
      {
        CountMatrix[[jj]]<-Evs_Counts
      }else{
        
        CountMatrix[[jj]]<-NULL
      }
      
      
    }else{
      
      
    }
    
    
  }
  
  
  Ids<-rep(c("_P1","_P2","_Ref"),length(Vec)/3)
  CountMatrix<-do.call(rbind,CountMatrix)
  rownames(CountMatrix)<-paste(Vec,Ids,sep="")
  
  PSI <- matrix(0, nrow = nrow(CountMatrix)/3, ncol = ncol(CountMatrix))
  colnames(PSI)  <- colnames(CountMatrix)
  rownames(PSI)  <- Vec[seq(1,length(Vec),by = 3)]
  
  for (n in 1:(nrow(CountMatrix)/3)) 
  {
    Signal1 <- CountMatrix[1+3*(n-1),]
    Signal2 <- CountMatrix[2+3*(n-1),]
    SignalR <- CountMatrix[3+3*(n-1),]
    Output <- estimateAbsoluteConc(Signal1, Signal2, SignalR, lambda = 1)
    psi <- Output$T1est / (Output$T1est + Output$T2est)
    PSI[n,] <- psi
  }
  
  return(PSI)
  
}