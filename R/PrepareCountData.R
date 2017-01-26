# Count/FPKM matrix for all alternative splicing events
#
#' Internal function called by EvenPointer
#'
#' Function to create a matrix with all the counts/fpkms per path and sample detected by EventPointer.
#'
#' @keywords internal
#'
#' @param Result All alternative splicing events detected by EventPointer
#'
#' @return Matrix of counts/fpkms
#'


PrepareCountData<-function(Result)
{
  CountMatrix<-vector("list",length=length(Result))
  
  for(jj in 1:length(Result))
  {
    # print(jj)
    A<-Result[[jj]]
    
    if(!is.null(A))
    {
      # Evs_Counts<-lapply(A,function(X){Res<-X$Counts;return(Res)})
      Evs_Counts<-lapply(A,function(X){Res<-X$FPKM;return(Res)})
      Evs_Counts<-lapply(Evs_Counts,function(X){X<-X[c("Ref","P1","P2"),];return(X)})
      Cols<-ncol(Evs_Counts[[1]])
      Mat<-matrix(unlist(Evs_Counts),ncol=length(A))
      colnames(Mat)<-paste(jj,"_",1:length(A),sep="")
      Samples<-colnames(Evs_Counts[[1]])
      rownames(Mat)<-paste(rep(Samples,each=3),c("_Ref","_P1","_P2"),sep="")
      colnames(Mat)<-paste(A[[1]]$Gene,1:length(A),sep="_")
      
      if(!any(is.na(Mat)))
      {
        CountMatrix[[jj]]<-Mat
      }else{
        
        CountMatrix[[jj]]<-NULL
      }
      
      
    }else{
      
      
    }
    
    
  }
  
  CountMatrix<-do.call(cbind,CountMatrix)
  return(CountMatrix)
  
  
}