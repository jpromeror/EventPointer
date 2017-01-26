# Prepare EventPointer output fot the user
#
#' Internal function called by EvenPointer
#'
#' Function to create and organize the data.frame that is given to the user as output
#' of EventPointer function
#'
#' @keywords internal
#'
#' @param Result Detected alternative splicing events for a particular gene
#' @param Final Statistical values for the tests made to the events
#'
#' @return data.frame with results from EventPointer analysis
#'


PrepareOutput<-function(Result,Final)
{

  iix<-matrix(as.numeric(unlist(strsplit(Final[,1],"_"))),ncol=2,byrow=TRUE)

  Types<-vector("list",length=nrow(Final))
  Positions<-vector("list",length=nrow(Final))
  GeneList<-vector("list",length=nrow(Final))
  GeneID<-vector("list",length=nrow(Final))

  GeneN<-c()
  Index<-c()
  GeneI<-c()
  for(jj in 1:length(Result))
  {

    GeneN<-c(GeneN,Result[[jj]][[1]]$GeneName)
    GeneI<-c(GeneI,Result[[jj]][[1]]$Gene)
    Index<-c(Index,jj)


  }

  InfoX<-data.frame(GeneName=GeneN,GeneID=as.numeric(GeneI),Index=as.numeric(Index),stringsAsFactors = FALSE)

  # browser()

  for(jj in 1:nrow(Final))
  {
    # print(jj)

    A<-InfoX[match(iix[jj,1],InfoX[,2]),3]
    B<-iix[jj,2]

    EventType<-Result[[A]][[B]]$Type

    Mat<-rbind(Result[[A]][[B]]$P1,Result[[A]][[B]]$P2)
    Mat<-Mat[Mat[,4]!=0,]
    Mat<-Mat[Mat[,5]!=0,]
    Chr<-Result[[A]][[B]]$P1[1,3]

    St<-min(Mat[,4])
    Sp<-max(Mat[,5])

    Position<-paste(Chr,":",St,"-",Sp,sep="")


    Types[[jj]]<-EventType
    Positions[[jj]]<-Position

    GeneList[[jj]]<-Final[jj,1]
    GeneID[[jj]]<-InfoX[A,1]
  }

  Res<-data.frame(Gene=unlist(GeneID),Event_Type=unlist(Types),Position=unlist(Positions),Pvalue=Final[,2],Zvalue=Final[,3],stringsAsFactors = FALSE)
  rownames(Res)<-unlist(GeneList)
  Pval_Order<-order(Res[,"Pvalue"])

  Res<-Res[Pval_Order,]

  return(Res)

}
