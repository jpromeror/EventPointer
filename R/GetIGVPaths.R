# IGV paths representation
#
#' Internal function called by EvenPointer
#'
#' Function to define the segments to be displayed in IGV to represent an alternative splicing event
#'
#' @keywords internal
#'
#' @param EventInfo Information related to the event that is being displayed
#' @param SG_Edges data.frame with the edges for the splicing graph
#'
#'
#' @return matrix with information to be passed to igv
#'

GetIGVPaths<-function(EventInfo,SG_Edges)
{

  Gene<-as.vector(EventInfo[1,1])

  Path1<-matrix(unlist(strsplit(as.vector(EventInfo[,"Path.1"]),"[,-]")),ncol=2,byrow=TRUE)
  Path2<-matrix(unlist(strsplit(as.vector(EventInfo[,"Path.2"]),"[,-]")),ncol=2,byrow=TRUE)
  PathR<-matrix(unlist(strsplit(as.vector(EventInfo[,"Path.Reference"]),"[,-]")),ncol=2,byrow=TRUE)

  SG_Edges_Orig<-SG_Edges

  SG_Edges[,1]<-gsub(".[ab]","",SG_Edges[,1])
  SG_Edges[,2]<-gsub(".[ab]","",SG_Edges[,2])

  P1.ix<-apply(Path1,1,function(x){ix<-which(SG_Edges[,1]==x[1] & SG_Edges[,2]==x[2]);return(ix)})
  P2.ix<-apply(Path2,1,function(x){ix<-which(SG_Edges[,1]==x[1] & SG_Edges[,2]==x[2]);return(ix)})
  PR.ix<-apply(PathR,1,function(x){ix<-which(SG_Edges[,1]==x[1] & SG_Edges[,2]==x[2]);return(ix)})

  Path1<-SG_Edges[P1.ix,]
  Path2<-SG_Edges[P2.ix,]
  PathR<-SG_Edges[PR.ix,]


  PlotPath1<-c()

  for(ii in 1:nrow(Path1))
  {
    Type<-as.vector(Path1[ii,"Type"])

    if(Type=="J")
    {

      Chr<-as.vector(rep(Path1[ii,"Chr"],2))
      St<-as.numeric(c(as.vector(Path1[ii,"Start"]),as.vector(Path1[ii,"End"])))
      Ed<-St
      Wd<-as.numeric(rep(0,2))
      Str<-as.vector(rep(Path1[ii,"Strand"],2))
      Gn<-as.vector(rep(Gene,2))
      Trs<-rep("A",2)
      Res<-data.frame(chromosome=Chr,start=St,end=Ed,width=Wd,strand=Str,gene=Gn,transcript=Trs,stringsAsFactors = FALSE)
      PlotPath1<-rbind(PlotPath1,Res)

    }else if(Type=="E")
    {
      Chr<-as.vector(Path1[ii,"Chr"])
      St<-as.numeric(as.vector(Path1[ii,"Start"]))
      Ed<-as.numeric(as.vector(Path1[ii,"End"]))
      Wd<-Ed-St
      Str<-as.vector(Path1[ii,"Strand"])
      Res<-data.frame(chromosome=Chr,start=St,end=Ed,width=Wd,strand=Str,gene=Gene,transcript="A",stringsAsFactors = FALSE)
      PlotPath1<-rbind(PlotPath1,Res)

    }

  }


  PlotPath2<-c()

  for(ii in 1:nrow(Path2))
  {
    Type<-as.vector(Path2[ii,"Type"])

    if(Type=="J")
    {

      Chr<-as.vector(rep(Path2[ii,"Chr"],2))
      St<-as.numeric(c(as.vector(Path2[ii,"Start"]),as.vector(Path2[ii,"End"])))
      Ed<-St
      Wd<-as.numeric(rep(0,2))
      Str<-as.vector(rep(Path2[ii,"Strand"],2))
      Gn<-as.vector(rep(Gene,2))
      Trs<-rep("B",2)
      Res<-data.frame(chromosome=Chr,start=St,end=Ed,width=Wd,strand=Str,gene=Gn,transcript=Trs,stringsAsFactors = FALSE)
      PlotPath2<-rbind(PlotPath2,Res)

    }else if(Type=="E")
    {
      Chr<-as.vector(Path2[ii,"Chr"])
      St<-as.numeric(as.vector(Path2[ii,"Start"]))
      Ed<-as.numeric(as.vector(Path2[ii,"End"]))
      Wd<-Ed-St
      Str<-as.vector(Path2[ii,"Strand"])
      Res<-data.frame(chromosome=Chr,start=St,end=Ed,width=Wd,strand=Str,gene=Gene,transcript="B",stringsAsFactors = FALSE)
      PlotPath2<-rbind(PlotPath2,Res)

    }

  }

  Ref.Group<-connectedComp(ftM2graphNEL(as.matrix(SG_Edges_Orig[rownames(PathR),1:2])))


  for(ii in 1:length(Ref.Group))
  {
    LL<-length(Ref.Group[[ii]])
    Ref.Group[[ii]]<-cbind(Ref.Group[[ii]][1:(LL-1)],Ref.Group[[ii]][2:LL])
    ixx<-row.match(as.data.frame(Ref.Group[[ii]]),SG_Edges_Orig[,1:2])
    RR<-SG_Edges_Orig[ixx,]
    RR[,1]<-gsub(".[ab]","",RR[,1])
    RR[,2]<-gsub(".[ab]","",RR[,2])
    Trs<-rep(paste("Ref",ii,sep=""),nrow(RR))
    RR<-cbind(RR,Trs)
    Ref.Group[[ii]]<-RR

  }

  Reference<-do.call(rbind,Ref.Group)
  PlotReference<-c()

  for(ii in 1:nrow(Reference))
  {

    Type<-as.vector(Reference[ii,"Type"])

    if(Type=="J")
    {

      Chr<-as.vector(rep(Reference[ii,"Chr"],2))
      St<-as.numeric(c(as.vector(Reference[ii,"Start"]),as.vector(Reference[ii,"End"])))
      Ed<-St
      Wd<-as.numeric(rep(0,2))
      Str<-as.vector(rep(Reference[ii,"Strand"],2))
      Gn<-as.vector(rep(Gene,2))
      Trs<-rep(Reference[ii,"Trs"],2)
      Res<-data.frame(chromosome=Chr,start=St,end=Ed,width=Wd,strand=Str,gene=Gn,transcript=Trs,stringsAsFactors = FALSE)
      PlotReference<-rbind(PlotReference,Res)

    }else if(Type=="E")
    {
      Chr<-as.vector(Reference[ii,"Chr"])
      St<-as.numeric(as.vector(Reference[ii,"Start"]))
      Ed<-as.numeric(as.vector(Reference[ii,"End"]))
      Wd<-Ed-St
      Str<-as.vector(Reference[ii,"Strand"])
      Res<-data.frame(chromosome=Chr,start=St,end=Ed,width=Wd,strand=Str,gene=Gene,transcript=Reference[ii,"Trs"],stringsAsFactors = FALSE)
      PlotReference<-rbind(PlotReference,Res)

    }


  }

  Plot<-rbind(PlotPath1,PlotPath2,PlotReference)
  # Plot[,1]<-paste("chr",Plot[,1],sep="")

  return(Plot)


}
