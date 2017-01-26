# annotate Events
#
#' Internal function called by EvenPointer
#'
#' Function to obtain all the information for the alternative splicing events for a particular gene
#'
#' @keywords internal
#'
#' @param Events Detected alternative splicing events for a particular gene
#' @param PSR_Gene Exon probes mapped to the gene that is being analyzes
#' @param Junc_Gene Junction probes mapped to the gene that is being analyzes
#' @param Gxx Identifier for the gene
#'
#'
#' @return Returns a list with the information required for the .flat file and EventsFound.txt
#'

annotateEvents<-function(Events,PSR_Gene,Junc_Gene,Gxx)
{

  GeneName<-Gxx
  ENSGID<-Gxx
  Chrom<-gsub("chr","",as.vector(Events[[1]]$P1[1,"Chr"]))
  Result<-vector("list")
  Flat<-vector("list")
  mm<-0

  for(ii in 1:length(Events))
  {
    Events[[ii]]$Probes_P1<-NULL
    Events[[ii]]$Probes_P2<-NULL
    Events[[ii]]$Probes_Ref<-NULL
    PSR_P1<-c()
    Junc_P1<-c()
    PSR_P2<-c()
    Junc_P2<-c()
    PSR_Ref<-c()
    Junc_Ref<-c()

    ExonsP1<-which(Events[[ii]]$P1[,"Type"]=="E")
    JunctionsP1<-which(Events[[ii]]$P1[,"Type"]=="J")

    if(length(ExonsP1)>0)
    {
      EPP1<-Events[[ii]]$P1[ExonsP1,]
      PSR_P1<-sapply(1:nrow(EPP1),function(x){which(as.numeric(PSR_Gene[,"Start"])>=as.numeric(EPP1[x,"Start"]) & as.numeric(PSR_Gene[,"Stop"])<=as.numeric(EPP1[x,"End"]))})
      PSR_P1<-PSR_Gene[unlist(PSR_P1),1]

    }

    if(length(JunctionsP1)>0)
    {
      JPP1<-Events[[ii]]$P1[JunctionsP1,]
      Junc_P1<-sapply(1:nrow(JPP1),function(x){which(as.numeric(Junc_Gene[,"Start"])==as.numeric(JPP1[x,"Start"]) & as.numeric(Junc_Gene[,"Stop"])==as.numeric(JPP1[x,"End"]))})
      Junc_P1<-Junc_Gene[unlist(Junc_P1),1]
    }


    ExonsP2<-which(Events[[ii]]$P2[,"Type"]=="E")
    JunctionsP2<-which(Events[[ii]]$P2[,"Type"]=="J")

    if(length(ExonsP2)>0)
    {
      EPP2<-Events[[ii]]$P2[ExonsP2,]
      PSR_P2<-sapply(1:nrow(EPP2),function(x){which(as.numeric(PSR_Gene[,"Start"])>=as.numeric(EPP2[x,"Start"]) & as.numeric(PSR_Gene[,"Stop"])<=as.numeric(EPP2[x,"End"]))})
      PSR_P2<-PSR_Gene[unlist(PSR_P2),1]

    }

    if(length(JunctionsP2)>0)
    {
      JPP2<-Events[[ii]]$P2[JunctionsP2,]
      Junc_P2<-sapply(1:nrow(JPP2),function(x){which(as.numeric(Junc_Gene[,"Start"])==as.numeric(JPP2[x,"Start"]) & as.numeric(Junc_Gene[,"Stop"])==as.numeric(JPP2[x,"End"]))})
      Junc_P2<-Junc_Gene[unlist(Junc_P2),1]
    }


    ExonsRef<-which(Events[[ii]]$Ref[,"Type"]=="E")
    JunctionsRef<-which(Events[[ii]]$Ref[,"Type"]=="J")

    if(length(ExonsRef)>0)
    {
      EPRef<-Events[[ii]]$Ref[ExonsRef,]
      PSR_Ref<-sapply(1:nrow(EPRef),function(x){which(as.numeric(PSR_Gene[,"Start"])>=as.numeric(EPRef[x,"Start"]) & as.numeric(PSR_Gene[,"Stop"])<=as.numeric(EPRef[x,"End"]))})
      PSR_Ref<-PSR_Gene[unlist(PSR_Ref),1]

    }

    if(length(JunctionsRef)>0)
    {
      JPRef<-Events[[ii]]$Ref[JunctionsRef,]
      Junc_Ref<-sapply(1:nrow(JPRef),function(x){which(as.numeric(Junc_Gene[,"Start"])==as.numeric(JPRef[x,"Start"]) & as.numeric(Junc_Gene[,"Stop"])==as.numeric(JPRef[x,"End"]))})
      Junc_Ref<-Junc_Gene[unlist(Junc_Ref),1]
    }


    #       PSR_P2<-sapply(1:nrow(Events[[ii]]$P2),function(x){which(as.numeric(PSR_Gene[,"Start"])>=as.numeric(Events[[ii]]$P2[x,"Start"]) & as.numeric(PSR_Gene[,"Stop"])<=as.numeric(Events[[ii]]$P2[x,"End"]))})
    #       PSR_P2<-PSR_Gene[unlist(PSR_P2),6]
    #       Junc_P2<-sapply(1:nrow(Events[[ii]]$P2),function(x){which(as.numeric(Junc_Gene[,"Start"])==as.numeric(Events[[ii]]$P2[x,"Start"]) & as.numeric(Junc_Gene[,"Stop"])==as.numeric(Events[[ii]]$P2[x,"End"]))})
    #       Junc_P2<-Junc_Gene[unlist(Junc_P2),1]
    #
    #       PSR_Ref<-sapply(1:nrow(Events[[ii]]$Ref),function(x){which(as.numeric(PSR_Gene[,"Start"])>=as.numeric(Events[[ii]]$Ref[x,"Start"]) & as.numeric(PSR_Gene[,"Stop"])<=as.numeric(Events[[ii]]$Ref[x,"End"]))})
    #       PSR_Ref<-PSR_Gene[unlist(PSR_Ref),6]
    #       Junc_Ref<-sapply(1:nrow(Events[[ii]]$Ref),function(x){which(as.numeric(Junc_Gene[,"Start"])==as.numeric(Events[[ii]]$Ref[x,"Start"]) & as.numeric(Junc_Gene[,"Stop"])==as.numeric(Events[[ii]]$Ref[x,"End"]))})
    #       Junc_Ref<-Junc_Gene[unlist(Junc_Ref),1]




    Events[[ii]]$Probes_P1<-c(PSR_P1,Junc_P1)
    Events[[ii]]$Probes_P2<-c(PSR_P2,Junc_P2)
    Events[[ii]]$Probes_Ref<-c(PSR_Ref,Junc_Ref)

    if(length(Events[[ii]]$Probes_P1)>0 & length(Events[[ii]]$Probes_P2)>0 & length(Events[[ii]]$Probes_Ref)>0)
    {
      mm<-mm+1

      EventNumber<-ii

      EventType<-Events[[ii]]$Type

      Positions<-rbind(Events[[ii]]$P1,Events[[ii]]$P2)[,4:5]
      Start<-as.numeric(Positions[,1])
      End<-as.numeric(Positions[,2])
      Start<-Start[which(Start!=0)]
      End<-End[which(End!=0)]

      # browser()
      minGPos<-min(Start)
      maxGPos<-max(End)
      GPos<-paste(Chrom,":",minGPos,"-",maxGPos,sep="")

      CP1s<-which(Events[[ii]]$P1[,1]=="S")
      CP1e<-which(Events[[ii]]$P1[,2]=="E")

      if(length(CP1s)>0|length(CP1e)>0)
      {
        CC<-c(CP1s,CP1e)
        Events[[ii]]$P1<-Events[[ii]]$P1[-CC,]
      }

      PS1<-as.numeric(gsub(".[ab]","",Events[[ii]]$P1[,1]))
      PE1<-as.numeric(gsub(".[ab]","",Events[[ii]]$P1[,2]))
      Path1<-as.matrix(cbind(PS1,PE1))
      Path1<-Path1[order(Path1[,1],Path1[,2]),,drop=FALSE]

      CP2s<-which(Events[[ii]]$P2[,1]=="S")
      CP2e<-which(Events[[ii]]$P2[,2]=="E")

      if(length(CP2s)>0|length(CP2e)>0)
      {
        CC<-c(CP2s,CP2e)
        Events[[ii]]$P2<-Events[[ii]]$P2[-CC,]
      }

      PS2<-as.numeric(gsub(".[ab]","",Events[[ii]]$P2[,1]))
      PE2<-as.numeric(gsub(".[ab]","",Events[[ii]]$P2[,2]))
      Path2<-as.matrix(cbind(PS2,PE2))
      Path2<-Path2[order(Path2[,1],Path2[,2]),,drop=FALSE]

      CPRs<-which(Events[[ii]]$Ref[,1]=="S")
      CPRe<-which(Events[[ii]]$Ref[,2]=="E")

      if(length(CPRs)>0|length(CPRe)>0)
      {
        CC<-c(CPRs,CPRe)
        Events[[ii]]$Ref<-Events[[ii]]$Ref[-CC,]
      }

      PSR<-as.numeric(gsub(".[ab]","",Events[[ii]]$Ref[,1]))
      PER<-as.numeric(gsub(".[ab]","",Events[[ii]]$Ref[,2]))
      PathR<-as.matrix(cbind(PSR,PER))
      PathR<-PathR[order(PathR[,1],PathR[,2]),,drop=FALSE]


      Path1<-paste(Path1[,1],"-",Path1[,2],sep="",collapse=",")
      Path2<-paste(Path2[,1],"-",Path2[,2],sep="",collapse=",")
      PathR<-paste(PathR[,1],"-",PathR[,2],sep="",collapse=",")

      ProbesP1<-paste(Events[[ii]]$Probes_P1,collapse=",")
      ProbesP2<-paste(Events[[ii]]$Probes_P2,collapse=",")
      ProbesR<-paste(Events[[ii]]$Probes_Ref,collapse=",")

      NEv<-data.frame(GeneName,ENSGID,EventNumber,EventType,GPos,Path1,Path2,PathR,ProbesP1,ProbesP2,ProbesR,stringsAsFactors = FALSE)
      Result[[mm]]<-NEv


      Tprobes<-rbind(PSR_Gene,Junc_Gene)
      ii.P1<-match(Events[[ii]]$Probes_P1,Tprobes[,1])
      ii.P2<-match(Events[[ii]]$Probes_P2,Tprobes[,1])
      ii.R<-match(Events[[ii]]$Probes_Ref,Tprobes[,1])

      lP1<-length(ii.P1)
      lP2<-length(ii.P2)
      lRef<-length(ii.R)

      xRef<-rep(paste(GeneName,"_",EventNumber,"_Ref",sep=""),lRef)
      xP1<-rep(paste(GeneName,"_",EventNumber,"_P1",sep=""),lP1)
      xP2<-rep(paste(GeneName,"_",EventNumber,"_P2",sep=""),lP2)
      xTot<-rep(paste(GeneName,"_",EventNumber,sep=""),lP1+lP2+lRef)

      AllProbes<-c(Events[[ii]]$Probes_Ref,Events[[ii]]$Probes_P1,Events[[ii]]$Probes_P2)
      flat_gene<-cbind(AllProbes,Tprobes[c(ii.R,ii.P1,ii.P2),c(2,3,9)],c(xRef,xP1,xP2),xTot)

      Flat[[mm]]<-flat_gene



    }

  }

  Result<-do.call(rbind,Result)
  Flat<-do.call(rbind,Flat)
  return(list(Events=Result,Flat=Flat))
}
