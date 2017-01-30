#' EventPointer Internal Functions
#'
#' Internal functions used by EventPointer in the different steps of the algorithm
#'
#' @keywords internal
#' @name InternalFunctions
#' 
#'
NULL

#' @rdname InternalFunctions
annotateEvents<-function(Events,PSR_Gene,Junc_Gene,Gxx)
{
  
  GeneName<-Gxx
  ENSGID<-Gxx
  Chrom<-gsub("chr","",as.vector(Events[[1]]$P1[1,"Chr"]))
  Result<-vector("list")
  Flat<-vector("list")
  mm<-0
  
  for(ii in seq_along(Events))
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

#' @rdname InternalFunctions
AnnotateEvents_RNASeq<-function(Events)
{
  Result<-vector("list",length=length(Events))
  for(ii in seq_along(Events))
  {
    
    
    GeneName<-as.vector(Events[[ii]]$GeneName)
    GeneID<-as.vector(Events[[ii]]$Gene)
    EventNumber<-ii
    EventID<-paste(GeneID,"_",EventNumber,sep="")
    EventType<-Events[[ii]]$Type
    Chrom<-as.vector(Events[[ii]]$P1[1,"Chr"])
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
    
    NEv<-data.frame(EventID,GeneName,EventNumber,EventType,GPos,Path1,Path2,PathR,stringsAsFactors = FALSE)
    Result[[ii]]<-NEv
    
  }
  
  Result<-do.call(rbind,Result)
  colnames(Result)<-c("EventID","Gene","Event Number","Event Type","Genomic Position","Path 1","Path 2","Path Reference")
  
  return(Result)
  
  
}


#' @rdname InternalFunctions
ClassifyEvents<-function(SG,Events)
{
  # Create a vector to store the type of event for all of the events found
  Type<-vector(length=length(Events))
  
  # Repeat the process as many events there are
  for(XX in seq_along(Events))
  {
    # Keep the components of Path 1 and Path 2
    P1<-Events[[XX]]$P1[,1:2]
    P2<-Events[[XX]]$P2[,1:2]
    Info<-rbind(P1,P2)
    
    # If there is an edge that leaves the Start node, we have an
    # Alternative First Exon
    if(any(Info[,1]=="S"))
    {
      Type[XX]<-"Alternative First Exon"
      next
    }
    
    # If there is an edge that enters the End node, we have an
    # Alternative Last Exon
    
    if(any(Info[,2]=="E"))
    {
      Type[XX]<-"Alternative Last Exon"
      next
    }
    
    # Create a Mini Adjacency Graph using only the elements
    # from  Path 1 and Path 2
    
    # Find the From Nodes and To Nodes in the complete Adjacency Matrix
    ii<-sort(match(Info[,1],rownames(SG$Adjacency)))
    jj<-sort(match(Info[,2],rownames(SG$Adjacency)))
    
    # The new one will have those rows and columsn
    MiniSGA<-SG$Adjacency[ii,jj]
    
    # Get unique values, as the elements might be repeated
    iixx<-unique(rownames(MiniSGA))
    jjxx<-unique(colnames(MiniSGA))
    MiniSGA<-MiniSGA[iixx,jjxx,drop=FALSE]
    
    
    # If MiniSGA has a dimension of 3x3, the event could be:
    # Cassette, Retained Intron, Alt 3' , Alt 5' or Complex
    #  1a --> 1b --> 2a --> 2b --> 3a --> 3b
    #          \------------------/
    
    #                 2a  2b  3a
    #              1b 1   .   1
    #              2a .   1   .
    #              2b .   .   1
    
    
    if(dim(MiniSGA)[1]==3 & dim(MiniSGA)[2]==3)
    {
      Type[XX]<-"Group1"
      next
    }
    
    # The case of Mutually Exclusive Exons
    
    #          /----> 2a --> 2b---\
    #  1a --> 1b                   4a--> 4b
    #          \----> 3a --> 3b---/
    
    #                 2a  2b  3a  3b  4a
    #              1b 1   .   1   .   .
    #              2a .   1   .   .   .
    #              2b .   .   *   .   1
    #              3a .   .   .   1   .
    #              3b .   .   .   .   1
    #
    #              * Element Must be Zero
    
    if(dim(MiniSGA)[1]==5 & dim(MiniSGA)[2]==5)
    {
      # Get Row and Column Names from MiniSGA and join them in a vector
      Names<-c(rownames(MiniSGA),colnames(MiniSGA))
      
      # Substitute remove .a./b from the Names Vector and keep unique values
      Names<-unique(gsub("[.ab]","",Names))
      
      # Transform from character to numeric
      Names<-as.numeric(Names)
      
      # The position 3,3 from the MiniSGA must be 0 and the elements must be
      # next to the other to be a mutually exclusive case.
      if(MiniSGA[3,3]==0 & all(diff(Names)==1))
      {
        Type[XX]<-"Mutually Exclusive Exons"
        next
      }
    }
    
    Type[XX]<-"Complex Event"
    
    
  }
  
  ## Certain conditions need to be obtained to classify all the events
  ## that belong to Group 1 (Cassettes,RetIntrons,Alt3,Alt5..)
  
  Gs1<-which(Type=="Group1")
  
  for(jj in Gs1)
  {
    # Get the nodes from Path 1 and order them in increasing order,
    # as the elements are character (11 > 21), we must sort them
    # numerically and then add again the .a and .b values
    P1<-c(Events[[jj]]$P1[,"From"],Events[[jj]]$P1[,"To"])
    P1<-unique(gsub("[.ab]","",P1))
    P1<-matrix(sort(c(paste(P1,".a",sep=""),paste(P1,".b",sep=""))),ncol=2,byrow=TRUE)
    
    # match the P1 matrix rows in the SG$Edges matrix to get te information from the edges
    iix<-which(outer(SG$Edges[,"From"],P1[,1],"==") & outer(SG$Edges[,"To"],P1[,2],"=="),arr.ind=TRUE)[,1]
    
    # Keep the edges
    P1<-SG$Edges[iix,]
    P1<-P1[rownames(P1)[order(as.numeric(rownames(P1)))],]
    
    # Conditions to indicate type of event
    Cond<-(as.numeric(P1[c(2:3),"Start"])-1)-(as.numeric(P1[c(1:2),"End"]))
    Cond2<-diff(as.numeric(rownames(P1)))
    
    if(any(Cond2>1) | any(is.na(Cond)))
    {
      Type[jj]<-"Complex Event"
      next
      
    }
    
    if(Cond[1]>0 & Cond[2]>0)
    {
      Type[jj]<-"Cassette Exon"
      next
    }
    
    if(Cond[1]==0 & Cond[2]==0)
    {
      Type[jj]<-"Retained Intron"
      next
    }
    
    if(Cond[1]>0 & Cond[2]==0)
    {
      Type[jj]<-"Alternative 5' Splice Site"
      next
    }
    
    if(Cond[1]==0 & Cond[2]>0)
    {
      Type[jj]<-"Alternative 3' Splice Site"
      next
    }
    
  }
  
  if(SG$Edges[1,"Strand"]=="-")
  {
    Type[which(Type=="Alternative 5' Splice Site")]<-"Alternative 3' Splice Site"
    Type[which(Type=="Alternative 3' Splice Site")]<-"Alternative 5' Splice Site"
    Type[which(Type=="Alternative First Exon")]<-"Alternative Last Exon"
    Type[which(Type=="Alternative Last Exon")]<-"Alternative First Exon"
    
  }
  
  
  for(ii in 1:length(Events))
  {
    Events[[ii]]$Type<-Type[ii]
    
  }
  
  
  return(Events)
}


#' @rdname InternalFunctions

##########################################################
# Given the signals of the three paths,
# estimate the concentrations of each of the isoforms
##########################################################

# lambda parameter included to regularize the affinities

estimateAbsoluteConc <- function(Signal1, Signal2, SignalR, lambda ) {
  # require(nnls)
  Signal1 <- as.numeric(Signal1)
  Signal2 <- as.numeric(Signal2)
  SignalR <- as.numeric(SignalR)
  
  A <- cbind(Signal1, Signal2) 
  b <- SignalR
  Salida <- nnls(A,b) # non negative least squares
  if (lambda == 0) {
    Salida <- Salida$x
    u <- Salida[1]
    v <- Salida[2]
    w <- 0
    offset <- w / (1-u-v) # some times the offset is way too large (1-u-v = 0)
    T1est <- Signal1 * u
    T2est <- Signal2 * v
    Relerror <- as.numeric(crossprod((A[,1:2])%*%c(u,v)-b)/crossprod(b))
    return(list(T1est = T1est, T2est=T2est, offset = offset, Relerror = Relerror))
  }
  # Add a new equation to make the values of u and v close to each other
  if (is.null(lambda)) lambda <- 1
  penalty <- sum(Salida$residuals)/abs(diff(Salida$x))*lambda*abs(log(Salida$x[1]+1e-3)-log(Salida$x[2]+1e-3))/sqrt(3)
  A <- rbind(A,c(penalty,-penalty),c(penalty,0),c(0,penalty))
  b <- c(SignalR,0,penalty,penalty)
  Salida <- nnls(A,b)$x # non negative least squares
  u <- Salida[1]
  v <- Salida[2]
  w <- 0 # No offset used in this model
  offset <- w / (1-u-v) # some times the offset is way too large (1-u-v = 0)
  T1est <- Signal1 * u
  T2est <- Signal2 * v
  Relerror <- as.numeric(crossprod(cbind(Signal1, Signal2)%*%c(u,v)-SignalR)/crossprod(SignalR))
  # if(Relerror==0){browser()}
  return(list(T1est = T1est, T2est=T2est, offset = offset, Relerror = Relerror))
}

#' @rdname InternalFunctions
findTriplets<-function(randSol,tol=1e-8)
{
  
  # Compute the Distance Matrix from the matrix of fluxes
  X<-as.matrix(dist(randSol))
  
  # Which distances are smaller than the tolerance (To ensure the flux is the same always)
  Inc<-(X<tol)
  
  # Create a graph from the adjacency matrix and find the connected components
  g<-graph_from_adjacency_matrix(Inc)
  Groups<-clusters(g)
  
  EdgG_Flux<-randSol[match(1:Groups$no,Groups$membership),]
  
  # All possible combination of two elements from the graph to create
  # all the posible sums to find the triplets of events
  Index <-combn(nrow(EdgG_Flux),2)
  flowsum <- EdgG_Flux[Index[1,],] + EdgG_Flux[Index[2,],]  # All the possible sums
  
  # Calculate the distance between all the possible sums of every element of the graph
  # and the flow matrix. The Events will be those in which the distance is smaller than
  # the tolerance (almost equal to 0). The tolerance is used to avoid rounding problems
  DistanceMat<-pdist2(flowsum,EdgG_Flux)
  
  x_Ref<-which(DistanceMat<tol, arr.ind = TRUE)
  
  # Obtain Paths
  P1<-Index[1,x_Ref[,1]]
  P2<-Index[2,x_Ref[,1]]
  Ref<-x_Ref[,2]
  
  GG<-Groups$membership
  
  triplets<-cbind(P1,P2,Ref)
  
  return(list(groups=GG,triplets=triplets))
  
  
}

#' @rdname InternalFunctions
GetCounts <- function(Events,sg_txiki, type = "counts") {
  readsC <- counts(sg_txiki)
  readsF <- FPKM(sg_txiki)
  countsEvents <- lapply(Events, getPathCounts,readsC)
  countsEvents<-lapply(countsEvents, getPathFPKMs,readsF)
  return(countsEvents)
}

getPathCounts <- function(x, readsC, widthinit) {
  reads <- rbind(colSums(readsC[x$P1$featureID,,drop = FALSE]),
                 colSums(readsC[x$P2$featureID,,drop = FALSE]),
                 colSums(readsC[x$Ref$featureID,,drop = FALSE]))
  
  rownames(reads) <- c("P1","P2","Ref")
  x$Counts<-reads
  return(x)
}

getPathFPKMs <- function(x, readsC, widthinit) {
  reads <- rbind(colSums(readsC[x$P1$featureID,,drop = FALSE]),
                 colSums(readsC[x$P2$featureID,,drop = FALSE]),
                 colSums(readsC[x$Ref$featureID,,drop = FALSE]))
  
  rownames(reads) <- c("P1","P2","Ref")
  x$FPKM<-reads
  return(x)
}

#' @rdname InternalFunctions
getEventPaths<-function(Events,SG)
{
  Exx<-vector("list",length=nrow(Events$triplets))
  
  for(ii in seq_len(nrow(Events$triplets)))
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

#' @rdname InternalFunctions
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
  
  for(ii in seq_len(nrow(Path1)))
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
  
  for(ii in seq_len(nrow(Path2)))
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
  
  
  for(ii in seq_along(Ref.Group))
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
  
  for(ii in seq_len(nrow(Reference)))
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

#' @rdname InternalFunctions
getPSI <- function(ExFit) {
  # Create matrix to fill with PSI values (1 per event and sample)
  PSI <- matrix(0, nrow = nrow(ExFit)/3, ncol = ncol(ExFit)-5)
  colnames(PSI)  <- colnames(ExFit[6:ncol(ExFit)])
  rownames(PSI)  <- ExFit[seq(1,nrow(ExFit),by = 3),1]
  
  NCols<-ncol(ExFit)
  # Perform the operations for every detectable alternative splicing event
  for (n in 1:(nrow(ExFit)/3)) {
    
    # Get expression signal from path 1
    Signal1 <- ExFit[1+3*(n-1),6:NCols]
    
    # Get expression signal from path 2
    Signal2 <- ExFit[2+3*(n-1),6:NCols]
    
    # Get expression signal from Reference
    SignalR <- ExFit[3+3*(n-1),6:NCols]
    
    # Function to estimate concentrations from the interrogated isoforms
    Output <- estimateAbsoluteConc(Signal1, Signal2, SignalR, lambda = 1)
    
    # Compute the actual PSI value (T1/T1+T2)
    psi <- Output$T1est / (Output$T1est + Output$T2est)
    PSI[n,] <- psi
  }
  return(PSI)
}

#' @rdname InternalFunctions
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

#' @rdname InternalFunctions
getRandomFlow <- function(Incidence, ncol = 1)
{
  # With the incidence matrix, it is possible to get its null-space and generate an
  # arbitrary flow on it. Using the flow it is possible to get the triplets of events.
  
  # The seed is set to ensure the order of events remains the same
  set.seed("0xABBA")
  
  # Solve the Null Space for the Incidence Matrix
  solh <- Null(t(Incidence))
  
  # COndition to ensure that everything that exits the Start node (-1),
  # exits at the End Node (1)
  solp <- ginv(Incidence) %*% c(-1,rep(0,nrow(Incidence)-2),1)
  
  # Matrix of fluxes, with as many columns as specified by the user
  v <- matrix(runif(ncol(solh)*ncol),ncol=ncol)
  randSol <- as.vector(solp) + solh %*% v
  
  return(randSol)
  
}

#' @rdname InternalFunctions
IHsummarization<-function(Pv1,t1,Pv2,t2, coherence = "Opposite")
{
  
  if (coherence == "Equal") {
    nPv1 <- (Pv1/2)*(t1>0)+(1-Pv1/2)*(t1<=0)
    nPv2 <- (Pv2/2)*(t2>0)+(1-Pv2/2)*(t2<=0)
    Psuma <- nPv1+nPv2
    PIH <- (Psuma^2)/2*(Psuma<1)+(1-(2-Psuma)^2/2)*(Psuma>=1)
    ZIH <- qnorm(PIH)
    PIH_2tail <- PIH*2 * (PIH<0.5) + ((1-PIH)*2) * (PIH>=0.5)
    
    # 1:  Pv1 > 0.5 ; Pv2 > 0.5 ; t1 >0 ; t2 >0
    
    Cambiar <- which(Pv1 < 0.5 & Pv2 < 0.5 & t1 <=0 & t2 <=0)
    nPv1[Cambiar] <- Pv1[Cambiar]/2
    nPv2[Cambiar] <- Pv2[Cambiar]/2
    Psuma[Cambiar] <- nPv1[Cambiar]+nPv2[Cambiar]
    PIH[Cambiar] <- (Psuma[Cambiar]^2)/2
    ZIH[Cambiar] <- -qnorm(PIH[Cambiar])
    PIH_2tail[Cambiar] <- (PIH[Cambiar])*2
    
    return(list(Pvalues=PIH_2tail,Tstats=ZIH))
  }
  else {
    return(IHsummarization(Pv1,t1,Pv2,-t2, coherence = "Equal"))
  }
}

#' @rdname InternalFunctions
pdist2<-function(X,Y)
{
  X1<-rowSums(X*X)
  Y1<-rowSums(Y*Y)
  Z<-outer(X1,Y1,"+")-2*X%*%t(Y)
}


#' @rdname InternalFunctions
PrepareCountData<-function(Result)
{
  CountMatrix<-vector("list",length=length(Result))
  
  for(jj in seq_along(Result))
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