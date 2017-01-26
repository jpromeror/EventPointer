# Classification of events into canonical categories
#
#' Internal function called by EvenPointer
#'
#' Function to classify the alternative splicing events into the canonical categories: cassette exon, alternative 3', alternative 5',
#' alternative start site, alternative termination site, retained intron, mutually exclusive exons or complex events
#'
#' @keywords internal
#'
#' @param Events Detected alternative splicing events for a particular gene
#' @param PSR_Gene Exon probes mapped to the gene that is being analyzes
#' @param Junc_Gene Junction probes mapped to the gene that is being analyzes
#' @param Gxx Identifier for the gene
#'
#'
#' @return Events list with Type slot
#'

ClassifyEvents<-function(SG,Events)
{
  # Create a vector to store the type of event for all of the events found
  Type<-vector(length=length(Events))

  # Repeat the process as many events there are
  for(XX in 1:length(Events))
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






