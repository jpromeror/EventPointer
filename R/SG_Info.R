# Splicing Graph representation for genes
#
#' Internal function called by EvenPointer
#'
#' Function to create the adjacency matrix, incidence matrix and edge relation to represent
#' the splicing graph of a particular gene.
#'
#' @keywords internal
#'
#' @param SG_Gene GRanges object with elements that build up the splicing graph.
#'
#' @return List with slots adjacency,incidence and Edges
#'

SG_Info<-function(SG_Gene)

  {

     SE_Cond<-FALSE
     TTS_Cond<-FALSE
     TSS_Cond<-FALSE

     # Obtain information of all the elements of the Graph

     Graph<-SGSeq:::exonGraph(SG_Gene,tx_view=FALSE)
     Graph_Nodes <- SGSeq:::nodes(Graph)
     Graph_Edges<- SGSeq:::edges(Graph)
#
#     Graph<-exonGraph(SG_Gene,tx_view=F)
#     Graph_Nodes <- nodes(Graph)
#     Graph_Edges<- edges(Graph)

     Nodes_Pos<-matrix(unlist(strsplit(Graph_Nodes[,2],"[:-]")),ncol=4,byrow=TRUE)
     Edges_Pos<-matrix(unlist(strsplit(Graph_Edges[,3],"[:-]")),ncol=4,byrow=TRUE)

     Graph_Nodes<-cbind(Graph_Nodes[,1],Nodes_Pos,Graph_Nodes[,3:4])
     Graph_Nodes<-as.data.frame(Graph_Nodes,stringsAsFactors=FALSE)
     colnames(Graph_Nodes)<-c("Name","Chr","Start","End","Strand","Type","featureID")
     Nodes<-as.numeric(as.vector(Graph_Nodes[,1]))

     Graph_Edges<-cbind(Graph_Edges[,1:2],Edges_Pos,Graph_Edges[,4:5])
     Graph_Edges<-as.data.frame(Graph_Edges,stringsAsFactors=FALSE)
     colnames(Graph_Edges)<-c("From","To","Chr","Start","End","Strand","Type","featureID")

     if(as.vector(strand(SG_Gene)@values)=="-")
     {
       Graph_Edges<-Graph_Edges[,c(2,1,3:8)]
       colnames(Graph_Edges)<-c("From","To","Chr","Start","End","Strand","Type","featureID")

     }

     # Determine Subexons

     SE.II<-as.numeric(as.vector(Graph_Nodes[2:nrow(Graph_Nodes),3]))-as.numeric(as.vector(Graph_Nodes[1:(nrow(Graph_Nodes)-1),4]))
     SE.II<-which(SE.II==1)
     L_SE.II<-length(SE.II)
     SE_chr<-rep(as.vector(Graph_Nodes[1,2]),L_SE.II)
     SE_St<-Graph_Nodes[SE.II,4]
     SE_Ed<-Graph_Nodes[SE.II+1,3]
     SE_std<-rep(as.vector(Graph_Nodes[1,5]),L_SE.II)
     SE_type<-rep("J",L_SE.II)
     SE_FID<-rep(0,L_SE.II)

     SubExons<-data.frame(SE.II,SE.II+1,SE_chr,SE_St,SE_Ed,SE_std,SE_type,SE_FID,stringsAsFactors = FALSE)
     colnames(SubExons)<-colnames(Graph_Edges)

     if(nrow(SubExons)!=0)
     {
       Graph_Edges<-rbind(Graph_Edges,SubExons)
       SE_Cond<-TRUE

     }

     # Determine Alternative Last

     TTS<-setdiff(Nodes,as.numeric(as.vector(Graph_Edges[,1])))

     if(length(TTS)>0)
     {
       TTS<-paste(TTS,".b",sep="")
       Ln_TTS<-length(TTS)

       EE<-rep("E",length(TTS))
       TTS<-data.frame(TTS,EE,rep(Graph_Edges[1,3],Ln_TTS),rep(0,Ln_TTS),rep(0,Ln_TTS),rep(as.vector(Graph_Edges[1,6]),Ln_TTS),rep("J",Ln_TTS),rep(0,Ln_TTS),stringsAsFactors = FALSE)
       colnames(TTS)<-colnames(Graph_Edges)

       TTS_Cond<-TRUE

     }


     # Determine Alternative First


     TSS<-setdiff(Nodes,as.numeric(as.vector(Graph_Edges[,2])))

     if(length(TSS)>0)
     {
       TSS<-paste(TSS,".a",sep="")
       Ln_TSS<-length(TSS)

       SS<-rep("S",Ln_TSS)
       TSS<-data.frame(SS,TSS,rep(Graph_Edges[1,3],Ln_TSS),rep(0,Ln_TSS),rep(0,Ln_TSS),rep(as.vector(Graph_Edges[1,6]),Ln_TSS),rep("J",Ln_TSS),rep(0,Ln_TSS),stringsAsFactors = FALSE)
       colnames(TSS)<-colnames(Graph_Edges)

       TSS_Cond<-TRUE

     }

     # Extend SG

     Graph_Edges[,1]<-paste(as.vector(Graph_Edges[,1]),".b",sep="")
     Graph_Edges[,2]<-paste(as.vector(Graph_Edges[,2]),".a",sep="")
     From<-paste(as.vector(Graph_Nodes[,1]),".a",sep="")
     To<-paste(as.vector(Graph_Nodes[,1]),".b",sep="")
     Extended<-cbind(From,To,Graph_Nodes[,2:7])

     Graph_Edges<-rbind(Graph_Edges,Extended)

     if(TSS_Cond)
     {
       Graph_Edges<-rbind(TSS,Graph_Edges)
     }

     if(TTS_Cond)
     {
       Graph_Edges<-rbind(Graph_Edges,TTS)
     }


     rownames(Graph_Edges)<-1:nrow(Graph_Edges)

     # Get Adjacency and Incidence Matrix

      GGraph<-graph_from_data_frame(Graph_Edges,directed=TRUE)
      Adjacency<-as_adj(GGraph)


     Incidence<-matrix(0,nrow=((length(Nodes)*2)+2),ncol=nrow(Graph_Edges))
     colnames(Incidence)<-rownames(Graph_Edges)
     rownames(Incidence)<-c("S",paste(rep(Nodes,each=2),c(".a",".b"),sep=""),"E")

     Incidence[cbind(as.vector(Graph_Edges[,"From"]),colnames(Incidence))]<--1
     Incidence[cbind(as.vector(Graph_Edges[,"To"]),colnames(Incidence))]<-1

      iijj<-match(rownames(Incidence),rownames(Adjacency))
      Adjacency<-Adjacency[iijj,iijj]
      Adjacency<-as(Adjacency,"dgTMatrix")

     # Return All Information

      Result<-list(Edges=Graph_Edges,Adjacency=Adjacency,Incidence=Incidence)
     # Result<-list(Edges=Graph_Edges,Incidence=Incidence)
     return(Result)



  }







