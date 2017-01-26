# Identfication of coupled fluxes and triplets
#
#' Internal function called by EvenPointer
#'
#' Function to identify the triplets that represent alternative splicing events.
#'
#' @keywords internal
#'
#' @param randSol Flux traversing the edges of the splicing graph
#' @param tol Tolerance used
#'
#'
#' @return list with all the triplets detected for a particular gene
#'


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
