#' Detect splicing multipath events using EventPointer methodology
#'
#' Identification of all the multipath alternative splicing events in the splicing graphs
#'
#' @param Input Output of the PrepareBam_EP function
#' @param cores Number of cores used for parallel processing
#' @param Path Directory where to write the EventsFound_RNASeq.txt file
#' @param paths Maximum number of paths of the events to find.
#'
#'
#' @return list with all the events found for all the genes present in the experiment.
#' It also generates a file called EventsFound_RNASeq.txt with the information each event.
#'
#' @examples
#'   # Run EventDetection function
#'    data(SG_RNASeq)
#'    TxtPath<-tempdir()
#'    AllEvents_RNASeq_MP<-EventDetectionMultipath(SG_RNASeq,cores=1,Path=TxtPath,paths=3)
#'
#' @export
EventDetectionMultipath <- function(Input, 
    cores, Path, paths = 2) {
    ################################### Detect AS Events Using EP Methodology
    
    if (is.null(Input)) {
        stop("Input field is empty")
    }
    
    # if (classsss(cores) != "numeric") {
    if (!is(cores,"numeric")) {
        stop("Number of cores incorrect")
    }
    
    if (is.null(Path)) {
        stop("Path field is empty")
    }
    
    SgF <- rowRanges(Input)
    SgFC <- Input
    
    
    
    # Obtain unique genes to find alternative
    # splicing events.
    
    genes <- unique(unlist(geneID(SgF)))
    
    GeneIDs <- as.data.frame(SgF)
    GeneIDs <- GeneIDs[, c("geneID", "geneName")]
    IdName <- apply(GeneIDs, 1, function(X) {
        A <- unlist(X[2])[1]
        return(A)
    })
    GeneIDs[, 2] <- IdName
    GeneIDs <- unique(GeneIDs)
    GeneIDs[is.na(GeneIDs[, 2]), 2] <- " "
    
    
    cat("\n Obtaining Events")
    
    pb <- txtProgressBar(min = 0, max = length(genes), 
        style = 3)
    
    # registerDoMC(cores=cores) browser()
    registerDoParallel(cores = cores)
    # Result <- vector(mode='list',length =
    # length(genes)) for(jj in
    # seq_along(genes))
    Result <- foreach(jj = seq_along(genes)) %dopar% 
        {
            setTxtProgressBar(pb, jj)
            # print(jj)
            Gene <- genes[jj]
            
            SG_Gene <- SgF[geneID(SgF) == 
                Gene]
            
            iix <- match(Gene, GeneIDs[, 
                1])
            
            if (!is.na(iix)) {
                GeneName <- GeneIDs[iix, 
                  2]
            } else {
                GeneName <- " "
            }
            
            
            if (nrow(as.data.frame(SG_Gene)) != 
                1) {
                featureID(SG_Gene) <- seq_len(nrow(as.data.frame(SG_Gene)))
                
                SG_Gene_Counts <- SgFC[geneID(SgFC) == 
                  Gene]
                
                # Function to obtain the information of
                # the SG, the information returned
                # Corresponds to the Edges, Adjacency
                # matrix and Incidence Matrix of the SG
                # that is being analyzed
                
                SG <- SG_creation_RNASeq(SG_Gene)
                
                # plot(ftM2graphNEL(as.matrix(SG$Edges[,1:2]),
                # edgemode='directed'))
                
                # Using Ax=b with A the incidence matrix
                # and b=0, we solve to obtain the null
                # space and obtain the flux through each
                # edge if we relate the SG with an
                # hydraulic installation
                
                A <- solve(diag(ncol(SG$Adjacency)) - 
                  SG$Adjacency) - diag(ncol(SG$Adjacency))
                A <- A[-ncol(A), -1]
                if (any(c(A[1, ], A[, ncol(A)]) == 
                  0)) {
                  warning("\nin Gene ", Gene, 
                    " the graph is not correct\n")
                } else {
                  randSol <- getRandomFlow(SG$Incidence, 
                    ncol = 10)
                  
                  # To find events, we need to have fluxes
                  # different from 1 in different edges
                  
                  if (any(round(randSol) != 
                    1)) {
                    # Identify the alternative splicing
                    # events using the fluxes obtained
                    # previously.  The number of fluxes
                    # correspond to the number of edges so we
                    # can relate events with edges
                    
                    # Events <- findTriplets(randSol) twop <-
                    # which(rowSums(Events$triplets !=0)==3)
                    
                    Events <- findTriplets2(SG$Incidence, 
                      paths = paths, randSol)
                    
                    # There should be al least one event to
                    # continue the algorithm
                    
                    # if (nrow(Events$triplets) > 0)
                    if (nrow(Events$multipaths) > 
                      0) {
                      
                      # Transform triplets, related with edges
                      # to the corresponding edges. Also the
                      # edges are divided into Path 1, Path 2
                      # and Reference
                      
                      twopaths <- which(rowSums(Events$multipaths != 
                        0) == 4)
                      Events <- getEventMultiPaths(Events, 
                        SG, twopaths, paths)
                      
                      # Events <- getEventPaths(Events, SG)
                      
                      # Condition to ensure that there is at
                      # least one event
                      
                      if (length(Events) > 
                        0) {
                        
                        # Classify the events into canonical
                        # categories using the adjacency matrix
                        
                        Events <- ClassifyEvents(SG, 
                          Events, twopaths)
                        
                        # A simple piece of code to get the
                        # number of counts for the three paths in
                        # every event within a gene
                        
                        Events <- GetCountsMP(Events, 
                          SG_Gene_Counts)
                        
                        Events <- Filter(Negate(is.null), 
                          Events)
                        
                        Events <- lapply(Events, 
                          function(x) {
                            x$GeneName <- GeneName
                            return(x)
                          })
                        
                        Events <- lapply(Events, 
                          function(x) {
                            x$Gene <- Gene
                            return(x)
                          })
                        
                        
                        
                        Info <- AnnotateEvents_RNASeq_MultiPath(Events, 
                          paths)
                        # browser()
                        
                        torome <- filterimagine(Info, 
                          paths)
                        
                        Info <- Info[-torome, 
                          ]
                        Events <- Events[-torome]
                        
                        if (dim(Info)[1] == 
                          0) {
                          Info <- NULL
                          Events <- NULL
                        }
                        
                        return(list(Events = Events, 
                          Info = Info))
                        # Result[[jj]] <- list(Events = Events,
                        # Info = Info)
                        
                        # Eventss[[jj]]<-Events
                        # Infoss[[jj]]<-Info
                      }
                    }
                  }
                }
            }
        }
    
    close(pb)
    cat("\n")
    
    Result <- Filter(Negate(is.null), Result)
    
    Eventss <- vector("list", length = length(Result))
    Infoss <- vector("list", length = length(Result))
    
    for (jj in seq_len(length(Result))) {
        Eventss[[jj]] <- Result[[jj]]$Events
        Infoss[[jj]] <- Result[[jj]]$Info
    }
    
    
    Events <- Filter(Negate(is.null), Eventss)
    TxtInfo <- Filter(Negate(is.null), Infoss)
    
    TxtInfo <- do.call(rbind, TxtInfo)
    iix <- which(TxtInfo[, 2] == " ")
    
    if (length(iix) > 0) {
        Info <- matrix(unlist(strsplit(as.vector(TxtInfo[iix, 
            1]), "_")), ncol = 2, byrow = TRUE)[, 
            1]
        TxtInfo[iix, 2] <- Info
    }
    
    write.table(TxtInfo, file = paste(Path, 
        "/EventsFound_RNASeq.txt", sep = ""), 
        sep = "\t", row.names = FALSE, col.names = TRUE, 
        quote = FALSE)
    
    return(Events)
}
