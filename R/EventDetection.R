#' Detect splicing events using EventPointer methodology
#'
#' Identification of all the alternative splicing events in the splicing graphs
#'
#' @param Input Output of the PrepareBam_EP function
#' @param cores Number of cores used for parallel processing
#' @param Path Directory where to write the EventsFound_RNASeq.txt file
#'
#'
#' @return list with all the events found for all the genes present in the experiment.
#' It also generates a file called EventsFound_RNASeq.txt with the information of each event.
#'
#' @examples
#'   # Run EventDetection function
#'    data(SG_RNASeq)
#'    TxtPath<-tempdir()
#'    AllEvents_RNASeq<-EventDetection(SG_RNASeq,cores=1,Path=TxtPath)
#'
#' @export

EventDetection <- function(Input, cores, Path) {
    ################################### Detect AS Events Using EP Methodology
    
  if(is.null(Input))
  {
    stop("Input field is empty")
  }
  
  if(class(cores)!="numeric")
  {
    stop("Number of cores incorrect")
  }
  
  if(is.null(Path))
  {
    stop("Path field is empty")
  }
  
    SgF <- rowRanges(Input)
    SgFC <- Input
    
    # Obtain unique genes to find alternative splicing events.
    
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
    
    pb <- txtProgressBar(min = 0, max = length(genes), style = 3)
    
    # registerDoMC(cores=cores)
    registerDoParallel(cores = cores)
    # for(jj in 1:length(genes))
    Result <- foreach(jj = seq_along(genes)) %dopar% 
    {
        setTxtProgressBar(pb, jj)
         # print(jj)
        Gene <- genes[jj]
        
        SG_Gene <- SgF[geneID(SgF) == Gene]
        
        iix <- match(Gene, GeneIDs[, 1])
        
        if (!is.na(iix)) {
            GeneName <- GeneIDs[iix, 2]
        } else {
            GeneName <- " "
        }
        
        
        if (nrow(as.data.frame(SG_Gene)) != 1) {
            featureID(SG_Gene) <- 1:nrow(as.data.frame(SG_Gene))
            
            SG_Gene_Counts <- SgFC[geneID(SgFC) == Gene]
            
            # Function to obtain the information of the SG, the information returned
            # Corresponds to the Edges, Adjacency matrix and Incidence Matrix of the SG that
            # is being analyzed
            SG <- SG_Info(SG_Gene)
            
            # Using Ax=b with A the incidence matrix and b=0, we solve to obtain the null
            # space and obtain the flux through each edge if we relate the SG with an
            # hydraulic installation
            
            randSol <- getRandomFlow(SG$Incidence, ncol = 10)
            
            # To find events, we need to have fluxes different from 1 in different edges
            
            if (any(round(randSol) != 1)) {
                # Identify the alternative splicing events using the fluxes obtained previously.
                # The number of fluxes correspond to the number of edges so we can relate events
                # with edges
                
                Events <- findTriplets(randSol)
                
                # There should be al least one event to continue the algorithm
                
                if (nrow(Events$triplets) > 0) {
                  twopaths <- which(rowSums(Events$triplets != 0)==3)
                  # Transform triplets, related with edges to the corresponding edges. Also the
                  # edges are divided into Path 1, Path 2 and Reference
                  
                  Events <- getEventPaths(Events, SG)
                  
                  # Condition to ensure that there is at least one event
                  
                  if (length(Events) > 0) {
                    
                    # Classify the events into canonical categories using the adjacency matrix
                    
                    Events <- ClassifyEvents(SG, Events,twopaths)
                    
                    # A simple piece of code to get the number of counts for the three paths 
                    # in every
                    # event within a gene
                    
                    Events <- GetCounts(Events, SG_Gene_Counts)
                    
                    Events <- Filter(Negate(is.null), Events)
                    
                    Events <- lapply(Events, function(x) {
                      x$GeneName <- GeneName
                      return(x)
                    })
                    
                    Events <- lapply(Events, function(x) {
                      x$Gene <- Gene
                      return(x)
                    })
                    
                    
                    
                    Info <- AnnotateEvents_RNASeq(Events)
                    # browser()
                    
                    return(list(Events = Events, Info = Info))
                    
                    
                  }
                  
                }
                
            }
            
        }
        
    }
    
    close(pb)
    cat("\n")
    
    Result <- Filter(Negate(is.null), Result)
    
    Events <- vector("list", length = length(Result))
    TxtInfo <- vector("list", length = length(Result))
    
    for (jj in 1:length(Result)) {
        Events[[jj]] <- Result[[jj]]$Events
        TxtInfo[[jj]] <- Result[[jj]]$Info
    }
    
    # browser()
    TxtInfo <- do.call(rbind, TxtInfo)
    iix <- which(TxtInfo[, 2] == " ")
    
    if (length(iix) > 0) {
        Info <- matrix(unlist(strsplit(as.vector(TxtInfo[iix, 1]), "_")), ncol = 2, 
            byrow = TRUE)[, 1]
        TxtInfo[iix, 2] <- Info
        
    }
    
    write.table(TxtInfo, file = paste(Path, "/EventsFound_RNASeq.txt", sep = ""), 
        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    return(Events)
    
}
