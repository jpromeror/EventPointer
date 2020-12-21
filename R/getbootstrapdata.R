#' getbootstrapdata
#' @description  Function to load the values of the bootstrap returned by kallisto or salmon pseudoaligners.
#'
#' @param PathSamples A vector with the complete directory to the folder of the output of kallisto/salmon.
#' @param type 'kallisto' or 'salmon'.
#'
#'
#' @return A list containing the quantification data with the bootstrap information.
#'
#'
#' @examples
#'
#'        PathSamples<-system.file("extdata",package="EventPointer")
#'        PathSamples <- paste0(PathSamples,"/output")
#'        PathSamples <- dir(PathSamples,full.names = TRUE)
#'        
#'        data_exp <- getbootstrapdata(PathSamples = PathSamples,type = "kallisto")
#'
#'
#' @export
#' @import rhdf5
#' @importFrom tximport tximport


getbootstrapdata <- function(PathSamples,
                             type){
  
  if(is.null(PathSamples)){
    stop("PathSamples field empty")
  }
  if(is.na(type)){
    stop("type field empty")
  }
  
  
  if(type == "salmon"){
    
    namesfiles <- gsub(".*/","",PathSamples)
    PathSamples <- paste0(PathSamples,"/quant.sf")
    names(PathSamples) <- namesfiles
    txi.inf.rep <- tximport(PathSamples, type = type, txOut = TRUE,countsFromAbundance="no",dropInfReps = FALSE)
    
    abundance <- txi.inf.rep$counts
    # dim(abundance)
    # dim(txi.inf.rep$length)
    totalSums <- colSums(abundance)
    abundance <- t(t(abundance)/totalSums)
    abundance <- abundance/(txi.inf.rep$length)*1e9
    # nb <- ncol(txi.inf.rep$infReps[[1]])
    data <- vector(mode = "list",length = ncol(abundance))
    names(data) <- colnames(abundance)
    for (i in seq_len(ncol(abundance))){
      # i <- 1  
      # identical(names(txi.inf.rep$infReps),colnames(abundance))
      bb <- txi.inf.rep$infReps[[i]]/txi.inf.rep$length[,i]/totalSums[i]*1e9
      data[[i]] <- cbind(abundance[,i],bb)
    }
    return(data)
    
  } else if(type == "kallisto"){
    
    namesfiles <- gsub(".*/","",PathSamples)
    PathSamples <- paste0(PathSamples,"/abundance.h5")
    names(PathSamples) <- namesfiles
    
    txi.inf.rep <- tximport(PathSamples, type = type, txOut = TRUE,countsFromAbundance="no",dropInfReps = FALSE)
    abundance <- txi.inf.rep$counts
    # dim(abundance)
    # dim(txi.inf.rep$length)
    totalSums <- colSums(abundance)
    abundance <- t(t(abundance)/totalSums)
    abundance <- abundance/(txi.inf.rep$length)*1e9
    # nb <- ncol(txi.inf.rep$infReps[[1]])
    data <- vector(mode = "list",length = ncol(abundance))
    names(data) <- colnames(abundance)
    for (i in seq_len(ncol(abundance))){
      # i <- 1  
      # identical(names(txi.inf.rep$infReps),colnames(abundance))
      bb <- txi.inf.rep$infReps[[i]]/txi.inf.rep$length[,i]/totalSums[i]*1e9
      data[[i]] <- cbind(abundance[,i],bb)
    }
    return(data)
    
  } else{
    stop("type must be 'kallisto' or 'salmon'.")
  }
  
  
}




