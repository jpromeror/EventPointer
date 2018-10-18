#' Author: JF
#' @description  Function to load the values of the bootstrap
#' 
#' Inputs:
#' @param pathValues A vector with the complete directory to the folder of the output of kallisto
#' @param nb number of bootstrap
#'
#'
#' @return A list containing the quantification data and with the bootrap information.
#'
#'
#' @examples 
#'    
#'    PathFiles <- system.file("extdata",package="EventPointer")
#'    PathFiles <- dir(paste0(PathFiles,"/output"),full.names = TRUE)
#'    
#'    #load the data
#'    
#'    mydatab <- getbootstrapkallisto(pathValues = PathFiles,nb = 20)
#'
#'
#' @export
#' @importFrom rhdf5 h5read 



getbootstrapkallisto <- function(pathValues=NA,nb){
  
  if(is.null(nb)){
    stop("nb field is empty")
  }
  
  if(nb == 0){
    stop("nb must be >0")
  }
  
  
  if(is.na(pathValues)){
    stop("not pathValues")
  }
  
  datafinal <- vector(mode="list",length = length(pathValues))
  
  for (jj in 1:length(datafinal)){
    
    name <- unlist(strsplit(pathValues[jj],"/"))
    name <- name[length(name)]
    
    annotationreference <- h5read(paste0(pathValues[jj],"/abundance.h5"),"aux/ids")
    
    llkb <- h5read(paste0(pathValues[jj],"/abundance.h5"),"aux/eff_lengths")
    abund <- h5read(paste0(pathValues[jj],"/abundance.h5"),"/est_counts")
    
    kkk <- abund/llkb
    ppx <- which(llkb==0)
    kkk[ppx]  <- 0
    abund_tpm <- 1/sum(kkk)*kkk*1e6
    
    #read bootstrap data:
    
    for (iix in 0:(nb-1)){
      command <- paste0("abund_b",iix,
                        " <- h5read(paste0(pathValues[jj],'/abundance.h5'),'bootstrap/bs",iix,"')")
      #print(command)
      eval(parse(text = command))
      
      
      command <- paste0("kkk <- abund_b",iix,"/llkb")
      #print(command)
      eval(parse(text = command))
      
      kkk[ppx] <- 0
      
      
      command <- paste0("abund_b",iix,"_tpm <- 1/sum(kkk)*kkk*1e6")
      #print(command)
      eval(parse(text = command))
    }
    command <- paste0("losboots <- cbind(",paste(paste0("abund_b",0:19,"_tpm"),
                                                 collapse = ","),")")
    eval(parse(text = command))
    
    #add the value of the m.l (bootstrap = 0)
    if (jj ==1){ 
      #en la primera iteracion inicializar la matriz mymatrix y ponerle el rowname 
      #(el resto tendran que tener el mismo orden)
      mymatrix <- matrix(0,nrow=nrow(abund_tpm),ncol = (nb+1)) # solo inicializar la matrix 1 vez
      rownames(mymatrix) <- annotationreference
      #en la primera columna poner el abundance (bootstrap 0)
      mymatrix[,1] <- abund_tpm
      #En el resto de columnas estan los bootstrap
      mymatrix[,2:(nb+1)] <- losboots
    }else{ #comprobar que estan en el mismo orden para el resto de muestras
      if(identical(rownames(mymatrix),annotationreference)){
        mymatrix[,1] <- abund_tpm
        mymatrix[,2:(nb+1)] <- losboots
      }else{ #comprobar que estan en el mismo orden
        mymatrix[,1] <- abund_tpm[match(rownames(mymatrix),annotationreference)]
        mymatrix[,2:(nb+1)] <- losboots[match(rownames(mymatrix),annotationreference)]
      }
    }
    
    datafinal[[jj]]<-mymatrix
    names(datafinal)[jj]<-name
    
  }
  
  return(datafinal)
}


