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
#'    PathFiles <- dir(paste0(PathFiles,"/output"),full.names = T)
#'    
#'    #load the data
#'    
#'    mydatab <- getbootstrapkallisto(pathValues = PathFiles,nb = 20)
#'
#'
#' @export



getbootstrapkallisto <- function(pathValues=NA,nb=0){
  
  datafinal <- vector(mode="list",length = length(pathValues))
  
  for (jj in 1:length(datafinal)){
    
    name <- unlist(strsplit(pathValues[jj],"/"))
    name <- name[length(name)]
    
    bootsdata <- read_kallisto_h5(paste0(pathValues[jj],"/abundance.h5"),read_bootstrap = TRUE,max_bootstrap = nb)
    if (jj ==1){ #en la primera iteracion inicializar la matriz mymatrix y ponerle el rowname
      mymatrix <- matrix(0,nrow=nrow(bootsdata$abundance),ncol = (nb+1)) # solo inicializar la matrix 1 vez
      rownames(mymatrix) <- bootsdata$abundance$target_id
      #en la primera columna poner el abundance (bootstrap 0)
      mymatrix[,1] <- bootsdata$abundance$tpm
    }else{ #poner el abundance cuando ya no es la primera iteracion, comprobar que estan en el mismo orden
      if(identical(rownames(mymatrix),bootsdata$abundance$target_id)){
        mymatrix[,1] <- bootsdata$abundance$tpm
      }else{ #comprobar que estan en el mismo orden
        mymatrix[,1] <- bootsdata$abundance$tpm[match(rownames(mymatrix),bootsdata$abundance$target_id)]
      }
    }
    #en cualquier iteracion hay q comprobar el orden de los datos del bootstrap
    
    mymatrix[,2:(nb+1)] <- sapply(bootsdata$bootstrap,function(X){
      if(identical(X$target_id,rownames(mymatrix))){
        return(X$tpm)
      }else{
        return(X$tpm[match(rownames(mymatrix),X$target_id)])
      }
    })
    
    
    datafinal[[jj]]<-mymatrix
    names(datafinal)[jj]<-name
    
  }
  
  return(datafinal)
}


