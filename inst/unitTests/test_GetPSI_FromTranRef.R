test_GetPSI_FromTranRef <- function() {
  
  obs <- tryCatch(GetPSI_FromTranRef(PathsxTranscript = NULL,
                                     Samples = c(1,2,3),
                                     Filter=TRUE,
                                     Qn=0.25), error=conditionMessage)
  
  checkIdentical("PathsxTranscript field is empty", obs)
  
  
  obs <- tryCatch(GetPSI_FromTranRef(PathsxTranscript = c(1,2,3),
                                     Samples = NULL,
                                     Filter=TRUE,
                                     Qn=0.25), error=conditionMessage)
  
  checkIdentical("Samples field is empty", obs)
  
  
  #check example:
      data(EventXtrans)
      PathFiles <- system.file("extdata",package="EventPointer")
      filesnames <- dir(paste0(PathFiles,"/output"))
      PathFiles <- dir(paste0(PathFiles,"/output"),full.names = TRUE)
      dirtoload <- paste0(PathFiles,"/","abundance.tsv")
      RNASeq <- read.delim(dirtoload[1],sep = "\t", colClasses = c(NA,"NULL","NULL","NULL",NA))
      for (n in 2:length(dirtoload)){
        RNASeq[,n+1] <- read.delim(dirtoload[n],sep = '\t',
                                   colClasses = c('NULL','NULL','NULL','NULL',NA))
      }
      rownames(RNASeq)<-RNASeq[,1]
      RNASeq<-RNASeq[,-1]
      colnames(RNASeq) <- filesnames
      rownames(RNASeq) <- sapply(strsplit(rownames(RNASeq),"\\|"),function(X) return(X[1]))
      RNASeq<-as.matrix(RNASeq) #must be a matrix variable
      
      #Obtain values of PSI
      
      PSIss2 <- GetPSI_FromTranRef(PathsxTranscript = EventXtrans,Samples = RNASeq,Filter = FALSE)
      
      data("PSIss")
      checkIdentical(PSIss,PSIss2)
}