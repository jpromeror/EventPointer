test_GetPSI_FromTranRef <- function() {
  
  obs <- tryCatch(GetPSI_FromTranRef(Samples="",
                                     PathsxTranscript=NULL,
                                     Bootstrap = FALSE,
                                     Filter = TRUE,
                                     Qn = 0.25), error=conditionMessage)
  
  checkIdentical("PathsxTranscript field is empty", obs)
  
  
  obs <- tryCatch(GetPSI_FromTranRef(Samples=NULL,
                                     PathsxTranscript="",
                                     Bootstrap = FALSE,
                                     Filter = TRUE,
                                     Qn = 0.25), error=conditionMessage)
  
  checkIdentical("Samples field is empty", obs)
  
  
  #check example:
  data(EventXtrans)
  # PathSamples <-system.file("extdata",package="EventPointer")
  PathSamples <-system.file("extdata",package="EventPointer")
  PathSamples <- paste0(PathSamples,"/output")
  PathSamples <- dir(PathSamples,full.names = TRUE)
  
  data_exp <- getbootstrapdata(PathSamples = PathSamples,type = "kallisto")
  
  #same annotation
  rownames(data_exp[[1]]) <- gsub("\\|.*","",rownames(data_exp[[1]]))
  #Obtain values of PSI
  
  PSIss2 <- GetPSI_FromTranRef(PathsxTranscript = EventXtrans,Samples = data_exp,Bootstrap = TRUE, Filter = FALSE)
  
  data("PSIss")
  checkIdentical(PSIss,PSIss2)
}