test_EventsGTFfromTrancriptomeGTF <- function() {
  
  
  obs <- tryCatch(EventsGTFfromTrancriptomeGTF(inputFile = NULL,
                                               Transcriptome = "text",
                                               Pathtxt="path",
                                               PathGTF="path"), error=conditionMessage)
  
  checkIdentical("not inputFile", obs)
  
  
  obs <- tryCatch(EventsGTFfromTrancriptomeGTF(inputFile = "text",
                                               Transcriptome = NULL,
                                               Pathtxt="path",
                                               PathGTF="path"), error=conditionMessage)
  
  checkIdentical("not Transcriptome", obs)
  
  obs <- tryCatch(EventsGTFfromTrancriptomeGTF(inputFile = "text",
                                               Transcriptome = "text",
                                               Pathtxt=NULL,
                                               PathGTF="path"), error=conditionMessage)
  
  checkIdentical("not Pathtxt", obs)
  
  
  obs <- tryCatch(EventsGTFfromTrancriptomeGTF(inputFile = "text",
                                               Transcriptome = "text",
                                               Pathtxt="text",
                                               PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("not PathGTF", obs)
  
  # check example:
  
  PathFiles<-system.file("extdata",package="EventPointer")
  inputFile <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep="")
  Transcriptome <- "Gencode24_2genes"
  Pathtxt <- tempdir()
  PathGTF <- tempdir()
  
  # Run the function 
  
  EventXtrans2 <- EventsGTFfromTrancriptomeGTF(inputFile = inputFile,
                                              Transcriptome = Transcriptome,
                                              Pathtxt=Pathtxt,PathGTF=PathGTF)
  
  data("EventXtrans")
  
  checkIdentical(EventXtrans, EventXtrans2)
  
}