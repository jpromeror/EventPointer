test_EventDetection_transcriptome <- function() {
  
  
  #if there is not input file
  obs <- tryCatch(EventDetection_transcriptome(inputFile = NULL,
                                               Transcriptome = "text",
                                               Pathtxt="path"), error=conditionMessage)
  
  checkIdentical("not inputFile", obs)
  
  
  
  #if there is not a Transcriptome reference
  obs <- tryCatch(EventDetection_transcriptome(inputFile = "text",
                                               Transcriptome = NULL,
                                               Pathtxt="path"), error=conditionMessage)
  
  checkIdentical("not Transcriptome", obs)
  
  
  
  #if there is not a path for the .txt file
  obs <- tryCatch(EventDetection_transcriptome(inputFile = "text",
                                               Transcriptome = "text",
                                               Pathtxt=NULL), error=conditionMessage)
  
  checkIdentical("not Pathtxt", obs)
  
  
  # # check example:
  # 
  # # PathFiles<-system.file("extdata",package="EventPointer")
  # PathFiles<-system.file("extdata",package="EventPointer")
  # inputFile <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep="")
  # Transcriptome <- "Gencode24_2genes"
  # Pathtxt <- tempdir()
  # 
  # 
  # # Run the function 
  # 
  # EventXtrans2 <- EventDetection_transcriptome(inputFile = inputFile,
  #                                             Transcriptome = Transcriptome,
  #                                             Pathtxt=Pathtxt,cores=1)
  # 
  # 
  # data("EventXtrans")
  # 
  # checkIdentical(EventXtrans, EventXtrans2)
  
}