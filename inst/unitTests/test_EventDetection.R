test_EventDetection <- function() {
  
  obs <- tryCatch(EventDetection(Input=NULL,cores=1,Path="C:/Users/"), error=conditionMessage)
  checkIdentical("Input field is empty", obs)
  
  obs <- tryCatch(EventDetection(Input="Test",cores="O",Path="C:/Users/"), error=conditionMessage)
  checkIdentical("Number of cores incorrect", obs)
  
  obs <- tryCatch(EventDetection(Input="Test",cores=1,Path=NULL), error=conditionMessage)
  checkIdentical("Path field is empty", obs)
  
  # data(SG_RNASeq)
  # TxtPath<-tempdir()
  # AllEvents_RNASeq2<-EventDetection(SG_RNASeq,cores=1,Path=TxtPath)
  # data(AllEvents_RNASeq)
  # checkIdentical(AllEvents_RNASeq2,AllEvents_RNASeq)
  
}