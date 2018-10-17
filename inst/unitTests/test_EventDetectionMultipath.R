test_EventDetectionMultipath <- function() {
  
  obs <- tryCatch(EventDetectionMultipath(Input=NULL,cores=1,Path="C:/Users/",paths = 2), error=conditionMessage)
  checkIdentical("Input field is empty", obs)
  
  obs <- tryCatch(EventDetectionMultipath(Input="Test",cores="O",Path="C:/Users/",paths = 2), error=conditionMessage)
  checkIdentical("Number of cores incorrect", obs)
  
  obs <- tryCatch(EventDetectionMultipath(Input="Test",cores=1,Path=NULL,paths=2), error=conditionMessage)
  checkIdentical("Path field is empty", obs)
  
  data(SG_RNASeq)
  TxtPath<-tempdir()
  AllEvents_RNASeq_MP_2<-EventDetectionMultipath(SG_RNASeq,cores=1,Path=TxtPath,paths = 3)
  data("AllEvents_RNASeq_MP")
  checkIdentical(AllEvents_RNASeq_MP_2,AllEvents_RNASeq_MP)
  
}