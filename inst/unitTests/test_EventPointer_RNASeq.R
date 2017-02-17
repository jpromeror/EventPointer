test_EventPointer_RNASeq <- function() {
  
  obs <- tryCatch(EventPointer_RNASeq(Events="NULL",
                                      Design=matrix(0,nrow=5),
                                      Contrast=matrix(0,nrow=5),
                                      Statistic="Test"), error=conditionMessage)
  
  checkIdentical("Wrong statistical test provided", obs)
  
  obs <- tryCatch(EventPointer_RNASeq(Events="NULL",
                                      Design=matrix(0,nrow=5),
                                      Contrast=rep(0,5),
                                      Statistic="LogFC"), error=conditionMessage)
  

  checkIdentical("Wrong Design and/or Contrast matrices", obs)
  
  obs <- tryCatch(EventPointer_RNASeq(Events=NULL,
                                      Design=matrix(0,nrow=5),
                                      Contrast=rep(0,5),
                                      Statistic="LogFC"), error=conditionMessage)
  
  
  checkIdentical("Missing alternative splicing events", obs)
  
  
}