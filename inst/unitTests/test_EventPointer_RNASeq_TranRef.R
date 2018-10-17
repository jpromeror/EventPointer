test_EventPointer_RNASeq_TranRef <- function() {
  
  obs <- tryCatch(EventPointer_RNASeq_TranRef(Count_Matrix = rep(0,5),
                                              Statistic = "LogFC",
                                              Design = matrix(0,nrow=5),
                                              Contrast = matrix(0,nrow=5)),error=conditionMessage)
    
  checkIdentical("Wrong Count_Matrix input: must be a list", obs)
  
  
  
  obs <- tryCatch(EventPointer_RNASeq_TranRef(Count_Matrix = list(),
                                              Statistic = "LogFC",
                                              Design = rep(0,5),
                                              Contrast = matrix(0,nrow=5)),error=conditionMessage)
  
  checkIdentical("Wrong Design and/or Contrast matrices", obs)
  
  
  
  obs <- tryCatch(EventPointer_RNASeq_TranRef(Count_Matrix = list(),
                                              Statistic = "Test",
                                              Design = matrix(0,nrow=5),
                                              Contrast = matrix(0,nrow=5)),error=conditionMessage)
  
  checkIdentical("Wrong statistical test given", obs)
  
  
}