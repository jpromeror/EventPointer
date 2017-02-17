test_EventPointer <- function() {
  
  obs <- tryCatch(EventPointer(Design=matrix(0,nrow=5),
                               Contrast=matrix(0,nrow=5),
                               ExFit=NULL,
                               Eventstxt=NULL,
                               Statistic="Test"), error=conditionMessage)
  
  checkIdentical("Wrong statistical test given", obs)
  
  obs <- tryCatch(EventPointer(Design=rep(0,5),
                               Contrast=matrix(0,nrow=5),
                               ExFit=NULL,
                               Eventstxt=NULL,
                               Statistic="Test"), error=conditionMessage)
  
  checkIdentical("Wrong Design and/or Contrast matrices", obs)
  
  
}