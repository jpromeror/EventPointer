test_getbootstrapkallisto <- function() {
  
  obs <- tryCatch(getbootstrapkallisto(pathValues=NA,
                                       nb=20), error=conditionMessage)
  
  checkIdentical("not pathValues", obs)
  
  obs <- tryCatch(getbootstrapkallisto(pathValues="text",
                                       nb=NULL), error=conditionMessage)
  
  checkIdentical("nb field is empty", obs)
  
  obs <- tryCatch(getbootstrapkallisto(pathValues="text",
                                       nb=0), error=conditionMessage)
  
  checkIdentical("nb must be >0", obs)
}