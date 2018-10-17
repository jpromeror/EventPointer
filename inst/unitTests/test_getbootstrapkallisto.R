test_getbootstrapkallisto <- function() {
  
  obs <- tryCatch(getbootstrapkallisto(pathValues=NA,
                                       nb=0), error=conditionMessage)
  
  
  
  checkIdentical("not pathValues", obs)
}