test_getbootstrapdata <- function() {
  
  obs <- tryCatch(getbootstrapdata(PathSamples=NULL,
                                       type=20), error=conditionMessage)
  
  checkIdentical("PathSamples field empty", obs)
  
  obs <- tryCatch(getbootstrapdata(PathSamples="text",
                                       type=NA), error=conditionMessage)
  
  checkIdentical("type field empty", obs)
  
  
  obs <- tryCatch(getbootstrapdata(PathSamples="text",
                                   type="aaa"), error=conditionMessage)
  
  checkIdentical("type must be 'kallisto' or 'salmon'.", obs)
}