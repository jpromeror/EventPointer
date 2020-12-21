test_ResulTable <- function() {
  
  obs <- tryCatch(ResulTable(EP_Result = NULL,
                             coef=1,
                             number=Inf), error=conditionMessage)
  
  checkIdentical("EP_Result field is empty", obs)
  
}