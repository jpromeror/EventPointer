test_PSI_Statistic <- function() {
  
  obs <- tryCatch(PSI_Statistic(PSI = NULL,
                                Design = matrix(0,nrow=5),
                                Contrast = matrix(0,nrow=5),
                                nboot = 100), error=conditionMessage)
  
  checkIdentical("PSI field is empty", obs)
  
  
  obs <- tryCatch(PSI_Statistic(PSI = matrix(0,nrow=5),
                                Design = rep(0,5),
                                Contrast = matrix(0,nrow=5),
                                nboot = 100), error=conditionMessage)
  
  checkIdentical("Wrong Design and/or Contrast matrices", obs)
  
  
  obs <- tryCatch(PSI_Statistic(PSI = matrix(0,nrow=5),
                                Design = matrix(0,nrow=5),
                                Contrast = matrix(0,nrow=5),
                                nboot = NULL), error=conditionMessage)
  
  checkIdentical("nboot field is empty", obs)
  
  
}