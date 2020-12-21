test_EventPointer_Bootstraps <- function() {
  
  obs <- tryCatch(EventPointer_Bootstraps(PSI = NULL,
                                          Design = "",
                                          Contrast = "",
                                          cores = "",
                                          ram = "",
                                          nBootstraps = "",
                                          UsePseudoAligBootstrap = "",
                                          Threshold = ""), error=conditionMessage)
  
  checkIdentical("PSI field is empty", obs)
  

  obs <- tryCatch(EventPointer_Bootstraps(PSI = "",
                                          Design = NULL,
                                          Contrast = "",
                                          cores = "",
                                          ram = "",
                                          nBootstraps = "",
                                          UsePseudoAligBootstrap = "",
                                          Threshold = ""), error=conditionMessage)
  
  checkIdentical("Design field is empty", obs)
  
  
  
  
  
  obs <- tryCatch(EventPointer_Bootstraps(PSI = "",
                                          Design = "",
                                          Contrast = NULL,
                                          cores = "",
                                          ram = "",
                                          nBootstraps = "",
                                          UsePseudoAligBootstrap = "",
                                          Threshold = ""), error=conditionMessage)
  
  checkIdentical("Contrast field is empty", obs)
  
}
