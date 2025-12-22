test_EventPointer_RNASeq_IGV <- function() {
  
  obs <- tryCatch(EventPointerBAM_IGV(SG_RNASeq=NULL,
                                      EventsCSV=NULL,
                                      PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Missing splicing graphs information", obs)
  
  obs <- tryCatch(EventPointerBAM_IGV(SG_RNASeq="",
                                      EventsCSV=NULL,
                                      PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Wrong or missing EventsCSV field", obs)
  
  obs <- tryCatch(EventPointerBAM_IGV(SG_RNASeq="",
                                      EventsCSV="",
                                      PathGTF=NULL), error=conditionMessage)
  
  
  checkIdentical("Wrong or missing PathGTF field", obs)
  
  
}