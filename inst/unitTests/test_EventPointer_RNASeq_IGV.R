test_EventPointer_RNASeq_IGV <- function() {
  
  obs <- tryCatch(EventPointer_RNASeq_IGV(Events=NULL,
                                   SG_RNASeq=NULL,
                                   EventsTxt=NULL,
                                   PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Missing alternative splicing events", obs)
  
  obs <- tryCatch(EventPointer_RNASeq_IGV(Events="Test",
                                   SG_RNASeq=NULL,
                                   EventsTxt=NULL,
                                   PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Missing splicing graphs information", obs)
  
  obs <- tryCatch(EventPointer_RNASeq_IGV(Events="Test",
                                   SG_RNASeq="Test",
                                   EventsTxt=NULL,
                                   PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Wrong or missing EventsTxt field", obs)
  
  obs <- tryCatch(EventPointer_RNASeq_IGV(Events="Test",
                                   SG_RNASeq="Test",
                                   EventsTxt=123,
                                   PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Wrong or missing EventsTxt field", obs)
  
  obs <- tryCatch(EventPointer_RNASeq_IGV(Events="Test",
                                   SG_RNASeq="Test",
                                   EventsTxt="Test",
                                   PathGTF=NULL), error=conditionMessage)
  
  checkIdentical("Wrong or missing PathGTF field", obs)
  
}