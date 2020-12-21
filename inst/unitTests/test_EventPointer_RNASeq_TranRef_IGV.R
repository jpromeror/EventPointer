test_EventPointer_RNASeq_TranRef_IGV <- function() {
  
  obs <- tryCatch(EventPointer_RNASeq_TranRef_IGV(SG_List = NULL,
                                                  pathtoeventstable = "",
                                                  PathGTF = ""), error=conditionMessage)
  
  checkIdentical("SG_List field is empty", obs)
  
  obs <- tryCatch(EventPointer_RNASeq_TranRef_IGV(SG_List = "",
                                                  pathtoeventstable = NULL,
                                                  PathGTF = ""), error=conditionMessage)
  
  checkIdentical("pathtoeventstable field is empty", obs)
  
  obs <- tryCatch(EventPointer_RNASeq_TranRef_IGV(SG_List = "",
                                                  pathtoeventstable = "",
                                                  PathGTF = NULL), error=conditionMessage)
  
  checkIdentical("PathGTF field is empty", obs)
  
}