test_EventPointer_IGV <- function() {
  
  obs <- tryCatch(EventPointer_IGV(Events=NULL,
                               input=NULL,
                               inputFile=NULL,
                               PSR=NULL,
                               Junc=NULL,
                               PathGTF=NULL,
                               EventsFile=NULL), error=conditionMessage)
  
  checkIdentical("No alternative splicing events provided", obs)
  
  obs <- tryCatch(EventPointer_IGV(Events="Test",
                                   input="Other",
                                   inputFile=NULL,
                                   PSR="Test",
                                   Junc="Test",
                                   PathGTF=NULL,
                                   EventsFile=NULL), error=conditionMessage)
  
  checkIdentical("Wrong input field", obs)
  
 
  
  obs <- tryCatch(EventPointer_IGV(Events="Test",
                                   input="GTF",
                                   inputFile=NULL,
                                   PSR="Test",
                                   Junc="Test",
                                   PathGTF=NULL,
                                   EventsFile=NULL), error=conditionMessage)
  
  checkIdentical("inputFile field empty", obs)
  
  obs <- tryCatch(EventPointer_IGV(Events="Test",
                                   input="GTF",
                                   inputFile=NULL,
                                   PSR="Test",
                                   Junc=NULL,
                                   PathGTF=NULL,
                                   EventsFile=NULL), error=conditionMessage)
  
  checkIdentical("Missing PSR and/or Junc files", obs)
  
}