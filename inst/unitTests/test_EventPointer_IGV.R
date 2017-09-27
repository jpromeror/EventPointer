test_EventPointer_IGV <- function() {
  
  obs <- tryCatch(EventPointer_IGV(Events=NULL,
                               input=NULL,
                               inputFile=NULL,
                               PSR=NULL,
                               Junc=NULL,
                               PathGTF=NULL,
                               EventsFile=NULL,
                               microarray="HTA-2_0"), error=conditionMessage)
  
  checkIdentical("No alternative splicing events provided", obs)
  
  obs <- tryCatch(EventPointer_IGV(Events="Test",
                                   input="Other",
                                   inputFile=NULL,
                                   PSR="Test",
                                   Junc="Test",
                                   PathGTF=NULL,
                                   EventsFile=NULL,
                                   microarray="RTA"), error=conditionMessage)
  
  checkIdentical("Unknown reference genome", obs)
  
 
  
  obs <- tryCatch(EventPointer_IGV(Events="Test",
                                   input="CustomGTF",
                                   inputFile=NULL,
                                   PSR="Test",
                                   Junc="Test",
                                   PathGTF=NULL,
                                   EventsFile=NULL,
                                   microarray="ClariomD"), error=conditionMessage)
  
  checkIdentical("inputFile parameter is empty", obs)
  
  obs <- tryCatch(EventPointer_IGV(Events="Test",
                                   input="AffyGTF",
                                   inputFile=NULL,
                                   PSR="Test",
                                   Junc=NULL,
                                   PathGTF=NULL,
                                   EventsFile=NULL,
                                   microarray="HTA-2_0"), error=conditionMessage)
  
  checkIdentical("Missing PSR and/or Junc files", obs)
  
}