test_FindPrimers <- function() {
  
  obs <- tryCatch(FindPrimers(SG = NULL,
                EventNum = 4,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10, 
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  
  checkIdentical("SG field is empty", obs)
  


  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 0,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("EventNum field is empty or is not > 0", obs)

  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = -45,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("EventNum field is empty or is not > 0", obs)

  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = NULL,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("EventNum field is empty or is not > 0", obs)
  
  
  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 4,
                Primer3Path = "",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("Primer3Path field should point to the .exe of Primers3.
         Check in the vignette for more information", obs)
  
  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 4,
                Primer3Path = NULL,
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("Primer3Path field is empty", obs)
  
  
  
  
  
  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 1,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "",
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("Dir = Complete path where primer3web_v4_0_0_default_settings.txt file and primer3_config directory are stored", obs)
  
  
  
  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 1,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = NULL,
                taqman = 1,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("Dir field is empty", obs)
  
  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 1,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = NA,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("taqman field is empty", obs)
  
  
  obs <- tryCatch(
    FindPrimers(SG = list(),
                EventNum = 1,
                Primer3Path = "C:\\PROGRA~2\\primer3\\PRIMER~1.EXE",
                Dir = "C:\\PROGRA~2\\primer3\\",
                taqman = 3,
                nProbes=1,
                nPrimerstwo=4,
                ncommonForward=4,
                ncommonReverse=4,
                nExons=10,
                nPrimers =5,
                maxLength = 1200), error=conditionMessage)
  checkIdentical("taqman variable should be equal to 1 or 0", obs)
  
  
}
