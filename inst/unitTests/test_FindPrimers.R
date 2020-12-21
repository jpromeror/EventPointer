test_FindPrimers <- function() {
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=1,
                              Primer3Path="",
                              Dir="",
                              mygenomesequence="",
                              taqman = NA), error=conditionMessage)
  
  checkIdentical("taqman field is empty", obs)
  
  obs <- tryCatch(FindPrimers(SG=NULL,
                              EventNum=1,
                              Primer3Path="",
                              Dir="",
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("SG field is empty", obs)
  
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=NULL,
                              Primer3Path="",
                              Dir="",
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("EventNum field is empty or is not > 0", obs)
  
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=0,
                              Primer3Path="",
                              Dir="",
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("EventNum field is empty or is not > 0", obs)
  
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=3,
                              Primer3Path=NULL,
                              Dir="",
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("Primer3Path field is empty", obs)
  
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=3,
                              Primer3Path="",
                              Dir="/",
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("Primer3Path field should point to the .exe of Primers3.
         Check in the vignette for more information", obs)
  
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=3,
                              Primer3Path="",
                              Dir=NULL,
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("Dir field is empty", obs)
  
  
  obs <- tryCatch(FindPrimers(SG=list(),
                              EventNum=3,
                              Primer3Path="",
                              Dir="",
                              mygenomesequence="",
                              taqman = TRUE), error=conditionMessage)
  
  checkIdentical("Dir = Complete path where primer3web_v4_0_0_default_settings.txt file and primer3_config directory are stored", obs)
  
}












