test_CDFfromGTF <- function() {
  obs <- tryCatch(CDFfromGTF(input="Test"), error=conditionMessage)
  checkIdentical("Microarray field empty", obs)
  
  obs <- tryCatch(CDFfromGTF(input = "GTF", inputFile = NULL, PSR="T1", Junc="T2"
                  , PathCDF="T3", microarray = "RTA"), error=conditionMessage)
  checkIdentical("inputFile parameter is empty", obs)

}