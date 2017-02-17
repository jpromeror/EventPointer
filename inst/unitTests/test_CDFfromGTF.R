test_CDFfromGTF <- function() {
  obs <- tryCatch(CDFfromGTF(input="Test"), error=conditionMessage)
  checkIdentical("Unknown reference genome", obs)
  
  obs <- tryCatch(CDFfromGTF(input="GTF"), error=conditionMessage)
  checkIdentical("inputFile parameter is empty", obs)

}