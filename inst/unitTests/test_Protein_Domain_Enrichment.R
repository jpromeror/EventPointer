test_Protein_Domain_Enrichment <- function() {
  
  obs <- tryCatch(Protein_Domain_Enrichment(PathsxTranscript = NULL,
                                            TxD = "",
                                            Diff_PSI = "",
                                            method = "spearman"), error=conditionMessage)
  
  checkIdentical("PathsxTranscript field is empty", obs)
  
  obs <- tryCatch(Protein_Domain_Enrichment(PathsxTranscript = "",
                                            TxD = NULL,
                                            Diff_PSI = "",
                                            method = "spearman"), error=conditionMessage)
  
  checkIdentical("TxD field is empty", obs)
  
  obs <- tryCatch(Protein_Domain_Enrichment(PathsxTranscript = "",
                                            TxD = "",
                                            Diff_PSI = NULL,
                                            method = "spearman"), error=conditionMessage)
  
  checkIdentical("Diff_PSI field is empty", obs)
  
  
}