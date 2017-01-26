#' Irwin Hall Summarization for Pvalues
#'
#' Summarization of pvalues according to the irwin-hall distribution
#'
#' @keywords internal
#'
#' @param Pv1 Pvalue
#' @param t1 t-statistic associated to Pv1
#' @param Pv2 Pvalue
#' @param t2 t-statistic associated to Pv2
#'
#'
#' @return Pvalues list with summarized pvalues and associated z values
#'


IHsummarization<-function(Pv1,t1,Pv2,t2, coherence = "Opposite")
{

  if (coherence == "Equal") {
    nPv1 <- (Pv1/2)*(t1>0)+(1-Pv1/2)*(t1<=0)
    nPv2 <- (Pv2/2)*(t2>0)+(1-Pv2/2)*(t2<=0)
    Psuma <- nPv1+nPv2
    PIH <- (Psuma^2)/2*(Psuma<1)+(1-(2-Psuma)^2/2)*(Psuma>=1)
    ZIH <- qnorm(PIH)
    PIH_2tail <- PIH*2 * (PIH<0.5) + ((1-PIH)*2) * (PIH>=0.5)

    # 1:  Pv1 > 0.5 ; Pv2 > 0.5 ; t1 >0 ; t2 >0

    Cambiar <- which(Pv1 < 0.5 & Pv2 < 0.5 & t1 <=0 & t2 <=0)
    nPv1[Cambiar] <- Pv1[Cambiar]/2
    nPv2[Cambiar] <- Pv2[Cambiar]/2
    Psuma[Cambiar] <- nPv1[Cambiar]+nPv2[Cambiar]
    PIH[Cambiar] <- (Psuma[Cambiar]^2)/2
    ZIH[Cambiar] <- -qnorm(PIH[Cambiar])
    PIH_2tail[Cambiar] <- (PIH[Cambiar])*2

    return(list(Pvalues=PIH_2tail,Tstats=ZIH))
  }
  else {
    return(IHsummarization(Pv1,t1,Pv2,-t2, coherence = "Equal"))
  }
}
