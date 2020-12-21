#' Splicing Factor Prediction
#' 
#' Methodology to predict context-specific splicing factors
#' 
#' @param P_value_PSI A data.frame with the p.values of the experiment.
#' @param ExS The ExS matrix biuldt in CreateExSmatrix function.
#' @param nSel Top ranked events to be considered as spliced events.
#' @param significance Threshold of P.value to consider which events are 
#' deferentially spliced. A vector of length equal to the number of contrasts.
#' If null it will consider the nSel top ranked events.
#' @param method methodology to apply: "Fisher" for Fisher's exact test 
#' (default) or "PoiBin" for Poisson Binomial test.
#' 
#' @return The function returs a list. This list has for each contrast a data.frame containing
#' the results of the prediction.
#' 
#' 
#' 
#' @import glmnet
#' @import poibin
#' @import Matrix
#' @importFrom stats binomial phyper qhyper pnorm dnorm gaussian coef
#' @importFrom speedglm control is.sparse
#' @importFrom IRanges IRanges



SF_Prediction <- function(P_value_PSI,ExS,nSel=1000,significance=NULL,method="Fisher"){
  
  if(method=="Fisher"){
    ExS <- ExS[rownames(P_value_PSI), ]
    myHypers <- vector(mode="list",length = ncol(P_value_PSI))
    N <- nrow(ExS) #the same for each Fisher's test
    # one iteratio for each contrast
    for(cSel in 1:ncol(P_value_PSI)){
      
      if(is.null(significance)){
        nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])[1:nSel]
        #name of the top nSel events differentially spliced in this contrast
      }else{
        if(is.numeric(significance)){
          
          if(is.na(significance[cSel])){
            warning(paste0("no threshold selected for contrast ",cSel,". Top nSel events taken"))
            nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])[1:nSel]
          }else{
            nmTopEv <- rownames(P_value_PSI)[which(P_value_PSI[,cSel] < significance[cSel])]
            nSel <- length(nmTopEv)
          }
          
        }else{
          stop("significance must be numeric")
        }
      }
      hyperM <- data.frame(RBP = colnames(ExS),
                           nHits = colSums(ExS),
                           Pvalue_hyp_PSI =NA,
                           N = N,
                           d = colSums(ExS),
                           n = nSel,
                           x = NA,
                           qhyp_0.5 = NA,
                           Fc = NA,
                           stringsAsFactors = FALSE)
      
      mynselevents <- matrix(0,nrow=nrow(ExS),ncol=1)
      rownames(mynselevents) <- rownames(ExS)
      mynselevents[nmTopEv,1] <- 1
      
      mix <- t(ExS) %*% mynselevents
      # identical(hyperM$RBP,rownames(mix))
      mid <- hyperM$nHits
      miN_D <- N-mid
      hyperM[, "Pvalue_hyp_PSI"] <- phyper(mix, mid, miN_D, nSel, lower.tail = FALSE)
      hyperM$x <- mix
      hyperM$qhyp_0.5 <- qhyper(0.5, mid, miN_D, nSel, lower.tail = FALSE)
      hyperM$Fc <- mix/hyperM$qhyp_0.5
      
      hyperM <- hyperM[order(hyperM$Pvalue_hyp_PSI), ]
      myHypers[[cSel]]<-hyperM
      
    }
    
    return(myHypers)
  }else 
    if(method == "PoiBin"){
      ExS <- ExS[rownames(P_value_PSI), ]
      myP <- getpij(ExS)
      myHypers <- vector(mode="list",length = ncol(P_value_PSI))
      N <- nrow(ExS) #the same for each Fisher's test
      
      for(cSel in 1:ncol(P_value_PSI)){
        # cSel <- 1
        
        
        if(is.null(significance)){
          nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])[1:nSel]
          #name of the top nSel events differentially spliced in this contrast
        }else{
          if(is.numeric(significance)){
            
            if(is.na(significance[cSel])){
              warning(paste0("no threshold selected for contrast ",cSel,". Top nSel events taken"))
              nmTopEv <- (rownames(P_value_PSI)[order(P_value_PSI[,cSel])])[1:nSel]
            }else{
              nmTopEv <- rownames(P_value_PSI)[which(P_value_PSI[,cSel] < significance[cSel])]
              nSel <- length(nmTopEv)
            }
            
          }else{
            stop("significance must be numeric")
          }
        }
        
        hyperM <- data.frame(RBP = colnames(ExS),
                             nHits = colSums(ExS),
                             Pvalue_hyp_PSI =NA,
                             N = N,
                             d = colSums(ExS),
                             n = nSel,
                             x = NA,
                             qhyp_0.5 = NA,
                             Fc = NA,
                             stringsAsFactors = FALSE)
        
        mynselevents <- matrix(0,nrow=nrow(ExS),ncol=1)
        rownames(mynselevents) <- rownames(ExS)
        mynselevents[nmTopEv,1] <- 1
        
        mix <- t(ExS) %*% mynselevents
        # identical(hyperM$RBP,rownames(mix))
        # mid <- hyperM$nHits
        # miN_D <- N-mid
        # for (i in 1:ncol(ExS) ) {
        # hyperM[i, "Pvalue_hyp_PSI"] <- (1-ppoibin(kk=mix[i], pp=myP[which(rownames(myP)%in%nmTopEv),i], method="RNA"))  
        # hyperM$qhyp_0.5[i] <- qpoibin(0.5,pp = myP[which(rownames(myP)%in%nmTopEv),i])
        # }
        
        myP2 <- myP[which(rownames(myP)%in%nmTopEv),]
        muk <- colSums(myP2)
        QQ <- 1 - myP2
        sigmak <- sqrt(diag(t(myP2)%*%QQ))
        MM <- QQ*(1-2*myP2)
        gammak <- diag(t(myP2)%*%MM)
        ind = gammak/(6 * sigmak^3 )
        
        kk1 = (mix + 0.5 - muk)/sigmak
        
        vkk.r = pnorm(kk1) + ind * (1 - kk1^2) * dnorm(kk1)
        vkk.r[vkk.r < 0] = 0
        vkk.r[vkk.r > 1] = 1
        
        hyperM$Pvalue_hyp_PSI <- 1-vkk.r
        
        
        hyperM$qhyp_0.5 <- t(myP) %*% mynselevents
        
        hyperM$x <- mix
        # hyperM$qhyp_0.5 <- qhyper(0.5, mid, miN_D, nSel, lower.tail = FALSE)
        
        hyperM$Fc <- mix/hyperM$qhyp_0.5
        
        hyperM <- hyperM[order(hyperM$Pvalue_hyp_PSI), ]
        myHypers[[cSel]]<-hyperM
        
      }
      return(myHypers)
    }
  
}