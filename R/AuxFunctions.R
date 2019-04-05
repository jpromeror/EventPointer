#' EventPointer Internal Functions
#'
#' Internal functions used by EventPointer in the different
#' steps of the algorithm
#'
#' @keywords internal
#' @name InternalFunctions
#' @return Internal outputs
#'
#'
NULL

#' @rdname InternalFunctions
annotateEvents <- function(Events, PSR_Gene, 
    Junc_Gene, Gxx) {
    GeneName <- Gxx
    ENSGID <- Gxx
    Chrom <- gsub("chr", "", as.vector(Events[[1]]$P1[1, 
        "Chr"]))
    Result <- vector("list")
    Flat <- vector("list")
    mm <- 0
    
    for (ii in seq_along(Events)) {
        Events[[ii]]$Probes_P1 <- NULL
        Events[[ii]]$Probes_P2 <- NULL
        Events[[ii]]$Probes_Ref <- NULL
        PSR_P1 <- c()
        Junc_P1 <- c()
        PSR_P2 <- c()
        Junc_P2 <- c()
        PSR_Ref <- c()
        Junc_Ref <- c()
        
        ExonsP1 <- which(Events[[ii]]$P1[, 
            "Type"] == "E")
        JunctionsP1 <- which(Events[[ii]]$P1[, 
            "Type"] == "J")
        
        if (length(ExonsP1) > 0) {
            EPP1 <- Events[[ii]]$P1[ExonsP1, 
                ]
            PSR_P1 <- sapply(seq_len(nrow(EPP1)), 
                function(x) {
                  which(as.numeric(PSR_Gene[, 
                    "Start"]) >= as.numeric(EPP1[x, 
                    "Start"]) & as.numeric(PSR_Gene[, 
                    "Stop"]) <= as.numeric(EPP1[x, 
                    "End"]))
                })
            PSR_P1 <- PSR_Gene[unlist(PSR_P1), 
                1]
        }
        
        if (length(JunctionsP1) > 0) {
            JPP1 <- Events[[ii]]$P1[JunctionsP1, 
                ]
            Junc_P1 <- sapply(seq_len(nrow(JPP1)), 
                function(x) {
                  which(as.numeric(Junc_Gene[, 
                    "Start"]) == as.numeric(JPP1[x, 
                    "Start"]) & as.numeric(Junc_Gene[, 
                    "Stop"]) == as.numeric(JPP1[x, 
                    "End"]))
                })
            Junc_P1 <- Junc_Gene[unlist(Junc_P1), 
                1]
        }
        
        
        ExonsP2 <- which(Events[[ii]]$P2[, 
            "Type"] == "E")
        JunctionsP2 <- which(Events[[ii]]$P2[, 
            "Type"] == "J")
        
        if (length(ExonsP2) > 0) {
            EPP2 <- Events[[ii]]$P2[ExonsP2, 
                ]
            PSR_P2 <- sapply(seq_len(nrow(EPP2)), 
                function(x) {
                  which(as.numeric(PSR_Gene[, 
                    "Start"]) >= as.numeric(EPP2[x, 
                    "Start"]) & as.numeric(PSR_Gene[, 
                    "Stop"]) <= as.numeric(EPP2[x, 
                    "End"]))
                })
            PSR_P2 <- PSR_Gene[unlist(PSR_P2), 
                1]
        }
        
        if (length(JunctionsP2) > 0) {
            JPP2 <- Events[[ii]]$P2[JunctionsP2, 
                ]
            Junc_P2 <- sapply(seq_len(nrow(JPP2)), 
                function(x) {
                  which(as.numeric(Junc_Gene[, 
                    "Start"]) == as.numeric(JPP2[x, 
                    "Start"]) & as.numeric(Junc_Gene[, 
                    "Stop"]) == as.numeric(JPP2[x, 
                    "End"]))
                })
            Junc_P2 <- Junc_Gene[unlist(Junc_P2), 
                1]
        }
        
        
        ExonsRef <- which(Events[[ii]]$Ref[, 
            "Type"] == "E")
        JunctionsRef <- which(Events[[ii]]$Ref[, 
            "Type"] == "J")
        
        if (length(ExonsRef) > 0) {
            EPRef <- Events[[ii]]$Ref[ExonsRef, 
                ]
            PSR_Ref <- sapply(seq_len(nrow(EPRef)), 
                function(x) {
                  which(as.numeric(PSR_Gene[, 
                    "Start"]) >= as.numeric(EPRef[x, 
                    "Start"]) & as.numeric(PSR_Gene[, 
                    "Stop"]) <= as.numeric(EPRef[x, 
                    "End"]))
                })
            PSR_Ref <- PSR_Gene[unlist(PSR_Ref), 
                1]
        }
        
        if (length(JunctionsRef) > 0) {
            JPRef <- Events[[ii]]$Ref[JunctionsRef, 
                ]
            Junc_Ref <- sapply(seq_len(nrow(JPRef)), 
                function(x) {
                  which(as.numeric(Junc_Gene[, 
                    "Start"]) == as.numeric(JPRef[x, 
                    "Start"]) & as.numeric(Junc_Gene[, 
                    "Stop"]) == as.numeric(JPRef[x, 
                    "End"]))
                })
            Junc_Ref <- Junc_Gene[unlist(Junc_Ref), 
                1]
        }
        
        
        Events[[ii]]$Probes_P1 <- c(PSR_P1, 
            Junc_P1)
        Events[[ii]]$Probes_P2 <- c(PSR_P2, 
            Junc_P2)
        Events[[ii]]$Probes_Ref <- c(PSR_Ref, 
            Junc_Ref)
        
        if (length(Events[[ii]]$Probes_P1) > 
            0 & length(Events[[ii]]$Probes_P2) > 
            0 & length(Events[[ii]]$Probes_Ref) > 
            0) {
            mm <- mm + 1
            
            EventNumber <- ii
            
            EventType <- Events[[ii]]$Type
            
            Positions <- rbind(Events[[ii]]$P1, 
                Events[[ii]]$P2)[, 4:5]
            Start <- as.numeric(Positions[, 
                1])
            End <- as.numeric(Positions[, 
                2])
            Start <- Start[which(Start != 
                0)]
            End <- End[which(End != 0)]
            
            # browser()
            minGPos <- min(Start)
            maxGPos <- max(End)
            GPos <- paste(Chrom, ":", minGPos, 
                "-", maxGPos, sep = "")
            
            CP1s <- which(Events[[ii]]$P1[, 
                1] == "S")
            CP1e <- which(Events[[ii]]$P1[, 
                2] == "E")
            
            if (length(CP1s) > 0 | length(CP1e) > 
                0) {
                CC <- c(CP1s, CP1e)
                Events[[ii]]$P1 <- Events[[ii]]$P1[-CC, 
                  ]
            }
            
            PS1 <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$P1[, 1]))
            PE1 <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$P1[, 2]))
            Path1 <- as.matrix(cbind(PS1, 
                PE1))
            Path1 <- Path1[order(Path1[, 
                1], Path1[, 2]), , drop = FALSE]
            
            CP2s <- which(Events[[ii]]$P2[, 
                1] == "S")
            CP2e <- which(Events[[ii]]$P2[, 
                2] == "E")
            
            if (length(CP2s) > 0 | length(CP2e) > 
                0) {
                CC <- c(CP2s, CP2e)
                Events[[ii]]$P2 <- Events[[ii]]$P2[-CC, 
                  ]
            }
            
            PS2 <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$P2[, 1]))
            PE2 <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$P2[, 2]))
            Path2 <- as.matrix(cbind(PS2, 
                PE2))
            Path2 <- Path2[order(Path2[, 
                1], Path2[, 2]), , drop = FALSE]
            
            CPRs <- which(Events[[ii]]$Ref[, 
                1] == "S")
            CPRe <- which(Events[[ii]]$Ref[, 
                2] == "E")
            
            if (length(CPRs) > 0 | length(CPRe) > 
                0) {
                CC <- c(CPRs, CPRe)
                Events[[ii]]$Ref <- Events[[ii]]$Ref[-CC, 
                  ]
            }
            
            PSR <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$Ref[, 1]))
            PER <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$Ref[, 2]))
            PathR <- as.matrix(cbind(PSR, 
                PER))
            PathR <- PathR[order(PathR[, 
                1], PathR[, 2]), , drop = FALSE]
            
            
            Path1 <- paste(Path1[, 1], "-", 
                Path1[, 2], sep = "", collapse = ",")
            Path2 <- paste(Path2[, 1], "-", 
                Path2[, 2], sep = "", collapse = ",")
            PathR <- paste(PathR[, 1], "-", 
                PathR[, 2], sep = "", collapse = ",")
            
            ProbesP1 <- paste(Events[[ii]]$Probes_P1, 
                collapse = ",")
            ProbesP2 <- paste(Events[[ii]]$Probes_P2, 
                collapse = ",")
            ProbesR <- paste(Events[[ii]]$Probes_Ref, 
                collapse = ",")
            
            NEv <- data.frame(GeneName, ENSGID, 
                EventNumber, EventType, GPos, 
                Path1, Path2, PathR, ProbesP1, 
                ProbesP2, ProbesR, stringsAsFactors = FALSE)
            Result[[mm]] <- NEv
            
            
            Tprobes <- rbind(PSR_Gene, Junc_Gene)
            ii.P1 <- match(Events[[ii]]$Probes_P1, 
                Tprobes[, 1])
            ii.P2 <- match(Events[[ii]]$Probes_P2, 
                Tprobes[, 1])
            ii.R <- match(Events[[ii]]$Probes_Ref, 
                Tprobes[, 1])
            
            lP1 <- length(ii.P1)
            lP2 <- length(ii.P2)
            lRef <- length(ii.R)
            
            xRef <- rep(paste(GeneName, "_", 
                EventNumber, "_Ref", sep = ""), 
                lRef)
            xP1 <- rep(paste(GeneName, "_", 
                EventNumber, "_P1", sep = ""), 
                lP1)
            xP2 <- rep(paste(GeneName, "_", 
                EventNumber, "_P2", sep = ""), 
                lP2)
            xTot <- rep(paste(GeneName, "_", 
                EventNumber, sep = ""), lP1 + 
                lP2 + lRef)
            
            AllProbes <- c(Events[[ii]]$Probes_Ref, 
                Events[[ii]]$Probes_P1, Events[[ii]]$Probes_P2)
            flat_gene <- cbind(AllProbes, 
                Tprobes[c(ii.R, ii.P1, ii.P2), 
                  c(2, 3, 9)], c(xRef, xP1, 
                  xP2), xTot)
            
            Flat[[mm]] <- flat_gene
        }
    }
    
    Result <- do.call(rbind, Result)
    Flat <- do.call(rbind, Flat)
    return(list(Events = Result, Flat = Flat))
}


#' @rdname InternalFunctions
annotateEventsMultipath <- function(Events, 
    PSR_Gene, Junc_Gene, Gxx, paths) {
    GeneName <- Gxx
    ENSGID <- Gxx
    Chrom <- gsub("chr", "", as.vector(Events[[1]]$P1[1, 
        "Chr"]))
    Result <- vector("list")
    Flat <- vector("list")
    mm <- 0
    
    for (ii in seq_along(Events)) {
        for (kk in seq_len(paths)) {
            command <- paste0("Events[[ii]]$Probes_P", 
                kk, "<-NULL")
            eval(parse(text = command))
            command <- paste0("PSR_P", kk, 
                "<-c()")
            eval(parse(text = command))
            command <- paste0("Junc_P", kk, 
                "<-c()")
            eval(parse(text = command))
        }
        Events[[ii]]$Probes_Ref <- NULL
        PSR_Ref <- c()
        Junc_Ref <- c()
        
        
        for (kk in seq_len(paths)) {
            command <- paste0("ExonsP", kk, 
                "<-which(Events[[ii]]$P", 
                kk, "[,'Type']=='E')")
            eval(parse(text = command))
            
            command <- paste0("JunctionsP", 
                kk, "<-which(Events[[ii]]$P", 
                kk, "[,'Type']=='J')")
            eval(parse(text = command))
            
            # ExonsP1,2,3 and JunctionP1,2,... etc
            
            command <- paste0("a <- length(ExonsP", 
                kk, ">0)")
            eval(parse(text = command))
            if (a > 0) {
                command <- paste0("EPP", 
                  kk, "<-Events[[ii]]$P", 
                  kk, "[ExonsP", kk, ",]")
                eval(parse(text = command))
                command <- paste0("PSR_P", 
                  kk, "<-sapply(1:nrow(EPP", 
                  kk, "),function(x){which(as.numeric(PSR_Gene[,'Start'])>=as.numeric(EPP", 
                  kk, "[x,'Start']) & as.numeric(PSR_Gene[,'Stop'])<=as.numeric(EPP", 
                  kk, "[x,'End']))})")
                eval(parse(text = command))
                command <- paste0("PSR_P", 
                  kk, "<-PSR_Gene[unlist(PSR_P", 
                  kk, "),1]")
                eval(parse(text = command))
            }
            
            command <- paste0("a <- length(JunctionsP", 
                kk, ">0)")
            eval(parse(text = command))
            if (a > 0) {
                command <- paste0("JPP", 
                  kk, "<-Events[[ii]]$P", 
                  kk, "[JunctionsP", kk, 
                  ",]")
                eval(parse(text = command))
                command <- paste0("Junc_P", 
                  kk, "<-sapply(1:nrow(JPP", 
                  kk, "),function(x){which(as.numeric(Junc_Gene[,'Start'])==as.numeric(JPP", 
                  kk, "[x,'Start']) & as.numeric(Junc_Gene[,'Stop'])==as.numeric(JPP", 
                  kk, "[x,'End']))})")
                eval(parse(text = command))
                command <- paste0("Junc_P", 
                  kk, "<-Junc_Gene[unlist(Junc_P", 
                  kk, "),1]")
                eval(parse(text = command))
            }
        }
        
        ExonsRef <- which(Events[[ii]]$Ref[, 
            "Type"] == "E")
        JunctionsRef <- which(Events[[ii]]$Ref[, 
            "Type"] == "J")
        
        if (length(ExonsRef) > 0) {
            EPRef <- Events[[ii]]$Ref[ExonsRef, 
                ]
            PSR_Ref <- sapply(seq_len(nrow(EPRef)), 
                function(x) {
                  which(as.numeric(PSR_Gene[, 
                    "Start"]) >= as.numeric(EPRef[x, 
                    "Start"]) & as.numeric(PSR_Gene[, 
                    "Stop"]) <= as.numeric(EPRef[x, 
                    "End"]))
                })
            PSR_Ref <- PSR_Gene[unlist(PSR_Ref), 
                1]
        }
        
        if (length(JunctionsRef) > 0) {
            JPRef <- Events[[ii]]$Ref[JunctionsRef, 
                ]
            Junc_Ref <- sapply(seq_len(nrow(JPRef)), 
                function(x) {
                  which(as.numeric(Junc_Gene[, 
                    "Start"]) == as.numeric(JPRef[x, 
                    "Start"]) & as.numeric(Junc_Gene[, 
                    "Stop"]) == as.numeric(JPRef[x, 
                    "End"]))
                })
            Junc_Ref <- Junc_Gene[unlist(Junc_Ref), 
                1]
        }
        
        
        for (kk in seq_len(paths)) {
            command <- paste0("Events[[ii]]$Probes_P", 
                kk, "<-c(PSR_P", kk, ",Junc_P", 
                kk, ")")
            eval(parse(text = command))
        }
        Events[[ii]]$Probes_Ref <- c(PSR_Ref, 
            Junc_Ref)
        
        
        # only the events in which all their
        # events are able to be measured are
        # shown. It is necesary to know the
        # number of paths of each Event
        for (kk in seq_len((Events[[ii]]$NumP + 
            1))) {
            if (kk == 1) {
                a <- paste0("a <- length(Events[[ii]]$Probes_P", 
                  kk, ")>0 & ")
            } else if (kk == (Events[[ii]]$NumP + 
                1)) {
                a <- paste0(a, "length(Events[[ii]]$Probes_Ref)>0")
            } else {
                a <- paste0(a, "length(Events[[ii]]$Probes_P", 
                  kk, ")>0 & ")
            }
        }
        eval(parse(text = a))
        
        
        
        
        if (a) {
            mm <- mm + 1
            
            EventNumber <- ii
            
            EventType <- Events[[ii]]$Type
            
            EventNumP <- Events[[ii]]$NumP
            for (kk in seq_len(EventNumP)) {
                if (kk == 1) {
                  Positions <- paste0("Positions <- rbind(Events[[ii]]$P", 
                    kk)
                } else if (kk == EventNumP) {
                  Positions <- paste0(Positions, 
                    ",Events[[ii]]$P", kk, 
                    ")[,4:5]")
                } else {
                  Positions <- paste0(Positions, 
                    ",Events[[ii]]$P", kk)
                }
            }
            eval(parse(text = Positions))
            # Positions<-rbind(Events[[ii]]$P1,Events[[ii]]$P2)[,4:5]
            Start <- as.numeric(Positions[, 
                1])
            End <- as.numeric(Positions[, 
                2])
            Start <- Start[which(Start != 
                0)]
            End <- End[which(End != 0)]
            
            # browser()
            minGPos <- min(Start)
            maxGPos <- max(End)
            GPos <- paste(Chrom, ":", minGPos, 
                "-", maxGPos, sep = "")
            
            for (kk in seq_len(EventNumP)) {
                command <- paste0("CP", kk, 
                  "s<-which(Events[[ii]]$P", 
                  kk, "[,1]=='S')")
                eval(parse(text = command))
                command <- paste0("CP", kk, 
                  "e<-which(Events[[ii]]$P", 
                  kk, "[,2]=='E')")
                eval(parse(text = command))
                a <- paste0("a<-length(CP", 
                  kk, "s)>0|length(CP", kk, 
                  "e)>0")
                eval(parse(text = a))
                if (a) {
                  command <- paste0("CC<-c(CP", 
                    kk, "s,CP", kk, "e)")
                  eval(parse(text = command))
                  command <- paste0("Events[[ii]]$P", 
                    kk, "<-Events[[ii]]$P", 
                    kk, "[-CC,]")
                  eval(parse(text = command))
                }
                command <- paste0("PS", kk, 
                  "<-as.numeric(gsub('.[ab]','',Events[[ii]]$P", 
                  kk, "[,1]))")
                eval(parse(text = command))
                command <- paste0("PE", kk, 
                  "<-as.numeric(gsub('.[ab]','',Events[[ii]]$P", 
                  kk, "[,2]))")
                eval(parse(text = command))
                command <- paste0("Path", 
                  kk, "<-as.matrix(cbind(PS", 
                  kk, ",PE", kk, "))")
                eval(parse(text = command))
                command <- paste0("Path", 
                  kk, "<-Path", kk, "[order(Path", 
                  kk, "[,1],Path", kk, "[,2]),,drop=FALSE]")
                eval(parse(text = command))
            }
            
            
            CPRs <- which(Events[[ii]]$Ref[, 
                1] == "S")
            CPRe <- which(Events[[ii]]$Ref[, 
                2] == "E")
            
            if (length(CPRs) > 0 | length(CPRe) > 
                0) {
                CC <- c(CPRs, CPRe)
                Events[[ii]]$Ref <- Events[[ii]]$Ref[-CC, 
                  ]
            }
            
            PSR <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$Ref[, 1]))
            PER <- as.numeric(gsub(".[ab]", 
                "", Events[[ii]]$Ref[, 2]))
            PathR <- as.matrix(cbind(PSR, 
                PER))
            PathR <- PathR[order(PathR[, 
                1], PathR[, 2]), , drop = FALSE]
            
            PathR <- paste(PathR[, 1], "-", 
                PathR[, 2], sep = "", collapse = ",")
            ProbesR <- paste(Events[[ii]]$Probes_Ref, 
                collapse = ",")
            NEv <- "NEv<-data.frame(GeneName,ENSGID,EventNumber,EventType,GPos,EventNumP,"
            for (kk in seq_len(paths)) {
                if (kk <= EventNumP) {
                  command <- paste0("Path", 
                    kk, "<-paste(Path", kk, 
                    "[,1],'-',Path", kk, 
                    "[,2],sep='',collapse=',')")
                  eval(parse(text = command))
                  command <- paste0("ProbesP", 
                    kk, "<-paste(Events[[ii]]$Probes_P", 
                    kk, ",collapse=',')")
                  eval(parse(text = command))
                } else {
                  command <- paste0("Path", 
                    kk, "<-'-'")
                  eval(parse(text = command))
                  command <- paste0("ProbesP", 
                    kk, "<-'-'")
                  eval(parse(text = command))
                }
                NEv <- paste0(NEv, "Path", 
                  kk, ",")
            }
            NEv <- paste0(NEv, "PathR,")
            for (kk in seq_len(paths)) {
                NEv <- paste0(NEv, "ProbesP", 
                  kk, ",")
            }
            NEv <- paste0(NEv, "ProbesR,stringsAsFactors = FALSE)")
            eval(parse(text = NEv))
            # NEv<-data.frame(GeneName,ENSGID,EventNumber,EventType,GPos,
            # Path1,Path2,PathR,ProbesP1,ProbesP2,ProbesR,stringsAsFactors
            # = FALSE)
            Result[[mm]] <- NEv
            
            
            Tprobes <- rbind(PSR_Gene, Junc_Gene)
            xTot <- "xTot<-rep(paste(GeneName,'_',EventNumber,sep=''),"
            AllProbes <- "AllProbes<-c(Events[[ii]]$Probes_Ref,"
            flat_gene <- "flat_gene<-cbind(AllProbes,Tprobes[c(ii.R,"
            for (kk in seq_len(EventNumP)) {
                command <- paste0("ii.P", 
                  kk, "<-match(Events[[ii]]$Probes_P", 
                  kk, ",Tprobes[,1])")
                eval(parse(text = command))
                command <- paste0("lP", kk, 
                  "<-length(ii.P", kk, ")")
                eval(parse(text = command))
                command <- paste0("xP", kk, 
                  "<-rep(paste(GeneName,'_',EventNumber,'_P", 
                  kk, "',sep=''),lP", kk, 
                  ")")
                eval(parse(text = command))
                xTot <- paste0(xTot, "lP", 
                  kk, "+")
                if (kk == EventNumP) {
                  AllProbes <- paste0(AllProbes, 
                    "Events[[ii]]$Probes_P", 
                    kk, ")")
                  flat_gene <- paste0(flat_gene, 
                    "ii.P", kk, "),c(2,3,9)],c(xRef,")
                  for (zz in seq_len(EventNumP)) {
                    if (zz == EventNumP) {
                      flat_gene <- paste0(flat_gene, 
                        "xP", zz, "),xTot)")
                    } else {
                      flat_gene <- paste0(flat_gene, 
                        "xP", zz, ",")
                    }
                  }
                } else {
                  AllProbes <- paste0(AllProbes, 
                    "Events[[ii]]$Probes_P", 
                    kk, ",")
                  flat_gene <- paste0(flat_gene, 
                    "ii.P", kk, ",")
                }
            }
            xTot <- paste0(xTot, "lRef)")
            
            ii.R <- match(Events[[ii]]$Probes_Ref, 
                Tprobes[, 1])
            lRef <- length(ii.R)
            xRef <- rep(paste(GeneName, "_", 
                EventNumber, "_Ref", sep = ""), 
                lRef)
            
            eval(parse(text = xTot))
            eval(parse(text = AllProbes))
            eval(parse(text = flat_gene))
            # AllProbes<-c(Events[[ii]]$Probes_Ref,Events[[ii]]$Probes_P1,
            # Events[[ii]]$Probes_P2)
            # flat_gene<-cbind(AllProbes,Tprobes[c(ii.R,ii.P1,ii.P2),c(2,3,9)],
            # c(xRef,xP1,xP2),xTot)
            colnames(flat_gene) <- c("Probe_ID", 
                "X", "Y", "Probe_Sequence", 
                "Group_ID", "Unit_ID")
            Flat[[mm]] <- flat_gene
        }
    }
    
    Result <- do.call(rbind, Result)
    Flat <- do.call(rbind, Flat)
    return(list(Events = Result, Flat = Flat))
}


#' @rdname InternalFunctions
AnnotateEvents_RNASeq <- function(Events) {
    Result <- vector("list", length = length(Events))
    for (ii in seq_along(Events)) {
        GeneName <- as.vector(Events[[ii]]$GeneName)
        GeneID <- as.vector(Events[[ii]]$Gene)
        EventNumber <- ii
        EventID <- paste(GeneID, "_", EventNumber, 
            sep = "")
        EventType <- Events[[ii]]$Type
        Chrom <- as.vector(Events[[ii]]$P1[1, 
            "Chr"])
        Positions <- rbind(Events[[ii]]$P1, 
            Events[[ii]]$P2)[, 4:5]
        Start <- as.numeric(Positions[, 1])
        End <- as.numeric(Positions[, 2])
        Start <- Start[which(Start != 0)]
        End <- End[which(End != 0)]
        
        # browser()
        minGPos <- min(Start)
        maxGPos <- max(End)
        GPos <- paste(Chrom, ":", minGPos, 
            "-", maxGPos, sep = "")
        
        CP1s <- which(Events[[ii]]$P1[, 1] == 
            "S")
        CP1e <- which(Events[[ii]]$P1[, 2] == 
            "E")
        
        if (length(CP1s) > 0 | length(CP1e) > 
            0) {
            CC <- c(CP1s, CP1e)
            Events[[ii]]$P1 <- Events[[ii]]$P1[-CC, 
                ]
        }
        
        PS1 <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$P1[, 1]))
        PE1 <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$P1[, 2]))
        Path1 <- as.matrix(cbind(PS1, PE1))
        Path1 <- Path1[order(Path1[, 1], 
            Path1[, 2]), , drop = FALSE]
        
        CP2s <- which(Events[[ii]]$P2[, 1] == 
            "S")
        CP2e <- which(Events[[ii]]$P2[, 2] == 
            "E")
        
        if (length(CP2s) > 0 | length(CP2e) > 
            0) {
            CC <- c(CP2s, CP2e)
            Events[[ii]]$P2 <- Events[[ii]]$P2[-CC, 
                ]
        }
        
        PS2 <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$P2[, 1]))
        PE2 <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$P2[, 2]))
        Path2 <- as.matrix(cbind(PS2, PE2))
        Path2 <- Path2[order(Path2[, 1], 
            Path2[, 2]), , drop = FALSE]
        
        CPRs <- which(Events[[ii]]$Ref[, 
            1] == "S")
        CPRe <- which(Events[[ii]]$Ref[, 
            2] == "E")
        
        if (length(CPRs) > 0 | length(CPRe) > 
            0) {
            CC <- c(CPRs, CPRe)
            Events[[ii]]$Ref <- Events[[ii]]$Ref[-CC, 
                ]
        }
        
        PSR <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$Ref[, 1]))
        PER <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$Ref[, 2]))
        PathR <- as.matrix(cbind(PSR, PER))
        PathR <- PathR[order(PathR[, 1], 
            PathR[, 2]), , drop = FALSE]
        
        
        Path1 <- paste(Path1[, 1], "-", Path1[, 
            2], sep = "", collapse = ",")
        Path2 <- paste(Path2[, 1], "-", Path2[, 
            2], sep = "", collapse = ",")
        PathR <- paste(PathR[, 1], "-", PathR[, 
            2], sep = "", collapse = ",")
        
        NEv <- data.frame(EventID, GeneName, 
            EventNumber, EventType, GPos, 
            Path1, Path2, PathR, stringsAsFactors = FALSE)
        Result[[ii]] <- NEv
    }
    
    Result <- do.call(rbind, Result)
    colnames(Result) <- c("EventID", "Gene", 
        "Event Number", "Event Type", "Genomic Position", 
        "Path 1", "Path 2", "Path Reference")
    
    return(Result)
}

#' @rdname InternalFunctions
AnnotateEvents_RNASeq_MultiPath <- function(Events, 
    paths) {
    Positions <- NULL
    Result <- vector("list", length = length(Events))
    for (ii in seq_along(Events)) {
        GeneName <- as.vector(Events[[ii]]$GeneName)
        GeneID <- as.vector(Events[[ii]]$Gene)
        EventNumber <- ii
        EventID <- paste(GeneID, "_", EventNumber, 
            sep = "")
        EventType <- Events[[ii]]$Type
        Chrom <- as.vector(Events[[ii]]$P1[1, 
            "Chr"])
        
        EventNumP <- Events[[ii]]$NumP
        command <- "Positions<-rbind(Events[[ii]]$P1,"
        for (kk in 2:EventNumP) {
            if (kk == EventNumP) {
                command <- paste0(command, 
                  "Events[[ii]]$P", kk, ")[,4:5]")
            } else {
                command <- paste0(command, 
                  "Events[[ii]]$P", kk, ",")
            }
        }
        eval(parse(text = command))
        # Positions<-rbind(Events[[ii]]$P1,Events[[ii]]$P2)[,4:5]
        Start <- as.numeric(Positions[, 1])
        End <- as.numeric(Positions[, 2])
        Start <- Start[which(Start != 0)]
        End <- End[which(End != 0)]
        
        # browser()
        minGPos <- min(Start)
        maxGPos <- max(End)
        GPos <- paste(Chrom, ":", minGPos, 
            "-", maxGPos, sep = "")
        
        # CP1s<-which(Events[[ii]]$P1[,1]=='S')
        # CP1e<-which(Events[[ii]]$P1[,2]=='E')
        # if(length(CP1s)>0|length(CP1e)>0) {
        # CC<-c(CP1s,CP1e)
        # Events[[ii]]$P1<-Events[[ii]]$P1[-CC,]
        # }
        # PS1<-as.numeric(gsub('.[ab]','',Events[[ii]]$P1[,1]))
        # PE1<-as.numeric(gsub('.[ab]','',Events[[ii]]$P1[,2]))
        # Path1<-as.matrix(cbind(PS1,PE1))
        # Path1<-Path1[order(Path1[,1],Path1[,2]),,drop=FALSE]
        
        for (kk in seq_len(EventNumP)) {
            command <- paste0("CP", kk, "s<-which(Events[[ii]]$P", 
                kk, "[,1]=='S')")
            eval(parse(text = command))
            command <- paste0("CP", kk, "e<-which(Events[[ii]]$P", 
                kk, "[,2]=='E')")
            eval(parse(text = command))
            a <- paste0("a<-length(CP", kk, 
                "s)>0|length(CP", kk, "e)>0")
            eval(parse(text = a))
            if (a) {
                command <- paste0("CC<-c(CP", 
                  kk, "s,CP", kk, "e)")
                eval(parse(text = command))
                command <- paste0("Events[[ii]]$P", 
                  kk, "<-Events[[ii]]$P", 
                  kk, "[-CC,]")
                eval(parse(text = command))
            }
            command <- paste0("PS", kk, "<-as.numeric(gsub('.[ab]','',Events[[ii]]$P", 
                kk, "[,1]))")
            eval(parse(text = command))
            command <- paste0("PE", kk, "<-as.numeric(gsub('.[ab]','',Events[[ii]]$P", 
                kk, "[,2]))")
            eval(parse(text = command))
            command <- paste0("Path", kk, 
                "<-as.matrix(cbind(PS", kk, 
                ",PE", kk, "))")
            eval(parse(text = command))
            command <- paste0("Path", kk, 
                "<-Path", kk, "[order(Path", 
                kk, "[,1],Path", kk, "[,2]),,drop=FALSE]")
            eval(parse(text = command))
        }
        
        
        
        
        
        
        # CP2s<-which(Events[[ii]]$P2[,1]=='S')
        # CP2e<-which(Events[[ii]]$P2[,2]=='E')
        # if(length(CP2s)>0|length(CP2e)>0) {
        # CC<-c(CP2s,CP2e)
        # Events[[ii]]$P2<-Events[[ii]]$P2[-CC,]
        # }
        # PS2<-as.numeric(gsub('.[ab]','',Events[[ii]]$P2[,1]))
        # PE2<-as.numeric(gsub('.[ab]','',Events[[ii]]$P2[,2]))
        # Path2<-as.matrix(cbind(PS2,PE2))
        # Path2<-Path2[order(Path2[,1],Path2[,2]),,drop=FALSE]
        
        CPRs <- which(Events[[ii]]$Ref[, 
            1] == "S")
        CPRe <- which(Events[[ii]]$Ref[, 
            2] == "E")
        
        if (length(CPRs) > 0 | length(CPRe) > 
            0) {
            CC <- c(CPRs, CPRe)
            Events[[ii]]$Ref <- Events[[ii]]$Ref[-CC, 
                ]
        }
        
        PSR <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$Ref[, 1]))
        PER <- as.numeric(gsub(".[ab]", "", 
            Events[[ii]]$Ref[, 2]))
        PathR <- as.matrix(cbind(PSR, PER))
        PathR <- PathR[order(PathR[, 1], 
            PathR[, 2]), , drop = FALSE]
        
        
        # Path1<-paste(Path1[,1],'-',Path1[,2],sep='',collapse=',')
        # Path2<-paste(Path2[,1],'-',Path2[,2],sep='',collapse=',')
        NEv <- "NEv<-data.frame(EventID,GeneName,EventNumber,EventType,GPos,EventNumP,"
        for (kk in seq_len(paths)) {
            if (kk <= EventNumP) {
                command <- paste0("Path", 
                  kk, "<-paste(Path", kk, 
                  "[,1],'-',Path", kk, "[,2],sep='',collapse=',')")
                eval(parse(text = command))
                # command <-
                # paste0('ProbesP',kk,'<-paste(Events[[ii]]$Probes_P',
                # kk,',collapse=',')') eval(parse(text =
                # command))
            } else {
                command <- paste0("Path", 
                  kk, "<-'-'")
                eval(parse(text = command))
                # command <- paste0('ProbesP',kk,'<-'-'')
                # eval(parse(text = command))
            }
            NEv <- paste0(NEv, "Path", kk, 
                ",")
        }
        
        PathR <- paste(PathR[, 1], "-", PathR[, 
            2], sep = "", collapse = ",")
        NEv <- paste0(NEv, "PathR,stringsAsFactors = FALSE)")
        eval(parse(text = NEv))
        rownames(NEv) <- NULL
        # NEv<-data.frame(EventID,GeneName,EventNumber,EventType,GPos,
        # Path1,Path2,PathR,stringsAsFactors =
        # FALSE)
        Result[[ii]] <- NEv
    }
    
    Result <- do.call(rbind, Result)
    command <- "colnames(Result)<-c('EventID','Gene','Event Number','Event Type','Genomic Position','Num of Paths','Path 1',"
    
    for (kk in 2:(paths + 1)) {
        if (kk == (paths + 1)) {
            command <- paste0(command, "'Path Ref')")
        } else {
            command <- paste0(command, "'Path ", 
                kk, "',")
        }
    }
    eval(parse(text = command))
    return(Result)
}

#' @rdname InternalFunctions
AnnotateEvents_KLL <- function(Events, Gxx, 
    GenI) {
    {
        # Gxx <- GeneName
        GeneName <- Gxx
        GeneID <- GenI
        Chrom <- gsub("chr", "", as.vector(Events[[1]]$P1[1, 
            "Chr"]))
        Result <- vector("list")
        # Flat<-vector('list')
        mm <- 0
        
        for (ii in seq_along(Events)) {
            if (!any(c(identical(unique(Events[[ii]]$P1$Type), 
                "V"), identical(unique(Events[[ii]]$P2$Type), 
                "V"), identical(unique(Events[[ii]]$Ref$Type), 
                "V")) == TRUE)) {
                mm <- mm + 1
                
                EventNumber <- ii
                
                EventType <- Events[[ii]]$Type
                
                Positions <- rbind(Events[[ii]]$P1, 
                  Events[[ii]]$P2)[, 4:5]
                Start <- as.numeric(Positions[, 
                  1])
                End <- as.numeric(Positions[, 
                  2])
                Start <- Start[which(Start != 
                  0)]
                End <- End[which(End != 0)]
                
                # browser()
                minGPos <- min(Start)
                maxGPos <- max(End)
                GPos <- paste(Chrom, ":", 
                  minGPos, "-", maxGPos, 
                  sep = "")
                
                CP1s <- which(Events[[ii]]$P1[, 
                  1] == "S")
                CP1e <- which(Events[[ii]]$P1[, 
                  2] == "E")
                
                if (length(CP1s) > 0 | length(CP1e) > 
                  0) {
                  CC <- c(CP1s, CP1e)
                  Events[[ii]]$P1 <- Events[[ii]]$P1[-CC, 
                    ]
                }
                
                PS1 <- as.numeric(gsub(".[ab]", 
                  "", Events[[ii]]$P1[, 1]))
                PE1 <- as.numeric(gsub(".[ab]", 
                  "", Events[[ii]]$P1[, 2]))
                Path1 <- as.matrix(cbind(PS1, 
                  PE1))
                Path1 <- Path1[order(Path1[, 
                  1], Path1[, 2]), , drop = FALSE]
                
                CP2s <- which(Events[[ii]]$P2[, 
                  1] == "S")
                CP2e <- which(Events[[ii]]$P2[, 
                  2] == "E")
                
                if (length(CP2s) > 0 | length(CP2e) > 
                  0) {
                  CC <- c(CP2s, CP2e)
                  Events[[ii]]$P2 <- Events[[ii]]$P2[-CC, 
                    ]
                }
                
                PS2 <- as.numeric(gsub(".[ab]", 
                  "", Events[[ii]]$P2[, 1]))
                PE2 <- as.numeric(gsub(".[ab]", 
                  "", Events[[ii]]$P2[, 2]))
                Path2 <- as.matrix(cbind(PS2, 
                  PE2))
                Path2 <- Path2[order(Path2[, 
                  1], Path2[, 2]), , drop = FALSE]
                
                CPRs <- which(Events[[ii]]$Ref[, 
                  1] == "S")
                CPRe <- which(Events[[ii]]$Ref[, 
                  2] == "E")
                
                if (length(CPRs) > 0 | length(CPRe) > 
                  0) {
                  CC <- c(CPRs, CPRe)
                  Events[[ii]]$Ref <- Events[[ii]]$Ref[-CC, 
                    ]
                }
                
                PSR <- as.numeric(gsub(".[ab]", 
                  "", Events[[ii]]$Ref[, 
                    1]))
                PER <- as.numeric(gsub(".[ab]", 
                  "", Events[[ii]]$Ref[, 
                    2]))
                PathR <- as.matrix(cbind(PSR, 
                  PER))
                PathR <- PathR[order(PathR[, 
                  1], PathR[, 2]), , drop = FALSE]
                
                
                Path1 <- paste(Path1[, 1], 
                  "-", Path1[, 2], sep = "", 
                  collapse = ",")
                Path2 <- paste(Path2[, 1], 
                  "-", Path2[, 2], sep = "", 
                  collapse = ",")
                PathR <- paste(PathR[, 1], 
                  "-", PathR[, 2], sep = "", 
                  collapse = ",")
                
                
                NEv <- data.frame(GeneName, 
                  GeneID, EventNumber, EventType, 
                  GPos, Path1, Path2, PathR, 
                  stringsAsFactors = FALSE)
                Result[[mm]] <- NEv
            }
        }
        
        Result <- do.call(rbind, Result)
        
        return(Result)
    }
}


#' @rdname InternalFunctions
ClassifyEvents <- function(SG, Events, twopaths) {
    Events <- lapply(seq_along(Events), function(XX) {
        if (XX %in% twopaths) {
            # Keep the components of Path 1 and Path
            # 2
            P1 <- Events[[XX]]$P1[, seq_len(2)]
            P2 <- Events[[XX]]$P2[, seq_len(2)]
            Info <- rbind(P1, P2)
            
            # If there is an edge that leaves the
            # Start node, we have an Alternative
            # First Exon
            if (any(Info[, 1] == "S")) {
                Events[[XX]]$Type <- "Alternative First Exon"
                # next
                return(Events[[XX]])
            }
            
            # If there is an edge that enters the End
            # node, we have an Alternative Last Exon
            
            if (any(Info[, 2] == "E")) {
                Events[[XX]]$Type <- "Alternative Last Exon"
                # next
                return(Events[[XX]])
            }
            
            # Create a Mini Adjacency Graph using
            # only the elements from Path 1 and Path
            # 2
            
            # Find the From Nodes and To Nodes in the
            # complete Adjacency Matrix
            ii <- sort(match(Info[, 1], rownames(SG$Adjacency)))
            jj <- sort(match(Info[, 2], rownames(SG$Adjacency)))
            
            # The new one will have those rows and
            # columns
            MiniSGA <- SG$Adjacency[ii, jj]
            
            # Get unique values, as the elements
            # might be repeated
            iixx <- unique(rownames(MiniSGA))
            jjxx <- unique(colnames(MiniSGA))
            MiniSGA <- MiniSGA[iixx, jjxx, 
                drop = FALSE]
            
            
            # If MiniSGA has a dimension of 3x3, the
            # event could be: Cassette, Retained
            
            
            if (dim(MiniSGA)[1] == 3 & dim(MiniSGA)[2] == 
                3) {
                
                # Get the nodes in order: of path 1 (path
                # 1 is the largest and the important one)
                MisExons <- data.frame(From = Events[[XX]]$P1$From, 
                  To = Events[[XX]]$P1$To)
                Info2 <- merge(MisExons, 
                  SG$Edges)
                Info2$dist <- as.numeric(as.vector(Info2$End)) - 
                  as.numeric(as.vector(Info2$Start))
                
                MisExons <- unique(as.numeric(gsub("[.ab]", 
                  "", c(as.vector(MisExons$From), 
                    as.vector(MisExons$To)))))
                
                misDiff <- diff(sort(MisExons))
                # if the misDiff == 1 1 -> we have
                # contiguous nodes.  This is required for
                # group 1 events.
                
                AA <- Info2$dist[Info2$Type == 
                  "J"]
                l_A <- length(AA)
                # The different among the group one
                # relies on the length of the edges.
                
                if (l_A == 2 & max(misDiff) == 
                  1) {
                  if (AA[1] > 1 & AA[2] > 
                    1) {
                    # skip a exon.
                    
                    Events[[XX]]$Type <- "Cassette Exon"
                    return(Events[[XX]])
                  }
                  
                  if (AA[1] > 1 & AA[2] == 
                    1) {
                    # altern 3'
                    Events[[XX]]$Type <- "Alternative 3' Splice Site"
                    return(Events[[XX]])
                  }
                  
                  if (AA[1] == 1 & AA[2] > 
                    1) {
                    # altern 5'
                    Events[[XX]]$Type <- "Alternative 5' Splice Site"
                    return(Events[[XX]])
                  }
                  
                  if (AA[1] == 1 & AA[2] == 
                    1) {
                    # Retained Intron
                    Events[[XX]]$Type <- "Retained Intron"
                    return(Events[[XX]])
                  }
                }
            }
            
            # The case of Mutually Exclusive Exons
            
            
            if (dim(MiniSGA)[1] == 5 & dim(MiniSGA)[2] == 
                5) {
                if (MiniSGA[3, 3] == 0) {
                  # path 1
                  MisExons_1 <- data.frame(From = Events[[XX]]$P1$From, 
                    To = Events[[XX]]$P1$To)
                  Info1 <- merge(MisExons_1, 
                    SG$Edges)
                  Info1$dist <- as.numeric(as.vector(Info1$End)) - 
                    as.numeric(as.vector(Info1$Start))
                  AA_1 <- Info1$dist[Info1$Type == 
                    "J"]
                  l_1 <- length(AA_1)
                  
                  MisExons_1 <- unique(as.numeric(gsub("[.ab]", 
                    "", c(as.vector(MisExons_1$From), 
                      as.vector(MisExons_1$To)))))
                  
                  
                  # path 2
                  MisExons_2 <- data.frame(From = Events[[XX]]$P2$From, 
                    To = Events[[XX]]$P2$To)
                  Info2 <- merge(MisExons_2, 
                    SG$Edges)
                  Info2$dist <- as.numeric(as.vector(Info2$End)) - 
                    as.numeric(as.vector(Info2$Start))
                  AA_2 <- Info2$dist[Info2$Type == 
                    "J"]
                  l_2 <- length(AA_2)
                  
                  MisExons_2 <- unique(as.numeric(gsub("[.ab]", 
                    "", c(as.vector(MisExons_2$From), 
                      as.vector(MisExons_2$To)))))
                  
                  misDiff <- diff(sort(unique(c(MisExons_1, 
                    MisExons_2))))
                  
                  if (l_2 == 2 & l_1 == 2 & 
                    max(misDiff) == 1) {
                    if (!any(c(AA_1, AA_2) == 
                      1)) {
                      Events[[XX]]$Type <- "Mutually Exclusive Exons"
                      return(Events[[XX]])
                    }
                  }
                }
            }
            
            if (is.null(Events[[XX]]$Type)) {
                Events[[XX]]$Type <- "Complex Event"
                return(Events[[XX]])
            }
        } else {
            Events[[XX]]$Type <- "Multipath"
            return(Events[[XX]])
        }
        
        
        return(Events[[XX]])
    })
    
    
    
    
    if (SG$Edges[1, "Strand"] == "-") {
        Types <- sapply(seq_along(Events), 
            function(x) {
                Events[[x]]$Type
            })
        
        if (any(Types == "Alternative 5' Splice Site" | 
            Types == "Alternative 3' Splice Site" | 
            Types == "Alternative First Exon" | 
            Types == "Alternative Last Exon")) {
            Events <- lapply(seq_along(Events), 
                function(XX) {
                  if (Events[[XX]]$Type == 
                    "Alternative 5' Splice Site") {
                    Events[[XX]]$Type <- "Alternative 3' Splice Site"
                    return(Events[[XX]])
                  }
                  if (Events[[XX]]$Type == 
                    "Alternative 3' Splice Site") {
                    Events[[XX]]$Type <- "Alternative 5' Splice Site"
                    return(Events[[XX]])
                  }
                  if (Events[[XX]]$Type == 
                    "Alternative First Exon") {
                    Events[[XX]]$Type <- "Alternative Last Exon"
                    return(Events[[XX]])
                  }
                  if (Events[[XX]]$Type == 
                    "Alternative Last Exon") {
                    Events[[XX]]$Type <- "Alternative First Exon"
                    return(Events[[XX]])
                  }
                  
                  return(Events[[XX]])
                })
        }
    }
    
    
    
    return(Events)
}


#' @rdname InternalFunctions

##### Given the signals of the three paths,
##### estimate the concentrations of each of
##### the isoforms


# lambda parameter included to regularize
# the affinities

estimateAbsoluteConc <- function(Signal1, 
    Signal2, SignalR, lambda) {
    # require(nnls)
    Signal1 <- as.numeric(Signal1)
    Signal2 <- as.numeric(Signal2)
    SignalR <- as.numeric(SignalR)
    
    cols <- length(Signal1)
    
    A <- cbind(Signal1, Signal2)
    b <- SignalR
    Salida <- nnls(A, b)  # non negative least squares
    if (lambda == 0) {
        resultado <- nnls(A, b)
        Salida <- resultado$x
        u <- Salida[1]
        v <- Salida[2]
        w <- 0
        offset <- w/(1 - u - v)  # some times the offset is way too large (1-u-v = 0)
        T1est <- Signal1 * u
        T2est <- Signal2 * v
        Relerror <- as.numeric(crossprod((A[, 
            seq_len(2)]) %*% c(u, v) - b)/crossprod(b))
        
        residuals <- resultado$residuals[seq_len(cols), 
            , drop = FALSE]
        
        return(list(T1est = T1est, T2est = T2est, 
            offset = offset, Relerror = Relerror, 
            Residuals = residuals))
    }
    # Add a new equation to make the values
    # of u and v close to each other
    if (is.null(lambda)) {
        lambda <- 0.1
    }
    
    penalty <- lambda * ((Salida$deviance)/(sum((Salida$x - 
        1)^2)))
    penalty <- sqrt(penalty)
    
    A <- rbind(A, c(penalty, 0), c(0, penalty))
    b <- c(SignalR, penalty, penalty)
    
    resultado <- nnls(A, b)
    Salida <- resultado$x  # non negative least squares
    
    u <- Salida[1]
    v <- Salida[2]
    w <- 0  # No offset used in this model
    offset <- w/(1 - u - v)  # some times the offset is way too large (1-u-v = 0)
    T1est <- Signal1 * u
    T2est <- Signal2 * v
    Relerror <- as.numeric(crossprod(cbind(Signal1, 
        Signal2) %*% c(u, v) - SignalR)/crossprod(SignalR))
    # if(Relerror==0){browser()}
    residuals <- resultado$residuals[seq_len(cols), 
        , drop = FALSE]
    residuals <- residuals/SignalR
    
    return(list(T1est = T1est, T2est = T2est, 
        offset = offset, Relerror = Relerror, 
        Residuals = residuals))
}

#' @rdname InternalFunctions
estimateAbsoluteConcmultipath <- function(datos, 
    lambda = 0.1) {
    
    # require(nnls)
    l <- dim(datos)[1]
    cols <- dim(datos)[2]
    Signal <- list()
    A <- c()
    for (k in seq_len((l - 1))) {
        Signal[[k]] <- as.numeric(datos[k, 
            ])
        A <- cbind(A, Signal[[k]])
    }
    Signal[[l]] <- as.numeric(datos[l, ])
    b <- Signal[[l]]
    Salida <- nnls(A, b)
    u <- c()
    Tset <- matrix(0, ncol = cols, nrow = (l - 
        1))
    if (lambda == 0) {
        resultado <- nnls(A, b)
        Salida <- resultado$x
        for (k in seq_len((l - 1))) {
            u <- c(u, Salida[k])
            Tset[k, ] <- Signal[[k]] * u[k]
        }
        w <- 0
        offset <- w/(1 - sum(u))
        Relerror <- as.numeric(crossprod((A[, 
            seq_len((l - 1))]) %*% u - b)/crossprod(b))
        residuals <- resultado$residuals[seq_len(cols), 
            , drop = FALSE]
        return(list(Tset = Tset, offset = offset, 
            Relerror = Relerror, Residuals = residuals))
    }
    if (is.null(lambda)) {
        lambda <- 0.1
    }
    
    
    penalty <- lambda * ((Salida$deviance)/(sum((Salida$x - 
        1)^2)))
    penalty <- sqrt(penalty)
    
    
    penal <- diag(penalty, nrow = (l - 1))
    A <- rbind(A, penal)
    b <- c(b, rep(penalty, (l - 1)))
    resultado <- nnls(A, b)
    Salida <- resultado$x
    
    for (k in seq_len((l - 1))) {
        u <- c(u, Salida[k])
        Tset[k, ] <- Signal[[k]] * u[k]
    }
    w <- 0
    offset <- w/(1 - sum(u))
    
    Relerror <- as.numeric(crossprod((A[seq_len(cols), 
        seq_len((l - 1))]) %*% u - b[seq_len(cols)])/crossprod(b[seq_len(cols)]))
    
    residuals <- resultado$residuals[seq_len(cols), 
        , drop = FALSE]
    
    
    return(list(Tset = Tset, offset = offset, 
        Relerror = Relerror, Residuals = residuals))
}



#' @rdname InternalFunctions
findTriplets <- function(randSol, tol = 1e-08) {
    
    # Compute the Distance Matrix from the
    # matrix of fluxes
    X <- as.matrix(dist(randSol))
    
    # Which distances are smaller than the
    # tolerance (To ensure the flux is the
    # same always)
    Inc <- (X < tol)
    
    # Create a graph from the adjacency
    # matrix and find the connected
    # components
    g <- graph_from_adjacency_matrix(Inc)
    Groups <- clusters(g)
    
    EdgG_Flux <- randSol[match(seq_len(Groups$no), 
        Groups$membership), ]
    
    # All possible combination of two
    # elements from the graph to create all
    # the posible sums to find the triplets
    # of events
    Index <- combn(nrow(EdgG_Flux), 2)
    flowsum <- EdgG_Flux[Index[1, ], , drop = FALSE] + 
        EdgG_Flux[Index[2, ], , drop = FALSE]  # All the possible sums
    
    # Calculate the distance between all the
    # possible sums of every element of the
    # graph and the flow matrix. The Events
    # will be those in which the distance is
    # smaller than the tolerance (almost
    # equal to 0). The tolerance is used to
    # avoid ounding problems
    DistanceMat <- pdist2(flowsum, EdgG_Flux)
    
    x_Ref <- which(DistanceMat < tol, arr.ind = TRUE)
    
    # Obtain Paths
    P1 <- Index[1, x_Ref[, 1]]
    P2 <- Index[2, x_Ref[, 1]]
    Ref <- x_Ref[, 2]
    
    GG <- Groups$membership
    
    triplets <- cbind(P1, P2, Ref)
    
    return(list(groups = GG, triplets = triplets))
}

#' @rdname InternalFunctions
findTriplets2 <- function(Incidence, paths = 2, 
    randSol) {
    # randSol <- getRandomFlow(Incidence,
    # ncol = 2)
    X <- as.matrix(dist(randSol))
    tol <- 1e-08
    Inc <- (X < tol)
    g <- graph_from_adjacency_matrix(Inc)
    Groups <- clusters(g)  # This approach looks to use too heavy weapons
    # to solve the problem...
    
    if (Groups$no == 2) {
        multipaths <- matrix(NA, nrow = 0, 
            ncol = 1)
        return(list(groups = Groups$membership, 
            multipaths = multipaths))
    }
    
    # TODO: Build the triplets using only the
    # unique fluxes NewIncidence <- Incidence
    # %*% EdgeXG NewIncidence <-
    # Incidence[,apply(EdgeXG, 2,which.max)]
    
    NewIncidence <- Incidence
    csI <- colCumsums(NewIncidence)
    # rownames(csI) <- rownames(Incidence)
    
    mg <- matrix(0, nrow = ncol(csI), ncol = Groups$no[1])
    mg[cbind(seq_len(length(Groups$membership)), 
        Groups$membership)] <- 1
    colnames(mg) <- seq_len(ncol(mg))
    # for (kk in seq_len(ncol(mg))){
    # mg[,kk]<-Groups$membership==kk }
    csI <- csI %*% mg
    csI <- uniquefast(csI)
    Index <- combn(nrow(csI), 2)
    
    BigDeltacsI <- csI[Index[2, ], ] - csI[Index[1, 
        ], ]
    
    BigDeltacsI <- uniquefast(BigDeltacsI)
    
    Ones <- rowSums(BigDeltacsI == 1)
    Several <- rowSums(BigDeltacsI != 0)
    MinusOnes <- rowSums(BigDeltacsI == -1)
    
    multipaths <- matrix(ncol = paths + 2)
    colnames(multipaths) <- c(paste(rep("p", 
        paths), sep = "", c(seq_len(paths))), 
        "Ref", "NumP")
    
    for (ii in 2:paths) {
        # Severalgood <- (Several >2) & (Several
        # == ii+1) & GoodOnes
        Severalgood <- (Several > 2) & (Several == 
            ii + 1)
        # Type One
        TypeOne <- which((Ones == 1) & Severalgood)
        if (length(TypeOne) > 0) {
            P12 <- apply(BigDeltacsI[TypeOne, 
                , drop = FALSE], 1, FUN = function(x) {
                which(x == -1)
            })
            PR <- apply(BigDeltacsI[TypeOne, 
                , drop = FALSE], 1, FUN = function(x) {
                which(x == 1)
            })
            
            comb <- cbind(t(P12), PR)
            # comb <- matrix(Groups$membership[comb],
            # ncol = ncol(comb)) comb <-
            # cbind(t(apply(comb[,seq_len(ii),drop=FALSE],1,
            # function(x){x[order(x)]})),comb[,ii+1])
            # comb <- unique(comb)
            
            A <- matrix(0, nrow = dim(comb)[1], 
                ncol = paths + 2)
            A[, seq_len(ii)] <- comb[, seq_len(ii)]
            A[, paths + 1] <- comb[, ii + 
                1]
            A[, paths + 2] <- ii
            multipaths <- rbind(multipaths, 
                A)
        }
        # Type MinusOne
        TypeMinusOne <- which((MinusOnes == 
            1) & Severalgood)
        if (length(TypeMinusOne) > 0) {
            P12 <- apply(BigDeltacsI[TypeMinusOne, 
                , drop = FALSE], 1, FUN = function(x) {
                which(x == 1)
            })
            PR <- apply(BigDeltacsI[TypeMinusOne, 
                , drop = FALSE], 1, FUN = function(x) {
                which(x == -1)
            })
            
            comb <- cbind(t(P12), PR)
            # comb <- matrix(Groups$membership[comb],
            # ncol = ncol(comb)) comb <-
            # cbind(t(apply(comb[,1:ii,drop=FALSE],1,
            # function(x){x[order(x)]})),comb[,ii+1])
            # comb <- unique(comb)
            
            A <- matrix(0, nrow = dim(comb)[1], 
                ncol = paths + 2)
            A[, seq_len(ii)] <- comb[, seq_len(ii)]
            A[, paths + 1] <- comb[, ii + 
                1]
            A[, paths + 2] <- ii
            multipaths <- rbind(multipaths, 
                A)
        }
    }
    multipaths <- unique(multipaths)
    multipaths <- multipaths[-1, ]
    if (is.null(nrow(multipaths))) {
        multipaths <- t(multipaths)
    }
    
    return(list(groups = Groups$membership, 
        multipaths = multipaths))
}


#' @rdname InternalFunctions
GetCounts <- function(Events, sg_txiki, type = "counts") {
    readsC <- counts(sg_txiki)
    readsF <- FPKM(sg_txiki)
    countsEvents <- lapply(Events, getPathCounts, 
        readsC)
    countsEvents <- lapply(countsEvents, 
        getPathFPKMs, readsF)
    return(countsEvents)
}

getPathCounts <- function(x, readsC, widthinit) {
    reads <- rbind(colSums(readsC[x$P1$featureID, 
        , drop = FALSE]), colSums(readsC[x$P2$featureID, 
        , drop = FALSE]), colSums(readsC[x$Ref$featureID, 
        , drop = FALSE]))
    
    rownames(reads) <- c("P1", "P2", "Ref")
    x$Counts <- reads
    return(x)
}

getPathFPKMs <- function(x, readsC, widthinit) {
    reads <- rbind(colSums(readsC[x$P1$featureID, 
        , drop = FALSE]), colSums(readsC[x$P2$featureID, 
        , drop = FALSE]), colSums(readsC[x$Ref$featureID, 
        , drop = FALSE]))
    
    rownames(reads) <- c("P1", "P2", "Ref")
    x$FPKM <- reads
    return(x)
}

#' @rdname InternalFunctions
GetCountsMP <- function(Events, sg_txiki, 
    type = "counts") {
    readsC <- counts(sg_txiki)
    readsF <- FPKM(sg_txiki)
    countsEvents <- lapply(Events, getPathCountsMP, 
        readsC)
    countsEvents <- lapply(countsEvents, 
        getPathFPKMsMP, readsF)
    return(countsEvents)
}

getPathCountsMP <- function(x, readsC, widthinit) {
    command <- "reads <- rbind(colSums(readsC[x$P1$featureID,,drop = FALSE]),"
    for (i in 2:(x$NumP + 1)) {
        if (i == (x$NumP + 1)) {
            command <- paste0(command, "colSums(readsC[x$Ref$featureID,,drop = FALSE]))")
        } else {
            command <- paste0(command, "colSums(readsC[x$P", 
                i, "$featureID,,drop = FALSE]),")
        }
    }
    eval(parse(text = command))
    # reads <-
    # rbind(colSums(readsC[x$P1$featureID,,drop
    # = FALSE]),
    # colSums(readsC[x$P2$featureID,,drop =
    # FALSE]),
    # colSums(readsC[x$Ref$featureID,,drop =
    # FALSE]))
    
    rownames(reads) <- c(paste("P", seq_len(x$NumP), 
        sep = ""), "Ref")
    x$Counts <- reads
    return(x)
}

getPathFPKMsMP <- function(x, readsC, widthinit) {
    command <- "reads <- rbind(colSums(readsC[x$P1$featureID,,drop = FALSE]),"
    for (i in 2:(x$NumP + 1)) {
        if (i == (x$NumP + 1)) {
            command <- paste0(command, "colSums(readsC[x$Ref$featureID,,drop = FALSE]))")
        } else {
            command <- paste0(command, "colSums(readsC[x$P", 
                i, "$featureID,,drop = FALSE]),")
        }
    }
    eval(parse(text = command))
    # reads <-
    # rbind(colSums(readsC[x$P1$featureID,,drop
    # = FALSE]),
    # colSums(readsC[x$P2$featureID,,drop =
    # FALSE]),
    # colSums(readsC[x$Ref$featureID,,drop =
    # FALSE]))
    
    rownames(reads) <- c(paste("P", seq_len(x$NumP), 
        sep = ""), "Ref")
    x$FPKM <- reads
    return(x)
}




#' @rdname InternalFunctions

getEventPaths <- function(Events, SG) {
    Groups <- Events$groups
    Triplets <- Events$triplets
    
    P1 <- lapply(seq_len(nrow(Triplets)), 
        function(x) {
            A <- SG$Edges[which(Groups == 
                Triplets[x, 1]), ]
            return(A)
        })
    P2 <- lapply(seq_len(nrow(Triplets)), 
        function(x) {
            A <- SG$Edges[which(Groups == 
                Triplets[x, 2]), ]
            return(A)
        })
    Ref <- lapply(seq_len(nrow(Triplets)), 
        function(x) {
            A <- SG$Edges[which(Groups == 
                Triplets[x, 3]), ]
            return(A)
        })
    
    Result <- lapply(seq_along(P1), function(X) {
        if (nrow(P1[[X]]) > nrow(P2[[X]])) {
            A <- list(P1 = P1[[X]], P2 = P2[[X]], 
                Ref = Ref[[X]])
        } else {
            A <- list(P1 = P2[[X]], P2 = P1[[X]], 
                Ref = Ref[[X]])
        }
        
        return(A)
    })
    
    return(Result)
}

#' @rdname InternalFunctions
getEventMultiPaths <- function(Events, SG, 
    twopaths, paths) {
    P1 <- NULL
    A <- NULL
    multipaths <- Events$multipaths
    Groups <- Events$groups
    
    
    for (ii in seq_len(paths)) {
        command <- paste(paste("P", ii, sep = ""), 
            "<-lapply(seq_len(nrow(multipaths)),function(x){A<-SG$Edges[which(Groups==multipaths[x,ii]),];return(A)})", 
            sep = "")
        eval(parse(text = command))
    }
    
    Ref <- lapply(seq_len(nrow(multipaths)), 
        function(x) {
            A <- SG$Edges[which(Groups == 
                multipaths[x, paths + 1]), 
                ]
            return(A)
        })
    NumP <- lapply(seq_len(nrow(multipaths)), 
        function(x) {
            A <- multipaths[x, paths + 2]
            return(A)
        })
    X <- 1
    Result <- lapply(seq_along(P1), function(X) {
        command <- "A <- list("
        for (i in seq_len((paths + 2))) {
            if (i == 1) {
                command <- paste(command, 
                  paste(paste(paste(paste("P", 
                    i, sep = ""), "=P", sep = ""), 
                    i, sep = ""), "[[X]]", 
                    sep = ""), sep = "")
            } else if (i == (paths + 1)) {
                command <- paste(command, 
                  "Ref=Ref[[X]]", sep = ",")
            } else if (i == (paths + 2)) {
                command <- paste(command, 
                  "NumP=NumP[[X]])", sep = ",")
            } else {
                command <- paste(command, 
                  paste(paste(paste(paste("P", 
                    i, sep = ""), "=P", sep = ""), 
                    i, sep = ""), "[[X]]", 
                    sep = ""), sep = ",")
            }
        }
        eval(parse(text = command))
        return(A)
    })
    
    if (length(twopaths) > 0) {
        for (j in seq_len(length(twopaths))) {
            if (nrow(Result[[twopaths[j]]]$P2) > 
                nrow(Result[[twopaths[j]]]$P1)) {
                d <- Result[[twopaths[j]]]$P2
                Result[[twopaths[j]]]$P2 <- Result[[twopaths[j]]]$P1
                Result[[twopaths[j]]]$P1 <- d
            }
        }
    }
    
    return(Result)
}

#' @rdname InternalFunctions
GetIGVPaths <- function(EventInfo, SG_Edges) {
    
    Gene <- as.vector(EventInfo[1, 1])
    
    Path1 <- matrix(unlist(strsplit(as.vector(EventInfo[, 
        "Path.1"]), "[,-]")), ncol = 2, byrow = TRUE)
    Path2 <- matrix(unlist(strsplit(as.vector(EventInfo[, 
        "Path.2"]), "[,-]")), ncol = 2, byrow = TRUE)
    PathR <- matrix(unlist(strsplit(as.vector(EventInfo[, 
        "Path.Reference"]), "[,-]")), ncol = 2, 
        byrow = TRUE)
    
    SG_Edges_Orig <- SG_Edges
    
    SG_Edges[, 1] <- gsub(".[ab]", "", SG_Edges[, 
        1])
    SG_Edges[, 2] <- gsub(".[ab]", "", SG_Edges[, 
        2])
    
    P1.ix <- apply(Path1, 1, function(x) {
        ix <- which(SG_Edges[, 1] == x[1] & 
            SG_Edges[, 2] == x[2])
        return(ix)
    })
    P2.ix <- apply(Path2, 1, function(x) {
        ix <- which(SG_Edges[, 1] == x[1] & 
            SG_Edges[, 2] == x[2])
        return(ix)
    })
    PR.ix <- apply(PathR, 1, function(x) {
        ix <- which(SG_Edges[, 1] == x[1] & 
            SG_Edges[, 2] == x[2])
        return(ix)
    })
    
    Path1 <- SG_Edges[P1.ix, ]
    Path2 <- SG_Edges[P2.ix, ]
    PathR <- SG_Edges[PR.ix, ]
    
    
    PlotPath1 <- c()
    
    for (ii in seq_len(nrow(Path1))) {
        Type <- as.vector(Path1[ii, "Type"])
        
        if (Type == "J") {
            Chr <- as.vector(rep(Path1[ii, 
                "Chr"], 2))
            St <- as.numeric(c(as.vector(Path1[ii, 
                "Start"]), as.vector(Path1[ii, 
                "End"])))
            Ed <- St
            Wd <- as.numeric(rep(0, 2))
            Str <- as.vector(rep(Path1[ii, 
                "Strand"], 2))
            Gn <- as.vector(rep(Gene, 2))
            Trs <- rep("A", 2)
            Res <- data.frame(chromosome = Chr, 
                start = St, end = Ed, width = Wd, 
                strand = Str, gene = Gn, 
                transcript = Trs, stringsAsFactors = FALSE)
            PlotPath1 <- rbind(PlotPath1, 
                Res)
        } else if (Type == "E") {
            Chr <- as.vector(Path1[ii, "Chr"])
            St <- as.numeric(as.vector(Path1[ii, 
                "Start"]))
            Ed <- as.numeric(as.vector(Path1[ii, 
                "End"]))
            Wd <- Ed - St
            Str <- as.vector(Path1[ii, "Strand"])
            Res <- data.frame(chromosome = Chr, 
                start = St, end = Ed, width = Wd, 
                strand = Str, gene = Gene, 
                transcript = "A", stringsAsFactors = FALSE)
            PlotPath1 <- rbind(PlotPath1, 
                Res)
        }
    }
    
    
    PlotPath2 <- c()
    
    for (ii in seq_len(nrow(Path2))) {
        Type <- as.vector(Path2[ii, "Type"])
        
        if (Type == "J") {
            Chr <- as.vector(rep(Path2[ii, 
                "Chr"], 2))
            St <- as.numeric(c(as.vector(Path2[ii, 
                "Start"]), as.vector(Path2[ii, 
                "End"])))
            Ed <- St
            Wd <- as.numeric(rep(0, 2))
            Str <- as.vector(rep(Path2[ii, 
                "Strand"], 2))
            Gn <- as.vector(rep(Gene, 2))
            Trs <- rep("B", 2)
            Res <- data.frame(chromosome = Chr, 
                start = St, end = Ed, width = Wd, 
                strand = Str, gene = Gn, 
                transcript = Trs, stringsAsFactors = FALSE)
            PlotPath2 <- rbind(PlotPath2, 
                Res)
        } else if (Type == "E") {
            Chr <- as.vector(Path2[ii, "Chr"])
            St <- as.numeric(as.vector(Path2[ii, 
                "Start"]))
            Ed <- as.numeric(as.vector(Path2[ii, 
                "End"]))
            Wd <- Ed - St
            Str <- as.vector(Path2[ii, "Strand"])
            Res <- data.frame(chromosome = Chr, 
                start = St, end = Ed, width = Wd, 
                strand = Str, gene = Gene, 
                transcript = "B", stringsAsFactors = FALSE)
            PlotPath2 <- rbind(PlotPath2, 
                Res)
        }
    }
    
    Ref.Group <- connectedComp(ftM2graphNEL(as.matrix(SG_Edges_Orig[rownames(PathR), 
        seq_len(2)])))
    
    
    for (ii in seq_along(Ref.Group)) {
        LL <- length(Ref.Group[[ii]])
        Ref.Group[[ii]] <- cbind(Ref.Group[[ii]][seq_len((LL - 
            1))], Ref.Group[[ii]][2:LL])
        ixx <- row.match(as.data.frame(Ref.Group[[ii]]), 
            SG_Edges_Orig[, c(1, 2)])
        RR <- SG_Edges_Orig[ixx, ]
        RR[, 1] <- gsub(".[ab]", "", RR[, 
            1])
        RR[, 2] <- gsub(".[ab]", "", RR[, 
            2])
        Trs <- rep(paste("Ref", ii, sep = ""), 
            nrow(RR))
        RR <- cbind(RR, Trs)
        Ref.Group[[ii]] <- RR
    }
    
    Reference <- do.call(rbind, Ref.Group)
    PlotReference <- c()
    
    for (ii in seq_len(nrow(Reference))) {
        Type <- as.vector(Reference[ii, "Type"])
        
        if (Type == "J") {
            Chr <- as.vector(rep(Reference[ii, 
                "Chr"], 2))
            St <- as.numeric(c(as.vector(Reference[ii, 
                "Start"]), as.vector(Reference[ii, 
                "End"])))
            Ed <- St
            Wd <- as.numeric(rep(0, 2))
            Str <- as.vector(rep(Reference[ii, 
                "Strand"], 2))
            Gn <- as.vector(rep(Gene, 2))
            Trs <- rep(Reference[ii, "Trs"], 
                2)
            Res <- data.frame(chromosome = Chr, 
                start = St, end = Ed, width = Wd, 
                strand = Str, gene = Gn, 
                transcript = Trs, stringsAsFactors = FALSE)
            PlotReference <- rbind(PlotReference, 
                Res)
        } else if (Type == "E") {
            Chr <- as.vector(Reference[ii, 
                "Chr"])
            St <- as.numeric(as.vector(Reference[ii, 
                "Start"]))
            Ed <- as.numeric(as.vector(Reference[ii, 
                "End"]))
            Wd <- Ed - St
            Str <- as.vector(Reference[ii, 
                "Strand"])
            Res <- data.frame(chromosome = Chr, 
                start = St, end = Ed, width = Wd, 
                strand = Str, gene = Gene, 
                transcript = Reference[ii, 
                  "Trs"], stringsAsFactors = FALSE)
            PlotReference <- rbind(PlotReference, 
                Res)
        }
    }
    
    Plot <- rbind(PlotPath1, PlotPath2, PlotReference)
    # Plot[,1]<-paste('chr',Plot[,1],sep='')
    
    return(Plot)
}

#' @rdname InternalFunctions
getPSI <- function(ExFit, lambda = 0.1) {
    # Create matrix to fill with PSI values
    # (1 per event and sample)
    PSI <- matrix(0, nrow = nrow(ExFit)/3, 
        ncol = ncol(ExFit) - 5)
    colnames(PSI) <- colnames(ExFit[6:ncol(ExFit)])
    rownames(PSI) <- ExFit[seq(1, nrow(ExFit), 
        by = 3), 1]
    
    NCols <- ncol(ExFit)
    ExFit2 <- as.matrix(ExFit[, 6:NCols])
    Residuals <- PSI
    # Perform the operations for every
    # detectable alternative splicing event
    for (n in seq_len(nrow(ExFit)/3)) {
        
        # Get expression signal from path 1
        Signal1 <- ExFit2[1 + 3 * (n - 1), 
            ]
        
        # Get expression signal from path 2
        Signal2 <- ExFit2[2 + 3 * (n - 1), 
            ]
        
        # Get expression signal from Reference
        SignalR <- ExFit2[3 + 3 * (n - 1), 
            ]
        
        # Function to estimate concentrations
        # from the interrogated isoforms
        Output <- estimateAbsoluteConc(Signal1, 
            Signal2, SignalR, lambda)
        
        # Compute the actual PSI value (T1/T1+T2)
        psi <- Output$T1est/(Output$T1est + 
            Output$T2est)
        PSI[n, ] <- psi
        Residuals[n, ] <- Output$Residuals
    }
    return(list(PSI = PSI, Residuals = Residuals))
}

getPSImultipath <- function(ExFit, lambda = 0.1) {
    
    # EventNames <- rownames(ExFit)
    EventNames <- ExFit$unitName
    
    EventNames <- sapply(strsplit(EventNames, 
        "_"), function(x) {
        a <- paste(x[1], x[2], sep = "")
        return(a)
    })
    EventNames <- as.matrix(table(EventNames))
    numrows <- sum(EventNames - 1)  # for each path we calculate the PSI
    # (PSI1, PSI2,..., PSIn)
    
    PSI <- matrix(0, nrow = numrows, ncol = ncol(ExFit) - 
        5)
    colnames(PSI) <- colnames(ExFit)[6:ncol(ExFit)]
    rownames(PSI) <- rep(rownames(EventNames), 
        EventNames[seq_len(nrow(EventNames))] - 
            1)
    
    Residuals <- matrix(0, ncol = ncol(PSI), 
        nrow = length(rownames(EventNames)))
    rownames(Residuals) <- rownames(EventNames)
    
    NCols <- ncol(ExFit)
    ExFit2 <- as.matrix(ExFit[, 6:NCols])
    
    s <- rbind(1, EventNames)
    s <- cumsum(s)
    s <- s[-length(s)]
    e <- as.numeric(s + EventNames - 1)
    
    sp <- rbind(1, (EventNames - 1))
    sp <- cumsum(sp)
    sp <- sp[-length(sp)]
    ep <- as.numeric(sp + EventNames - 2)
    
    for (n in seq_len(length(e))) {
        datos <- ExFit2[s[n]:e[n], ]
        Output <- estimateAbsoluteConcmultipath(datos, 
            lambda)
        Tset <- Output$Tset
        TR <- apply(Tset, 2, sum)
        
        datospsi <- apply(Tset, 1, function(X) {
            return(X/TR)
        })
        datospsi <- t(datospsi)
        PSI[sp[n]:ep[n], ] <- datospsi
        Relerror <- Output$Relerror
        Residuals[n, ] <- Output$Residuals
    }
    result <- list(PSI = PSI, Residuals = Residuals)
    return(result)
}


#' @rdname InternalFunctions
getPSI_RNASeq <- function(Result, lambda = 0.1) {
    CountMatrix <- vector("list", length = length(Result))
    Vec <- c()
    
    for (jj in seq_along(Result)) {
        # print(jj)
        A <- Result[[jj]]
        
        if (!is.null(A)) {
            Evs_Counts <- lapply(A, function(X) {
                Res <- X$FPKM
                return(Res)
            })
            names(Evs_Counts) <- seq_len(length(Evs_Counts))
            Ids <- paste(A[[1]]$Gene, "_", 
                names(Evs_Counts), sep = "")
            Ids <- rep(Ids, each = 3)
            Vec <- c(Vec, Ids)
            Evs_Counts <- do.call(rbind, 
                Evs_Counts)
            
            if (!any(is.na(Evs_Counts))) {
                CountMatrix[[jj]] <- Evs_Counts
            } else {
                CountMatrix[[jj]] <- NULL
            }
        } else {
            
            
        }
    }
    
    
    Ids <- rep(c("_P1", "_P2", "_Ref"), length(Vec)/3)
    CountMatrix <- do.call(rbind, CountMatrix)
    rownames(CountMatrix) <- paste(Vec, Ids, 
        sep = "")
    
    PSI <- matrix(0, nrow = nrow(CountMatrix)/3, 
        ncol = ncol(CountMatrix))
    colnames(PSI) <- colnames(CountMatrix)
    rownames(PSI) <- Vec[seq(1, length(Vec), 
        by = 3)]
    
    Residuals <- PSI
    
    for (n in seq_len(nrow(CountMatrix)/3)) {
        Signal1 <- CountMatrix[1 + 3 * (n - 
            1), ]
        Signal2 <- CountMatrix[2 + 3 * (n - 
            1), ]
        SignalR <- CountMatrix[3 + 3 * (n - 
            1), ]
        Output <- estimateAbsoluteConc(Signal1, 
            Signal2, SignalR, lambda)
        psi <- Output$T1est/(Output$T1est + 
            Output$T2est)
        PSI[n, ] <- psi
        Residuals[n, ] <- Output$Residuals
    }
    
    return(list(PSI = PSI, Residuals = Residuals))
}

#' @rdname InternalFunctions
getPSI_RNASeq_MultiPath <- function(Result, 
    lambda = 0.1) {
    CountMatrix <- vector("list", length = length(Result))
    Vec <- c()
    indices <- c()
    # seq_along(Result)
    for (jj in seq_along(Result)) {
        # jj<-1 print(jj)
        A <- Result[[jj]]
        
        if (!is.null(A)) {
            Evs_Counts <- lapply(A, function(X) {
                Res <- X$FPKM
                return(Res)
            })
            nump <- sapply(A, function(X) {
                Res <- X$NumP
                return(Res)
            })
            nump <- nump + 1
            indices <- c(indices, nump)
            names(Evs_Counts) <- seq_len(length(Evs_Counts))
            Ids <- paste(A[[1]]$Gene, "_", 
                names(Evs_Counts), sep = "")
            Ids <- rep(Ids, nump)
            Vec <- c(Vec, Ids)
            Evs_Counts <- do.call(rbind, 
                Evs_Counts)
            
            if (!any(is.na(Evs_Counts))) {
                CountMatrix[[jj]] <- Evs_Counts
            } else {
                CountMatrix[[jj]] <- NULL
            }
        } else {
            
        }
    }
    
    # Ids<-rep(c('_P1','_P2','_Ref'),length(Vec)/3)
    Ids <- vector(mode = "character", length = length(Vec))
    # indices<-table(Vec)
    maxp <- max(indices)
    EventNames <- indices
    Ids[cumsum(indices)] <- "_Ref"
    indices <- c(1, indices)
    indices <- indices[-length(indices)]
    Ids[cumsum(indices)] <- "_P1"
    Ids[1 + cumsum(indices)] <- "_P2"
    if (maxp > 3) {
        for (kk in 2:(maxp - 2)) {
            Ids[kk + cumsum(indices)[which(Ids[kk + 
                cumsum(indices)] == "")]] <- paste0("_P", 
                kk + 1)
        }
    }
    # 
    
    CountMatrix <- do.call(rbind, CountMatrix)
    rownames(CountMatrix) <- paste(Vec, Ids, 
        sep = "")
    
    names(EventNames) <- NULL
    
    
    # PSI <- matrix(0, nrow =
    # nrow(CountMatrix)/3, ncol =
    # ncol(CountMatrix))
    numrows <- sum(EventNames - 1)  # for each path we calculate the PSI
    # (PSI1, PSI2,..., PSIn)
    PSI <- matrix(0, nrow = numrows, ncol = ncol(CountMatrix))
    
    colnames(PSI) <- colnames(CountMatrix)
    # rownames(PSI) <-
    # Vec[seq(1,length(Vec),by = 3)]
    rownames(PSI) <- rownames(CountMatrix)[-cumsum(EventNames)]
    
    namesrowres <- Vec[cumsum(EventNames)]
    Residuals <- matrix(0, nrow = length(namesrowres), 
        ncol = ncol(PSI))
    rownames(Residuals) <- namesrowres
    
    s <- c(1, EventNames)
    s <- cumsum(s)
    s <- s[-length(s)]
    e <- as.numeric(s + EventNames - 1)
    
    sp <- c(1, (EventNames - 1))
    sp <- cumsum(sp)
    sp <- sp[-length(sp)]
    ep <- as.numeric(sp + EventNames - 2)
    
    for (n in seq_len(length(e))) {
        datos <- CountMatrix[s[n]:e[n], , 
            drop = FALSE]
        Output <- estimateAbsoluteConcmultipath(datos, 
            lambda)
        Tset <- Output$Tset
        TR <- apply(Tset, 2, sum)
        
        datospsi <- apply(Tset, 1, function(X) {
            return(X/TR)
        })
        datospsi <- t(datospsi)
        PSI[sp[n]:ep[n], ] <- datospsi
        Relerror <- Output$Relerror
        
        Residuals[n, ] <- Output$Residuals
    }
    result <- list(PSI = PSI, Residuals = Residuals)
    return(result)
}



#' @rdname InternalFunctions
getRandomFlow <- function(Incidence, ncol = 1) {
    # With the incidence matrix, it is
    # possible to get its null-space and
    # generate an arbitrary flow on it. Using
    # the flow it is possible to get the
    # triplets of events.
    
    # The seed is set to ensure the order of
    # events remains the same
    set.seed("0xABBA")
    
    # Solve the Null Space for the Incidence
    # Matrix
    solh <- Null(t(Incidence))
    
    # Condition to ensure that everything
    # that exits the Start node (-1), exits
    # at the End Node (1)
    solp <- ginv(Incidence) %*% c(-1, rep(0, 
        nrow(Incidence) - 2), 1)
    
    # Matrix of fluxes, with as many columns
    # as specified by the user
    v <- matrix(runif(ncol(solh) * ncol), 
        ncol = ncol)
    randSol <- as.vector(solp) + solh %*% 
        v
    
    return(randSol)
}

#' @rdname InternalFunctions
IHsummarization <- function(Pv1, t1, Pv2, 
    t2, coherence = "Opposite") {
    if (coherence == "Equal") {
        nPv1 <- (Pv1/2) * (t1 > 0) + (1 - 
            Pv1/2) * (t1 <= 0)
        nPv2 <- (Pv2/2) * (t2 > 0) + (1 - 
            Pv2/2) * (t2 <= 0)
        Psuma <- nPv1 + nPv2
        PIH <- (Psuma^2)/2 * (Psuma < 1) + 
            (1 - (2 - Psuma)^2/2) * (Psuma >= 
                1)
        ZIH <- qnorm(PIH)
        PIH_2tail <- PIH * 2 * (PIH < 0.5) + 
            ((1 - PIH) * 2) * (PIH >= 0.5)
        
        # 1: Pv1 > 0.5 ; Pv2 > 0.5 ; t1 >0 ; t2
        # >0
        
        Cambiar <- which(Pv1 < 0.5 & Pv2 < 
            0.5 & t1 <= 0 & t2 <= 0)
        nPv1[Cambiar] <- Pv1[Cambiar]/2
        nPv2[Cambiar] <- Pv2[Cambiar]/2
        Psuma[Cambiar] <- nPv1[Cambiar] + 
            nPv2[Cambiar]
        PIH[Cambiar] <- (Psuma[Cambiar]^2)/2
        ZIH[Cambiar] <- -qnorm(PIH[Cambiar])
        PIH_2tail[Cambiar] <- (PIH[Cambiar]) * 
            2
        
        return(list(Pvalues = PIH_2tail, 
            Tstats = ZIH))
    } else {
        return(IHsummarization(Pv1, t1, Pv2, 
            -t2, coherence = "Equal"))
    }
}

#' @rdname InternalFunctions
pdist2 <- function(X, Y) {
    X1 <- rowSums(X * X)
    Y1 <- rowSums(Y * Y)
    Z <- outer(X1, Y1, "+") - 2 * X %*% t(Y)
    return(Z)
}


#' @rdname InternalFunctions
PrepareCountData <- function(Result) {
    CountMatrix <- vector("list", length = length(Result))
    
    for (jj in seq_along(Result)) {
        # print(jj)
        A <- Result[[jj]]
        
        if (!is.null(A)) {
            # Evs_Counts<-lapply(A,function(X){Res<-X$Counts;return(Res)})
            Evs_Counts <- lapply(A, function(X) {
                Res <- X$FPKM
                return(Res)
            })
            Evs_Counts <- lapply(Evs_Counts, 
                function(X) {
                  X <- X[c("Ref", "P1", "P2"), 
                    ]
                  return(X)
                })
            Cols <- ncol(Evs_Counts[[1]])
            Mat <- matrix(unlist(Evs_Counts), 
                ncol = length(A))
            colnames(Mat) <- paste(jj, "_", 
                seq_len(length(A)), sep = "")
            Samples <- colnames(Evs_Counts[[1]])
            rownames(Mat) <- paste(rep(Samples, 
                each = 3), c("_Ref", "_P1", 
                "_P2"), sep = "")
            colnames(Mat) <- paste(A[[1]]$Gene, 
                seq_len(length(A)), sep = "_")
            
            if (!any(is.na(Mat))) {
                CountMatrix[[jj]] <- Mat
            } else {
                CountMatrix[[jj]] <- NULL
            }
        } else {
            
            
        }
    }
    
    CountMatrix <- do.call(cbind, CountMatrix)
    return(CountMatrix)
}

#' @rdname InternalFunctions
PrepareProbes <- function(Probes, Class) {
    if (Class == "PSR") {
        Probes <- Probes[, c(seq_len(4), 
            8:11, 7)]
        colnames(Probes) <- c("Probe ID", 
            "X Coord", "Y Coord", "Gene", 
            "Chr", "Start", "Stop", "Strand", 
            "Probe Sequence")
    } else if (Class == "Junction") {
        
        # There are some probes that have more
        # than 2 alignments (3,4,5), for now we
        # will discard those probes. We should
        # ask Affy.
        
        Probes <- Probes[, c(seq_len(4), 
            9:11, 8)]
        Probes[, 5] <- paste("chr", Probes[, 
            5], sep = "")
        ix <- str_count(Probes[, 6], ",")
        ix <- which(ix == 1)
        Probes <- Probes[ix, ]
        ProbSS <- matrix(as.numeric(unlist(strsplit(Probes[, 
            6], "[,-]"))), ncol = 4, byrow = 2)[, 
            2:3]
        Probes <- cbind(Probes[, seq_len(5)], 
            ProbSS, Probes[, 7:8])
        colnames(Probes) <- c("Probe ID", 
            "X Coord", "Y Coord", "Gene", 
            "Chr", "Start", "Stop", "Strand", 
            "Probe Sequence")
        # Probes<-Probes[ix,]
    }
    
    
    return(Probes)
}

#' @rdname InternalFunctions
PrepareOutput <- function(Result, Final) {
    GeneN <- sapply(seq_along(Result), function(x) {
        Result[[x]][[1]]$GeneName
    })
    GeneI <- sapply(seq_along(Result), function(x) {
        Result[[x]][[1]]$Gene
    })
    Index <- seq_along(Result)
    
    iix <- matrix(as.numeric(unlist(strsplit(Final[, 
        1], "_"))), ncol = 2, byrow = TRUE)
    
    Types <- vector("list", length = nrow(Final))
    Positions <- vector("list", length = nrow(Final))
    GeneList <- vector("list", length = nrow(Final))
    GeneID <- vector("list", length = nrow(Final))
    
    InfoX <- data.frame(GeneName = GeneN, 
        GeneID = as.numeric(GeneI), Index = as.numeric(Index), 
        stringsAsFactors = FALSE)
    
    Output <- lapply(seq_len(nrow(Final)), 
        function(x) {
            A <- InfoX[match(iix[x, 1], InfoX[, 
                2]), 3]
            B <- iix[x, 2]
            
            EventType <- Result[[A]][[B]]$Type
            
            Mat <- rbind(Result[[A]][[B]]$P1, 
                Result[[A]][[B]]$P2)
            Mat <- Mat[Mat[, 4] != 0, ]
            Mat <- Mat[Mat[, 5] != 0, ]
            Chr <- Result[[A]][[B]]$P1[1, 
                3]
            
            St <- min(Mat[, 4])
            Sp <- max(Mat[, 5])
            
            Position <- paste(Chr, ":", St, 
                "-", Sp, sep = "")
            
            
            Res <- data.frame(Gene = InfoX[A, 
                1], Event_Type = EventType, 
                Position = Position, Pvalue = Final[x, 
                  2], Zvalue = Final[x, 3], 
                stringsAsFactors = FALSE)
            rownames(Res) <- Final[x, 1]
            return(Res)
        })
    
    Output <- do.call(rbind, Output)
    Pval_Order <- order(Output[, "Pvalue"])
    Output <- Output[Pval_Order, ]
    
    return(Output)
}


#' @rdname InternalFunctions
SG_Info <- function(SG_Gene) {
    SE_Cond <- FALSE
    TTS_Cond <- FALSE
    TSS_Cond <- FALSE
    
    # Obtain information of all the elements
    # of the Graph
    
    Graph <- SGSeq:::exonGraph(SG_Gene, tx_view = FALSE)
    Graph_Nodes <- SGSeq:::nodes(Graph)
    Graph_Edges <- SGSeq:::edges(Graph)
    # Graph<-exonGraph(SG_Gene,tx_view=FALSE)
    # Graph_Nodes <- nodes(Graph)
    # Graph_Edges<- edges(Graph)
    
    Nodes_Pos <- matrix(unlist(strsplit(Graph_Nodes[, 
        2], "[:-]")), ncol = 4, byrow = TRUE)
    Edges_Pos <- matrix(unlist(strsplit(Graph_Edges[, 
        3], "[:-]")), ncol = 4, byrow = TRUE)
    
    Graph_Nodes <- cbind(Graph_Nodes[, 1], 
        Nodes_Pos, Graph_Nodes[, 3:4])
    Graph_Nodes <- as.data.frame(Graph_Nodes, 
        stringsAsFactors = FALSE)
    colnames(Graph_Nodes) <- c("Name", "Chr", 
        "Start", "End", "Strand", "Type", 
        "featureID")
    Nodes <- as.numeric(as.vector(Graph_Nodes[, 
        1]))
    
    Graph_Edges <- cbind(Graph_Edges[, c(1, 
        2)], Edges_Pos, Graph_Edges[, 4:5])
    Graph_Edges <- as.data.frame(Graph_Edges, 
        stringsAsFactors = FALSE)
    colnames(Graph_Edges) <- c("From", "To", 
        "Chr", "Start", "End", "Strand", 
        "Type", "featureID")
    
    if (as.vector(strand(SG_Gene)@values) == 
        "-") {
        Graph_Edges <- Graph_Edges[, c(2, 
            1, 3:8)]
        colnames(Graph_Edges) <- c("From", 
            "To", "Chr", "Start", "End", 
            "Strand", "Type", "featureID")
    }
    
    # Determine Subexons
    
    SE.II <- as.numeric(as.vector(Graph_Nodes[2:nrow(Graph_Nodes), 
        3])) - as.numeric(as.vector(Graph_Nodes[seq_len((nrow(Graph_Nodes) - 
        1)), 4]))
    SE.II <- which(SE.II == 1)
    L_SE.II <- length(SE.II)
    SE_chr <- rep(as.vector(Graph_Nodes[1, 
        2]), L_SE.II)
    SE_St <- Graph_Nodes[SE.II, 4]
    SE_Ed <- Graph_Nodes[SE.II + 1, 3]
    SE_std <- rep(as.vector(Graph_Nodes[1, 
        5]), L_SE.II)
    SE_type <- rep("J", L_SE.II)
    SE_FID <- rep(0, L_SE.II)
    
    SubExons <- data.frame(SE.II, SE.II + 
        1, SE_chr, SE_St, SE_Ed, SE_std, 
        SE_type, SE_FID, stringsAsFactors = FALSE)
    colnames(SubExons) <- colnames(Graph_Edges)
    
    if (nrow(SubExons) != 0) {
        Graph_Edges <- rbind(Graph_Edges, 
            SubExons)
        SE_Cond <- TRUE
    }
    
    # Determine Alternative Last
    
    TTS <- setdiff(Nodes, as.numeric(as.vector(Graph_Edges[, 
        1])))
    
    if (length(TTS) > 0) {
        TTS <- paste(TTS, ".b", sep = "")
        Ln_TTS <- length(TTS)
        
        EE <- rep("E", length(TTS))
        TTS <- data.frame(TTS, EE, rep(Graph_Edges[1, 
            3], Ln_TTS), rep(0, Ln_TTS), 
            rep(0, Ln_TTS), rep(as.vector(Graph_Edges[1, 
                6]), Ln_TTS), rep("J", Ln_TTS), 
            rep(0, Ln_TTS), stringsAsFactors = FALSE)
        colnames(TTS) <- colnames(Graph_Edges)
        
        TTS_Cond <- TRUE
    }
    
    
    # Determine Alternative First
    
    
    TSS <- setdiff(Nodes, as.numeric(as.vector(Graph_Edges[, 
        2])))
    
    if (length(TSS) > 0) {
        TSS <- paste(TSS, ".a", sep = "")
        Ln_TSS <- length(TSS)
        
        SS <- rep("S", Ln_TSS)
        TSS <- data.frame(SS, TSS, rep(Graph_Edges[1, 
            3], Ln_TSS), rep(0, Ln_TSS), 
            rep(0, Ln_TSS), rep(as.vector(Graph_Edges[1, 
                6]), Ln_TSS), rep("J", Ln_TSS), 
            rep(0, Ln_TSS), stringsAsFactors = FALSE)
        colnames(TSS) <- colnames(Graph_Edges)
        
        TSS_Cond <- TRUE
    }
    
    # Extend SG
    
    Graph_Edges[, 1] <- paste(as.vector(Graph_Edges[, 
        1]), ".b", sep = "")
    Graph_Edges[, 2] <- paste(as.vector(Graph_Edges[, 
        2]), ".a", sep = "")
    From <- paste(as.vector(Graph_Nodes[, 
        1]), ".a", sep = "")
    To <- paste(as.vector(Graph_Nodes[, 1]), 
        ".b", sep = "")
    Extended <- cbind(From, To, Graph_Nodes[, 
        2:7])
    
    Graph_Edges <- rbind(Graph_Edges, Extended)
    
    if (TSS_Cond) {
        Graph_Edges <- rbind(TSS, Graph_Edges)
    }
    
    if (TTS_Cond) {
        Graph_Edges <- rbind(Graph_Edges, 
            TTS)
    }
    
    
    rownames(Graph_Edges) <- seq_len(nrow(Graph_Edges))
    
    # Get Adjacency and Incidence Matrix
    
    GGraph <- graph_from_data_frame(Graph_Edges, 
        directed = TRUE)
    Adjacency <- as_adj(GGraph)
    
    
    Incidence <- matrix(0, nrow = ((length(Nodes) * 
        2) + 2), ncol = nrow(Graph_Edges))
    colnames(Incidence) <- rownames(Graph_Edges)
    rownames(Incidence) <- c("S", paste(rep(Nodes, 
        each = 2), c(".a", ".b"), sep = ""), 
        "E")
    
    Incidence[cbind(as.vector(Graph_Edges[, 
        "From"]), colnames(Incidence))] <- -1
    Incidence[cbind(as.vector(Graph_Edges[, 
        "To"]), colnames(Incidence))] <- 1
    
    iijj <- match(rownames(Incidence), rownames(Adjacency))
    Adjacency <- Adjacency[iijj, iijj]
    Adjacency <- as(Adjacency, "dgTMatrix")
    
    # Return All Information
    
    Result <- list(Edges = Graph_Edges, Adjacency = Adjacency, 
        Incidence = Incidence)
    # Result<-list(Edges=Graph_Edges,Incidence=Incidence)
    return(Result)
}

#' @rdname InternalFunctions
SG_creation <- function(SG_Gene) {
    SG_Gene_SoloE <- SG_Gene[type(SG_Gene) == 
        "E"]
    Ts <- unique(unlist(txName(SG_Gene_SoloE)))
    
    # Build Adjacency Matrix
    ncolAdj <- length(SG_Gene_SoloE) * 2 + 
        2
    Adj <- matrix(0, ncol = ncolAdj, nrow = ncolAdj)
    nombres <- rep(paste(seq_len((ncolAdj/2 - 
        1)), ".b", sep = ""), each = 2)
    nombresa <- rep(paste(seq_len((ncolAdj/2 - 
        1)), ".a", sep = ""), each = 2)
    nombres[seq(1, (ncolAdj - 2), by = 2)] <- nombresa[seq(1, 
        (ncolAdj - 2), by = 2)]
    rownames(Adj) <- colnames(Adj) <- c("S", 
        nombres, "E")
    for (Trans in Ts) {
        dummy <- sapply(txName(SG_Gene_SoloE) == 
            Trans, any)
        dummy2 <- rep(which(dummy) * 2, each = 2)
        dummy2[seq(2, length(dummy2), by = 2)] <- dummy2[seq(2, 
            length(dummy2), by = 2)] + 1
        dummy2 <- c(1, dummy2, ncolAdj)
        Adj[cbind(dummy2[seq_len((length(dummy2) - 
            1))], dummy2[2:length(dummy2)])] <- 1
    }
    
    PosAdj <- c(0, as.numeric(rbind(start(ranges(SG_Gene_SoloE)), 
        end(ranges(SG_Gene_SoloE)))), 0)
    
    # Build incidence matrix
    Inc <- matrix(0, nrow = ncol(Adj), ncol = sum(Adj > 
        0))
    Inc[cbind(which(Adj > 0, arr.ind = 1)[, 
        2], seq_len(ncol(Inc)))] <- 1
    Inc[cbind(which(Adj > 0, arr.ind = 1)[, 
        1], seq_len(ncol(Inc)))] <- -1
    
    rownames(Inc) <- rownames(Adj)
    colnames(Inc) <- seq_len(ncol(Inc))
    
    NuevoOrden <- which(Inc[1, ] == -1)
    NuevoOrden <- c(NuevoOrden, setdiff(unlist(apply(Inc[grep("b", 
        rownames(Inc)), ], 1, function(x) {
        A <- which(x == -1)
    })), which(Inc[nrow(Inc), ] == 1)))
    NuevoOrden <- c(NuevoOrden, unlist(apply(Inc[grep("a", 
        rownames(Inc)), ], 1, function(x) {
        which(x == -1)
    })))
    NuevoOrden <- c(NuevoOrden, which(Inc[nrow(Inc), 
        ] == 1))
    
    Inc <- Inc[, NuevoOrden]
    
    
    Adjacency <- Matrix(Adj)
    
    # Build edges data.frame
    Edges <- data.frame()
    # From <- colnames(Adj)[which(Adj>0,
    # arr.ind=1)[,1]]
    From <- rownames(Inc)[apply(Inc, 2, function(X) {
        which(X == -1)
    })]
    # To <- colnames(Adj)[which(Adj>0,
    # arr.ind=1)[,2]]
    To <- rownames(Inc)[apply(Inc, 2, function(X) {
        which(X == 1)
    })]
    Edges <- data.frame(From, To)
    Edges$Chr <- as.vector(seqnames(SG_Gene)@values)
    Edges$Start <- 0
    Edges$End <- 0
    Exons <- grep("a", Edges$From)
    # Edges$Strand <-
    # as.character(strand(SG_Gene_SoloE[Exons[1]]))
    Edges$Strand <- strand(SG_Gene)@values
    Edges$Type <- "J"
    Virtual <- c(grep("S", Edges$From), grep("E", 
        Edges$To))
    Edges[Exons, "Type"] <- "E"
    Edges[Virtual, "Type"] <- "V"  # If J is required, comment this line.
    Edges$Start <- PosAdj[match(Edges$From, 
        colnames(Adj))]
    Edges$End <- PosAdj[match(Edges$To, colnames(Adj))]
    Edges$Chr <- max(Edges$Chr)
    
    # Put featuresID
    matched <- match(Edges$End + 1e+06 * 
        Edges$Start, end(ranges(SG_Gene)) + 
        1e+06 * start(ranges(SG_Gene)))
    IndexEdge <- which(!is.na(matched))
    Edges$featureID <- 0
    Edges$featureID[IndexEdge] <- featureID(SG_Gene[matched[IndexEdge]])
    return(list(Edges = Edges, Adjacency = Adjacency, 
        Incidence = Inc))
}


#' @rdname InternalFunctions
SG_creation_RNASeq <- function(SG_Gene) {
    # which((start(SG_Gene)-end(SG_Gene))==0)
    
    SG_Gene <- SG_Gene[which((start(SG_Gene) - 
        end(SG_Gene)) != 0)]
    SG_Gene_SoloE <- SG_Gene[type(SG_Gene) == 
        "E"]
    nodos <- sort(unique(c(start(SG_Gene_SoloE), 
        end(SG_Gene_SoloE))))
    ncolAdj <- length(nodos) + 2
    Adj <- matrix(0, nrow = length(nodos) + 
        2, ncol = length(nodos) + 2)
    colnames(Adj) <- rownames(Adj) <- c("S", 
        nodos, "E")
    edges <- cbind(match(start(SG_Gene), 
        colnames(Adj)), match(end(SG_Gene), 
        colnames(Adj)))
    
    # subexones adyacentes
    adyacentes1 <- which(diff(nodos) == 1)
    edges <- rbind(edges, cbind(1 + adyacentes1, 
        adyacentes1 + 2))
    Adj[edges] <- 1
    diag(Adj) <- 0
    
    
    # Include new start and end sites based
    # on annotation
    
    if (!sum(sapply(txName(SG_Gene_SoloE), 
        length)) == 0) {
        Ts <- unique(unlist(txName(SG_Gene_SoloE)))
        Salida <- lapply(txName(SG_Gene_SoloE), 
            match, Ts)
        veces <- sapply(Salida, length)
        y <- rep(seq_len(length(veces)), 
            veces)
        x <- unlist(Salida)
        Transcripts <- matrix(FALSE, nrow = max(x), 
            ncol = max(y))
        Transcripts[cbind(x, y)] <- TRUE
        initial <- unique(max.col(Transcripts, 
            ties.method = c("first")))
        final <- unique(max.col(Transcripts, 
            ties.method = c("last")))
        colstarts <- match(start(SG_Gene_SoloE[initial]), 
            colnames(Adj))
        Adj[cbind(1, colstarts)] <- 1
        rowends <- match(end(SG_Gene_SoloE[final]), 
            rownames(Adj))
        Adj[cbind(rowends, ncol(Adj))] <- 1
    }
    
    
    # Fill orphan and widows nodes
    orphans <- which(colSums(Adj) == 0)
    orphans <- orphans[-1]
    if (length(orphans) > 0) {
        if (any(names(orphans) == "E")) {
            orphans <- orphans[-length(orphans)]
        }
        if (length(orphans) > 0) {
            Adj[cbind(1, orphans)] <- 1
        }
    }
    
    widows <- which(rowSums(Adj) == 0)
    widows <- widows[-length(widows)]
    if (length(widows) > 0) {
        Adj[cbind(widows, ncol(Adj))] <- 1
    }
    diag(Adj) <- 0
    
    # Change the name of the rows and columns
    # to 1.a 1.b 2.a 2.b,...
    
    nombres <- rep(paste(seq_len((ncolAdj/2 - 
        1)), ".b", sep = ""), each = 2)
    nombresa <- rep(paste(seq_len((ncolAdj/2 - 
        1)), ".a", sep = ""), each = 2)
    nombres[seq(1, (ncolAdj - 2), by = 2)] <- nombresa[seq(1, 
        (ncolAdj - 2), by = 2)]
    rownames(Adj) <- colnames(Adj) <- c("S", 
        nombres, "E")
    
    
    
    PosAdj <- c(0, as.numeric(rbind(start(ranges(SG_Gene_SoloE)), 
        end(ranges(SG_Gene_SoloE)))), 0)
    
    # Build incidence matrix
    Inc <- matrix(0, nrow = ncol(Adj), ncol = sum(Adj > 
        0))
    Inc[cbind(which(Adj > 0, arr.ind = 1)[, 
        2], seq_len(ncol(Inc)))] <- 1
    Inc[cbind(which(Adj > 0, arr.ind = 1)[, 
        1], seq_len(ncol(Inc)))] <- -1
    
    rownames(Inc) <- rownames(Adj)
    colnames(Inc) <- seq_len(ncol(Inc))
    
    NuevoOrden <- which(Inc[1, ] == -1)
    NuevoOrden <- c(NuevoOrden, setdiff(unlist(apply(Inc[grep("b", 
        rownames(Inc)), , drop = FALSE], 
        1, function(x) {
            A <- which(x == -1)
        })), which(Inc[nrow(Inc), ] == 1)))
    NuevoOrden <- c(NuevoOrden, unlist(apply(Inc[grep("a", 
        rownames(Inc)), , drop = FALSE], 
        1, function(x) {
            which(x == -1)
        })))
    NuevoOrden <- c(NuevoOrden, which(Inc[nrow(Inc), 
        ] == 1))
    
    Inc <- Inc[, NuevoOrden]
    
    
    Adjacency <- Matrix(Adj)
    
    # Build edges data.frame
    Edges <- data.frame()
    # From <- colnames(Adj)[which(Adj>0,
    # arr.ind=1)[,1]]
    From <- rownames(Inc)[apply(Inc, 2, function(X) {
        which(X == -1)
    })]
    # To <- colnames(Adj)[which(Adj>0,
    # arr.ind=1)[,2]]
    To <- rownames(Inc)[apply(Inc, 2, function(X) {
        which(X == 1)
    })]
    Edges <- data.frame(From, To)
    Edges$Chr <- as.vector(seqnames(SG_Gene)@values)
    Edges$Start <- 0
    Edges$End <- 0
    Exons <- grep("a", Edges$From)
    # Edges$Strand <-
    # as.character(strand(SG_Gene_SoloE[Exons[1]]))
    Edges$Strand <- strand(SG_Gene)@values
    Edges$Type <- "J"
    Virtual <- c(grep("S", Edges$From), grep("E", 
        Edges$To))
    Edges[Exons, "Type"] <- "E"
    Edges[Virtual, "Type"] <- "V"  # If J is required, comment this line.
    Edges$Start <- PosAdj[match(Edges$From, 
        colnames(Adj))]
    Edges$End <- PosAdj[match(Edges$To, colnames(Adj))]
    Edges$Chr <- max(Edges$Chr)
    
    # Put featuresID
    matched <- match(Edges$End + 1e+06 * 
        Edges$Start, end(ranges(SG_Gene)) + 
        1e+06 * start(ranges(SG_Gene)))
    IndexEdge <- which(!is.na(matched))
    Edges$featureID <- 0
    Edges$featureID[IndexEdge] <- featureID(SG_Gene[matched[IndexEdge]])
    return(list(Edges = Edges, Adjacency = Adjacency, 
        Incidence = Inc))
}


#' @rdname InternalFunctions
######## function to create Event plots in IGV
WriteGTF <- function(PATH, Data, Probes, 
    Paths) {
    STRAND <- unique(Paths[, 5])
    FILE.probes <- paste(PATH, "/probes.gtf", 
        sep = "")
    ### Error en los grep.. si no encuentra
    ### regresa 0 y no NA
    PATHS <- as.vector(unique(Probes[, 6]))
    
    
    for (i in seq_len(nrow(Probes))) {
        if (STRAND == "+") {
            START <- Probes[i, 3]
            END <- Probes[i, 3] + Probes[i, 
                4]
        } else if (STRAND == "-") {
            START <- Probes[i, 3] - Probes[i, 
                4]
            END <- Probes[i, 3]
        }
        
        # browser()
        if (Probes[i, 6] == "Ref") {
            COL <- "#B0B0B0"
        } else if (Probes[i, 6] == "Path1") {
            COL <- "#D00000"
        } else if (Probes[i, 6] == "Path2") {
            COL <- "#00CC33"
        }
        PROBES <- paste(Probes[i, 2], "\t", 
            "microarray", "\t", "probe", 
            "\t", Probes[i, 3], "\t", END, 
            "\t", "0", "\t", "*", "\t", "0", 
            "\t", "event=", Data[1, 4], "; path=", 
            Probes[i, 6], "; color=", COL, 
            ";", "probeID=", Probes[i, 1], 
            ";", sep = "")
        cat(file = FILE.probes, PROBES, sep = "\n", 
            append = TRUE)
    }
    
    
    
    # Paths GTF
    
    II <- order(Paths[, 7])
    Paths <- Paths[II, ]
    PATHS <- unique(Paths[, 7])
    FILE.paths <- paste(PATH, "/paths.gtf", 
        sep = "")
    
    
    GENESSSSS <- paste(Paths[, 1], "\t", 
        "EventPointer", "\t", "gene", "\t", 
        min(Paths[, 2]), "\t", max(Paths[, 
            3]), "\t", "0", "\t", Paths[, 
            5], "\t", "0", "\t", paste("gene_id ", 
            Data[1, 2], "_", Data[1, 3], 
            "; ", "type ", shQuote(as.vector(Data[1, 
                4]), type = "cmd"), "; color=", 
            "#000000", ";", sep = ""), sep = "")
    GENESSSSS <- unique(GENESSSSS)
    
    cat(file = FILE.paths, GENESSSSS, sep = "\n", 
        append = TRUE)
    
    
    # browser()
    for (i in seq_along(PATHS)) {
        ii <- which(Paths[, 7] == PATHS[i])
        # if (all(!is.na(grep('Ref',PATHS[i])))){
        if (length(!is.na(grep("Ref", PATHS[i]))) != 
            0) {
            COL <- "#B0B0B0"
            aaaaaa <- 3
        }
        # if (all(!is.na(match('A',PATHS[i])))){
        if (length(!is.na(grep("A", PATHS[i]))) != 
            0) {
            COL <- "#D00000"
            aaaaaa <- 2
        }
        # if (all(!is.na(match('B',PATHS[i])))){
        if (length(!is.na(grep("B", PATHS[i]))) != 
            0) {
            COL <- "#00CC33"
            aaaaaa <- 1
        }
        # if
        # (all(!is.na(match('Empty',PATHS[i])))){
        if (length(!is.na(grep("Empty", PATHS[i]))) != 
            0) {
            COL <- "#FFFFFF"
        }
        
        
        # browser() TRANS <-
        # paste(Paths[ii,1],'\t','EventPointer','\t','transcript','\t',
        # min(Paths[,2])-10*as.numeric(as.matrix(Data[2]))-aaaaaa
        # MINGENE-as.numeric(as.matrix(Data[2])),
        TRANS <- paste(Paths[ii, 1], "\t", 
            "EventPointer", "\t", "transcript", 
            "\t", min(Paths[ii, 2]), "\t", 
            max(Paths[ii, 3]), "\t", "0", 
            "\t", Paths[ii, 5], "\t", "0", 
            "\t", paste("gene_id ", Data[1, 
                2], "_", Data[1, 3], "; ", 
                "transcript_id ", shQuote(Paths[ii, 
                  7], type = "cmd"), "_", 
                gsub(" ", "_", Data[1, 4]), 
                "_", unique(Paths[ii, 6]), 
                "_", Data[1, 3], "; ", "type ", 
                shQuote(Data[1, 4], type = "cmd"), 
                "; color=", COL, ";", sep = ""), 
            sep = "")
        TRANS <- unique(TRANS)
        
        GTF <- paste(Paths[ii, 1], "\t", 
            "EventPointer", "\t", "exon", 
            "\t", Paths[ii, 2], "\t", Paths[ii, 
                3], "\t", "0", "\t", Paths[ii, 
                5], "\t", "0", "\t", paste("gene_id ", 
                Data[1, 2], "_", Data[1, 
                  3], "; ", "transcript_id ", 
                shQuote(Paths[ii, 7], type = "cmd"), 
                "_", gsub(" ", "_", Data[1, 
                  4]), "_", unique(Paths[ii, 
                  6]), "_", Data[1, 3], "; ", 
                "type ", shQuote(Data[1, 
                  4], type = "cmd"), "; exon_number ", 
                seq_len(length(ii)), "; color=", 
                COL, ";", sep = ""), sep = "")
        # if (i == 1){
        # cat(file=FILE.paths,TRANS,sep='\n')
        # }else{
        cat(file = FILE.paths, TRANS, sep = "\n", 
            append = TRUE)
        # }
        cat(file = FILE.paths, GTF, sep = "\n", 
            append = TRUE)
    }
}

#' @rdname InternalFunctions
####### function to create Event plots in IGV
WriteGTF_RNASeq <- function(PATH, Data, Paths) {
    
    # browser() Paths GTF
    STRAND <- unique(Paths[, 5])
    II <- order(Paths[, 7])
    Paths <- Paths[II, ]
    PATHS <- unique(Paths[, 7])
    FILE.paths <- paste(PATH, "/paths_RNASeq.gtf", 
        sep = "")
    
    GENESSSSS <- paste(Paths[, 1], "\t", 
        "EventPointer", "\t", "gene", "\t", 
        min(Paths[, 2]), "\t", max(Paths[, 
            3]), "\t", "0", "\t", Paths[, 
            5], "\t", "0", "\t", paste("gene_id ", 
            Data[1, 1], "; ", "type ", shQuote(as.vector(Data[1, 
                4]), type = "cmd"), "; color=", 
            "#000000", ";", sep = ""), sep = "")
    
    GENESSSSS <- unique(GENESSSSS)
    
    # browser()
    cat(file = FILE.paths, GENESSSSS, sep = "\n", 
        append = TRUE)
    
    
    # browser()
    for (i in seq_along(PATHS)) {
        ii <- which(Paths[, 7] == PATHS[i])
        # if (all(!is.na(grep('Ref',PATHS[i])))){
        if (length(!is.na(grep("Ref", PATHS[i]))) != 
            0) {
            COL <- "#B0B0B0"
            aaaaaa <- 3
        }
        # if (all(!is.na(match('A',PATHS[i])))){
        if (length(!is.na(grep("A", PATHS[i]))) != 
            0) {
            COL <- "#D00000"
            aaaaaa <- 2
        }
        # if (all(!is.na(match('B',PATHS[i])))){
        if (length(!is.na(grep("B", PATHS[i]))) != 
            0) {
            COL <- "#00CC33"
            aaaaaa <- 1
        }
        # if
        # (all(!is.na(match('Empty',PATHS[i])))){
        if (length(!is.na(grep("Empty", PATHS[i]))) != 
            0) {
            COL <- "#FFFFFF"
        }
        
        
        # browser() TRANS <-
        # paste(Paths[ii,1],'\t','EventPointer','\t','transcript','\t',min(Paths[,2])-10*as.numeric(as.matrix(Data[2]))-aaaaaa
        # MINGENE-as.numeric(as.matrix(Data[2])),
        TRANS <- paste(Paths[ii, 1], "\t", 
            "EventPointer", "\t", "transcript", 
            "\t", min(Paths[ii, 2]), "\t", 
            max(Paths[ii, 3]), "\t", "0", 
            "\t", Paths[ii, 5], "\t", "0", 
            "\t", paste("gene_id ", Data[1, 
                1], "; ", "transcript_id ", 
                shQuote(Paths[ii, 7], type = "cmd"), 
                "_", gsub(" ", "_", Data[1, 
                  4]), "_", Data[1, 1], "; ", 
                "type ", shQuote(Data[1, 
                  4], type = "cmd"), "; color=", 
                COL, ";", sep = ""), sep = "")
        TRANS <- unique(TRANS)
        
        GTF <- paste(Paths[ii, 1], "\t", 
            "EventPointer", "\t", "exon", 
            "\t", Paths[ii, 2], "\t", Paths[ii, 
                3], "\t", "0", "\t", Paths[ii, 
                5], "\t", "0", "\t", paste("gene_id ", 
                Data[1, 1], "; ", "transcript_id ", 
                shQuote(Paths[ii, 7], type = "cmd"), 
                "_", gsub(" ", "_", Data[1, 
                  4]), "_", Data[1, 1], "; ", 
                "type ", shQuote(Data[1, 
                  4], type = "cmd"), "; exon_number ", 
                seq_len(length(ii)), "; color=", 
                COL, ";", sep = ""), sep = "")
        
        
        # browser() if (i == 1){
        # cat(file=FILE.paths,TRANS,sep='\n')
        # }else{
        cat(file = FILE.paths, TRANS, sep = "\n", 
            append = TRUE)
        # }
        cat(file = FILE.paths, GTF, sep = "\n", 
            append = TRUE)
    }
}



#' @rdname InternalFunctions
#################################################################### flat2Cdf
#---------
# this function takes a 'flat' file and
# converts it to a binary CDF file it was
# downloaded from:
# http://www.aroma-project.org/howtos/create_CDF_from_scratch/
# for further details see that link.
# example:
# flat2Cdf(file='hjay.r1.flat',chipType='hjay',tag='r1,TC')
# file: assumes header...better perhaps
# to have ...  that passes to
# read.table?; requires header X, Y ucol:
# unit column gcol: group column
# col.class: column classes of file (see
# read.table); NOTE: needs check that
# right number?  splitn: parameter that
# controls the number of initial chunks
# that are unwrapped (number of
# characters of unit names used to keep
# units together for initial chunks)
# rows: cols:
flat2Cdf <- function(file, chipType, tags = NULL, 
    rows = 2560, cols = 2560, verbose = 10, 
    xynames = c("X", "Y"), gcol = 5, ucol = 6, 
    splitn = 4, col.class = c("integer", 
        "character")[c(1, 1, 1, 2, 2, 2)], 
    Directory = getwd(), ...) {
    split.quick <- function(r, ucol, splitn = 3, 
        verbose = TRUE) {
        rn3 <- substr(r[, ucol], 1, splitn)
        split.matrix <- split.data.frame
        rr <- split(r, factor(rn3))
        if (verbose) {
            cat(" split into", length(rr), 
                "initial chunks ...")
        }
        rr <- unlist(lapply(rr, FUN = function(u) split(u, 
            u[, ucol])), recursive = FALSE)
        if (verbose) {
            cat(" unwrapped into", length(rr), 
                "chunks ...")
        }
        names(rr) <- substr(names(rr), splitn + 
            2, nchar(rr))
        rr
        ## rr<-unlist(lapply(rr,FUN=function(u)
        ## split(u,u[,ucol])),recursive=FALSE,use.names=FALSE)
        ## namrr<-sapply(rr,function(u){nam<-unique(u[,ucol]);
        ## if(length(nam)>1) stop('Programming
        ## Error (units', nam,'). Please report')
        ## else return(nam)},USE.NAMES=FALSE)
        ## names(rr)<-namrr
        ## rr<-unlist(lapply(rr,FUN=function(u)
        ## split(u[,-ucol],u[,ucol])),recursive=FALSE)
        ## names(rr)<-sapply(strsplit(names(rr),'\\.'),.subset2,2)
    }
    
    if (verbose) {
        cat("Reading TXT file ...")
    }
    file <- read.table(file, header = TRUE, 
        colClasses = col.class, stringsAsFactors = FALSE, 
        comment.char = "", ...)
    if (verbose) {
        cat(" Done.\n")
    }
    
    if (verbose) {
        cat("Splitting TXT file into units ...")
    }
    gxys <- split.quick(file, ucol, splitn)
    rm(file)
    gc()
    if (verbose) {
        cat(" Done.\n")
    }
    
    l <- vector("list", length(gxys))
    if (verbose) {
        cat("Creating structure for", length(gxys), 
            "units (dot=250):\n")
    }
    for (i in seq_len(length(gxys))) {
        sp <- split(gxys[[i]], factor(gxys[[i]][, 
            gcol]))
        e <- vector("list", length(sp))
        for (j in seq_len(length(sp))) {
            np <- nrow(sp[[j]])
            e[[j]] <- list(x = sp[[j]][, 
                xynames[1]], y = sp[[j]][, 
                xynames[2]], pbase = rep("A", 
                np), tbase = rep("T", np), 
                atom = 0:(np - 1), indexpos = 0:(np - 
                  1), groupdirection = "sense", 
                natoms = np, ncellsperatom = 1)
        }
        names(e) <- names(sp)
        l[[i]] <- list(unittype = 1, unitdirection = 1, 
            groups = e, natoms = nrow(gxys[[i]]), 
            ncells = nrow(gxys[[i]]), ncellsperatom = 1, 
            unitnumber = i)
        if (verbose) {
            if (i%%250 == 0) {
                cat(".")
            }
            if (i%%5000 == 0) {
                cat("(", i, ")\n", sep = "")
            }
        }
    }
    cat("\n")
    names(l) <- names(gxys)
    if (!is.null(tags) && tags != "") {
        filename <- paste(chipType, tags, 
            sep = ",")
    } else {
        filename <- chipType
    }
    filename <- paste0(Directory, "/", filename, 
        ".cdf")
    hdr <- list(probesets = length(l), qcprobesets = 0, 
        reference = "", chiptype = chipType, 
        filename = filename, nqcunits = 0, 
        nunits = length(l), rows = rows, 
        cols = cols, refseq = "", nrows = rows, 
        ncols = cols)
    writeCdf(hdr$filename, cdfheader = hdr, 
        cdf = l, cdfqc = NULL, overwrite = TRUE, 
        verbose = verbose)
    invisible(list(cdfList = l, cdfHeader = hdr))
}


#' @rdname InternalFunctions
uniquefast <- function(X) {
    b <- X %*% rnorm(dim(X)[2])
    b <- duplicated(b)
    X <- X[!b, ]
    return(X)
}

#' @rdname InternalFunctions
filterimagine <- function(Info, paths) {
    l <- dim(Info)[1]
    tofilter <- vector(length = l)
    for (ii in seq_len(l)) {
        command <- paste0("p <- c(Info$`Path 1`[", 
            ii, "],")
        for (kk in 2:(Info$`Num of Paths`[ii] + 
            1)) {
            if (kk == (Info$`Num of Paths`[ii] + 
                1)) {
                command <- paste0(command, 
                  "Info$`Path Ref`[", ii, 
                  "])")
            } else {
                command <- paste0(command, 
                  "Info$`Path ", kk, "`[", 
                  ii, "],")
            }
        }
        eval(parse(text = command))
        p <- any(p == "-")
        tofilter[ii] <- p
    }
    
    return(which(tofilter == TRUE))
}



#' @rdname InternalFunctions
transfromedge <- function(SG, SG_Gene) {
    # La salida es un data.frame que tiene el
    # numero del edge, su start y end y a los
    # transcritos a los que pertenece
    
    
    A <- SG$Edges  # lo que viene de SG_Creation
    
    SG_Gene_SoloE <- SG_Gene[type(SG_Gene) == 
        "E"]
    B <- SG_Gene_SoloE  # El correspondiente a un Gen
    
    nn <- length(txName(B))
    aa <- c()
    for (i in seq_len(nn)) {
        aa <- c(aa, txName(B)[[i]])
    }
    aa <- unique(aa)
    
    edge <- seq_len(dim(A))[1]
    transcripts <- vector(mode = "character", 
        length = length(edge))
    
    iixe <- which(A$Type == "E")
    
    matrixexons <- matrix(0, nrow = length(aa), 
        ncol = length(iixe))
    colnames(matrixexons) <- iixe
    rownames(matrixexons) <- aa
    
    for (ii in iixe) {
        s <- A$Start[ii]
        e <- A$End[ii]
        ix <- which(start(B) == s & end(B) == 
            e)
        trans <- txName(B)[[ix]]
        n <- length(trans)
        if (n == 0) {
            transcripts[ii] <- "no mapeado en la referencia"
        } else if (n == 1) {
            matrixexons[trans, as.character(ii)] <- 1
            transcripts[ii] <- trans
        } else {
            matrixexons[trans, as.character(ii)] <- 1
            tt <- trans
            trans <- trans[1]
            for (kk in 2:n) {
                trans <- paste(trans, tt[kk], 
                  sep = "|")
            }
            transcripts[ii] <- trans
        }
    }
    
    
    iix <- which(A$Type == "V")
    transcripts[iix] <- ""
    
    
    
    # matrixexons
    iix <- which(A$Type == "J")
    for (ii in iix) {
        s <- as.character(A$From[ii])
        e <- as.character(A$To[ii])
        
        # ex1 <-
        # as.character(iixe[which(A$To[iixe]==s)])
        # ex2 <-
        # as.character(iixe[which(A$From[iixe]==e)])
        
        ex1 <- which(A$To[iixe] == s)
        ex2 <- which(A$From[iixe] == e)
        
        
        trans <- rownames(matrixexons)[which(rowSums(matrixexons[, 
            ex1:ex2]) == 2 & matrixexons[, 
            ex1] == 1 & matrixexons[, ex2] == 
            1)]
        
        n <- length(trans)
        if (n == 0) {
            transcripts[ii] <- "no mapeado en la referencia"
        } else if (n == 1) {
            transcripts[ii] <- trans
        } else {
            tt <- trans
            trans <- trans[1]
            for (kk in 2:n) {
                trans <- paste(trans, tt[kk], 
                  sep = "|")
            }
            transcripts[ii] <- trans
        }
    }
    
    From <- A$From
    To <- A$To
    Start <- A$Start
    End <- A$End
    D <- data.frame(edge = edge, From = From, 
        To = To, Start = Start, End = End, 
        transcripts = transcripts)
    return(D)
}

#' @rdname InternalFunctions
sacartranscritos <- function(edgetr, events) {
    p1 <- events$Path1
    p2 <- events$Path2
    pref <- events$PathR
    
    
    p1 <- sapply(strsplit(p1, ","), function(X) {
        n <- length(X)
        trans <- ""
        for (i in seq_len(n)) {
            ss <- strsplit(X[i], "-")[[1]][1]
            ee <- strsplit(X[i], "-")[[1]][2]
            if (ss == ee) {
                ss <- paste0(ss, ".a")
                ee <- paste0(ee, ".b")
            } else {
                ss <- paste0(ss, ".b")
                ee <- paste0(ee, ".a")
            }
            if (i == 1) {
                trans <- paste0(trans, edgetr$transcripts[which(edgetr$From == 
                  ss & edgetr$To == ee)])
            } else {
                trans <- paste(trans, edgetr$transcripts[which(edgetr$From == 
                  ss & edgetr$To == ee)], 
                  sep = "|")
            }
        }
        trans <- sapply(strsplit(trans, "\\|"), 
            function(X) {
                tt <- unique(X)
                if (length(tt) == 1) {
                  return(tt)
                } else {
                  tran <- tt[1]
                  for (j in 2:length(tt)) {
                    tran <- paste(tran, tt[j], 
                      sep = "|")
                  }
                  return(tran)
                }
            })
        
        return(trans)
    })
    
    p2 <- sapply(strsplit(p2, ","), function(X) {
        n <- length(X)
        trans <- ""
        for (i in seq_len(n)) {
            ss <- strsplit(X[i], "-")[[1]][1]
            ee <- strsplit(X[i], "-")[[1]][2]
            if (ss == ee) {
                ss <- paste0(ss, ".a")
                ee <- paste0(ee, ".b")
            } else {
                ss <- paste0(ss, ".b")
                ee <- paste0(ee, ".a")
            }
            if (i == 1) {
                trans <- paste0(trans, edgetr$transcripts[which(edgetr$From == 
                  ss & edgetr$To == ee)])
            } else {
                trans <- paste(trans, edgetr$transcripts[which(edgetr$From == 
                  ss & edgetr$To == ee)], 
                  sep = "|")
            }
        }
        
        trans <- sapply(strsplit(trans, "\\|"), 
            function(X) {
                tt <- unique(X)
                if (length(tt) == 1) {
                  return(tt)
                } else {
                  tran <- tt[1]
                  for (j in 2:length(tt)) {
                    tran <- paste(tran, tt[j], 
                      sep = "|")
                  }
                  return(tran)
                }
            })
        
        return(trans)
    })
    
    pref <- sapply(strsplit(pref, ","), function(X) {
        n <- length(X)
        trans <- ""
        for (i in seq_len(n)) {
            ss <- strsplit(X[i], "-")[[1]][1]
            ee <- strsplit(X[i], "-")[[1]][2]
            if (ss == ee) {
                ss <- paste0(ss, ".a")
                ee <- paste0(ee, ".b")
            } else {
                ss <- paste0(ss, ".b")
                ee <- paste0(ee, ".a")
            }
            if (i == 1) {
                trans <- paste0(trans, edgetr$transcripts[which(edgetr$From == 
                  ss & edgetr$To == ee)])
            } else {
                trans <- paste(trans, edgetr$transcripts[which(edgetr$From == 
                  ss & edgetr$To == ee)], 
                  sep = "|")
            }
        }
        
        trans <- sapply(strsplit(trans, "\\|"), 
            function(X) {
                tt <- unique(X)
                if (length(tt) == 1) {
                  return(tt)
                } else {
                  tran <- tt[1]
                  for (j in 2:length(tt)) {
                    tran <- paste(tran, tt[j], 
                      sep = "|")
                  }
                  return(tran)
                }
            })
        
        
        return(trans)
    })
    
    
    return(data.frame(p1 = p1, p2 = p2, ref = pref))
}

#' @rdname InternalFunctions
comprobaciontranscritos2 <- function(Result) {
    Events <- Result[, c("tran_P1", "tran_P2", 
        "tran_Ref")]
    
    comprobacion <- apply(Events, 1, function(X) {
        # condicion 0: los transcritos de un path
        # no pueden ser unicamente
        #' no mapeado en
        # la referencia' (hay q eliminarlos)
        p1 <- unlist(strsplit(as.character(X[1]), 
            "\\|"))
        p2 <- unlist(strsplit(as.character(X[2]), 
            "\\|"))
        p3 <- unlist(strsplit(as.character(X[3]), 
            "\\|"))
        
        iix1 <- which(p1 == "no mapeado en la referencia")
        if (length(iix1) > 0) {
            p1 <- p1[-iix1]
        }
        
        iix2 <- which(p2 == "no mapeado en la referencia")
        if (length(iix2) > 0) {
            p2 <- p2[-iix2]
        }
        
        iix3 <- which(p3 == "no mapeado en la referencia")
        if (length(iix3) > 0) {
            p3 <- p3[-iix3]
        }
        
        if (length(p1) > 0 & length(p2) > 
            0 & length(p3) > 0) {
            
            
            # condicion 1: los que estan en p1 no
            # pueden estar en p2
            
            cond1 <- (any(p1 %in% p2 == TRUE) | 
                any(p2 %in% p1 == TRUE))
            # Falso si se cumple la condicion (True
            # si no se cumple)
            
            # condicion 2: la suma de los que estan
            # en p1 y p2 tienen q ser los mismos que
            # los q estan en p3
            
            # p1p2 <- paste(X[1],X[2],sep='|') p1p2
            # <-
            # sapply(strsplit(p1p2,'\\|'),function(XX)
            # return(XX))
            
            # p3<-as.character(X[3]) p3 <-
            # sapply(strsplit(p3,'\\|'),
            # function(XX) return(XX))
            
            p1p2 <- c(p1, p2)
            
            cond2 <- (any(p1p2 %in% p3 == 
                FALSE) | any(p3 %in% p1p2 == 
                FALSE))
            ## Falso si se cumple la condicon (True si
            ## no se cumple)
            
            return(!(cond1 | cond2))
        } else {
            return(FALSE)
        }
    })
}




# The functios of the correction of
# SGSeq:
#' @rdname InternalFunctions
convertToSGFeatures2 <- function(x, coerce = FALSE, 
    merge = FALSE) {
    if (!is(x, "TxFeatures")) {
        stop("x must be a TxFeatures object")
    }
    if (length(x) == 0) {
        return(SGFeatures())
    }
    if (coerce) {
        features <- granges(x)
        mcols(features)$type <- as.character(type(x))
        splice5p <- mcols(features)$type %in% 
            c("I", "L")
        splice3p <- mcols(features)$type %in% 
            c("I", "F")
        splice5p[mcols(features)$type == 
            "J"] <- NA
        splice3p[mcols(features)$type == 
            "J"] <- NA
        mcols(features)$type[mcols(features)$type != 
            "J"] <- "E"
        mcols(features)$splice5p <- splice5p
        mcols(features)$splice3p <- splice3p
        mcols(features)$txName <- txName(x)
        mcols(features)$geneName <- geneName(x)
    } else {
        features <- processFeatures2(x, merge = FALSE)
    }
    features <- SGSeq:::addFeatureID(features)
    features <- SGSeq:::addGeneID(features)
    features <- SGSeq::SGFeatures(features)
    if (!coerce) {
        features <- annotate2(features, x)
    }
    return(features)
}
#' @rdname InternalFunctions
processFeatures2 <- function(features, coerce = FALSE, 
    merge = FALSE) {
    junctions <- granges(features)[type(features) == 
        "J"]
    junctions_D <- flank(junctions, -1, TRUE)
    junctions_A <- flank(junctions, -1, FALSE)
    mcols(junctions)$type <- rep("J", length(junctions))
    if (is(features, "TxFeatures")) {
        exons <- features[type(features) %in% 
            c("I", "F", "L", "U")]
        exons_D <- flank(features[type(features) %in% 
            c("I", "F")], -1, FALSE)
        exons_A <- flank(features[type(features) %in% 
            c("I", "L")], -1, TRUE)
    } else if (is(features, "SGFeatures")) {
        exons <- features[type(features) == 
            "E"]
        exons_D <- flank(features[splice3p(features)], 
            -1, FALSE)
        exons_A <- flank(features[splice5p(features)], 
            -1, TRUE)
    }
    exons <- granges(exons)
    exons_D <- granges(exons_D)
    exons_A <- granges(exons_A)
    D <- unique(c(junctions_D, exons_D))
    mcols(D)$type <- rep("D", length(D))
    A <- unique(c(junctions_A, exons_A))
    mcols(A)$type <- rep("A", length(A))
    splicesites <- c(D, A)
    other <- c(junctions, splicesites)
    exons <- disjoin(exons)
    exons_start <- flank(exons, -1, TRUE)
    exons_end <- flank(exons, -1, FALSE)
    i_q <- which(!exons_end %over% splicesites)
    i_s <- which(!exons_start %over% splicesites)
    ol <- findOverlaps(suppressWarnings(flank(exons[i_q], 
        1, FALSE)), exons_start[i_s])
    if ((length(ol) > 0)) {
        qH <- i_q[queryHits(ol)]
        sH <- i_s[subjectHits(ol)]
        i_to_be_merged <- union(qH, sH)
        d <- data.frame(from = qH, to = sH)
        v <- data.frame(name = i_to_be_merged)
        g <- graph.data.frame(d = d, directed = TRUE, 
            vertices = v)
        k <- clusters(g)$membership
        exons_to_be_merged <- split(exons[i_to_be_merged], 
            k)
        exons_merged <- unlist(reduce(exons_to_be_merged))
        if (length(exons_to_be_merged) != 
            length(exons_merged)) {
            stop("cannot merge non-adjacent exons")
        }
        if (merge) {
            exons <- c(exons[-i_to_be_merged], 
                exons_merged)
        }
    }
    exons_start <- flank(exons, -1, TRUE)
    exons_end <- flank(exons, -1, FALSE)
    splice5p <- rep(FALSE, length(exons))
    i_spliced <- unique(queryHits(findOverlaps(exons_start, 
        A)))
    i_adjacent <- unique(queryHits(findOverlaps(suppressWarnings(flank(exons, 
        1, TRUE)), exons)))
    splice5p[setdiff(i_spliced, i_adjacent)] <- TRUE
    splice3p <- rep(FALSE, length(exons))
    i_spliced <- unique(queryHits(findOverlaps(exons_end, 
        D)))
    i_adjacent <- unique(queryHits(findOverlaps(suppressWarnings(flank(exons, 
        1, FALSE)), exons)))
    splice3p[setdiff(i_spliced, i_adjacent)] <- TRUE
    mcols(exons)$type <- rep("E", length(exons))
    mcols(exons)$splice5p <- splice5p
    mcols(exons)$splice3p <- splice3p
    mcols(other)$splice5p <- rep(NA, length(other))
    mcols(other)$splice3p <- rep(NA, length(other))
    features <- setNames(c(exons, other), 
        NULL)
    features <- sort(features)
    return(features)
}
#' @rdname InternalFunctions
annotate2 <- function(query, subject) {
    # query <- features subject <- a
    if (!is(subject, "TxFeatures")) {
        stop("subject must be a TxFeatures object")
    }
    if (is(query, "SGFeatures")) {
        query <- annotateFeatures2(query, 
            subject)
    } else if (is(query, "SGVariants")) {
        query <- updateObject(query, verbose = TRUE)
        query_class <- class(query)
        query_mcols <- mcols(query)
        query_unlisted <- unlist(query, use.names = FALSE)
        extended <- addDummySpliceSites(query_unlisted)
        extended <- annotate(extended, subject)
        i <- match(featureID(query_unlisted), 
            featureID(extended))
        query_unlisted <- extended[i]
        query <- relist(query_unlisted, query)
        mcols(query) <- query_mcols
        query <- new(query_class, query)
        query <- annotatePaths(query)
    } else if (is(query, "Counts")) {
        rd <- rowRanges(query)
        rd <- annotate(rd, subject)
        rowRanges(query) <- rd
    }
    return(query)
}
#' @rdname InternalFunctions
annotateFeatures2 <- function(query, subject) {
    # query <- features subject <- a
    i <- which(type(subject) %in% c("F", 
        "L"))
    if (length(i) > 0) {
        subject <- c(subject[-i], mergeExonsTerminal2(subject[i], 
            1))
    }
    if (is(query, "TxFeatures")) {
        hits <- matchTxFeatures(query, subject)
    } else if (is(query, "SGFeatures")) {
        hits <- SGSeq:::matchSGFeatures(query, 
            subject)
    }
    qH <- queryHits(hits)
    sH <- subjectHits(hits)
    for (option in c("tx", "gene")) {
        q_id <- factor(slot(query, "featureID"))
        s_ann <- slot(subject, paste0(option, 
            "Name"))
        id_ann <- SGSeq:::splitCharacterList(s_ann[sH], 
            q_id[qH])
        q_ann <- setNames(id_ann[match(q_id, 
            names(id_ann))], NULL)
        slot(query, paste0(option, "Name")) <- q_ann
    }
    if (is(query, "SGFeatures")) {
        query2 <- SGSeq:::propagateAnnotation(query)
    }
    return(query)
}
#' @rdname InternalFunctions
mergeExonsTerminal2 <- function(features, 
    min_n_sample = 1) {
    # features <- subject[i] min_n_sample <-
    # 1
    index <- which(type(features) %in% c("F", 
        "L"))
    if (length(index) > 0) {
        features <- features[index]
        splicesite <- SGSeq:::feature2name(features, 
            collapse_terminal = FALSE)
        # here is where is merged the starts and
        # ends.  collapse_terminal = TRUE y ahora
        # es FALSE
        splicesite_n <- table(splicesite)
        i <- which(splicesite %in% names(which(splicesite_n >= 
            min_n_sample)))
        features <- features[i]
        splicesite <- splicesite[i]
        splicesite <- factor(splicesite)
        splicesite_i <- split(seq_along(features), 
            splicesite)
        splicesite_w <- split(width(features), 
            splicesite)
        splicesite_i <- mapply(function(i, 
            w) {
            i[which.max(w)]
        }, i = splicesite_i, w = splicesite_w, 
            SIMPLIFY = TRUE)
        exons <- features[splicesite_i]
        for (ann in c("txName", "geneName")) {
            exons_ann <- SGSeq:::splitCharacterList(slot(features, 
                ann), splicesite)
            slot(exons, ann) <- setNames(exons_ann, 
                NULL)
        }
    } else {
        si <- seqinfo(features)
        exons <- TxFeatures()
        seqinfo(exons) <- si
    }
    return(exons)
}

########### function for the bootstrap statistic----




#' @rdname InternalFunctions
get_beta <- function(combboots, incrPSI_original, 
    ncontrastes) {
    newcombboots <- rep(list(matrix(NA, dim(combboots[[1]])[1], 
        dim(combboots[[1]])[2])), ncontrastes)
    
    newcombboots <- lapply(combboots, function(X) {
        return((1 + X)/2)  # Make values between 0-1
    })
    
    # Data for beta distribution
    rmedia <- sapply(newcombboots, rowMeans)
    rvar <- sapply(newcombboots, rowVars) + 
        1e-05
    alpha <- ((1 - rmedia) * (rmedia^2)/rvar) - 
        rmedia
    beta <- alpha * (1 - rmedia)/rmedia
    
    # Example of plots of an event (any
    # event) ---- n<-391 #Index of event
    # hist(newcombboots[[1]][n,],seq(0,1,by =
    # 0.01),freq = F, main = paste('Histogram
    # of increase in PSI'), xlab =
    # expression(paste('Increase in PSI')))
    # Histogram x <- seq(0,1, by =.001)
    # lines(x,dbeta(x,alpha[n],beta[n]),type
    # ='l', lwd = 1, col ='orange') #Density
    # function with calculated alpha and beta
    # lines(density(newcombboots[[1]][n,]),
    # type ='l', lwd = 1, col ='purple')
    # #Empirical density function
    
    # Obtain p-values for all events ----
    deltaPSI <- t(incrPSI_original)
    pvalues <- matrix(NA, nrow = dim(deltaPSI)[1], 
        ncol = dim(deltaPSI)[2])
    
    pvalues <- pbeta(deltaPSI, alpha, beta)
    positionMa <- which(pvalues > 0.5)
    pvalues[positionMa] <- pbeta(deltaPSI[positionMa], 
        alpha[positionMa], beta[positionMa], 
        lower.tail = FALSE)
    pvalues <- pvalues * 2
    
    
    deltaPSI <- (deltaPSI * 2) - 1
    result <- list(deltaPSI = deltaPSI, pvalues = pvalues)
    
    return(result)
}

#' @rdname InternalFunctions
get_table <- function(PSI_arrayP, nevents, 
    totchunk, chunk, nsamples, incrPSI_original, 
    V, nboot, nbootin, ncontrastes) {
    
    # Obtain the part of the incrPSI_original
    # needed for the minichunk and its length
    indexincr <- match(rownames(PSI_arrayP), 
        colnames(incrPSI_original))
    incrPSI_originalChunk <- incrPSI_original[, 
        indexincr, drop = FALSE]
    l <- length(indexincr)  # Length of the minichunk
    
    combboots <- rep(list(matrix(NA, l, nboot * 
        nbootin)), ncontrastes)
    # Intialize matrix for the increase in
    # PSI
    
    I <- as.integer(rep(seq_len(l), nsamples))
    J <- as.integer(rep(seq_len(nsamples), 
        each = l))
    CTEind <- I + (J - 1L) * l - 1L * nsamples * 
        l
    # Constant needed for function get_YB
    output <- matrix(NA, l, ncontrastes)
    gc()
    for (boot in seq_len(nboot)) {
        Yb <- get_YB(PSI_arrayP, l, nsamples, 
            I, J, CTEind)
        # Obtain the combination of bootstraps,
        # the matrix contains the values of PSI
        
        for (boot2 in seq_len(nbootin)) {
            Yb1 <- Yb[, sample(ncol(Yb), 
                replace = TRUE)]
            # Samples the Yb (mixes data)
            output <- tcrossprod(V, Yb1)  # Obtain the increase in PSI
            for (boot3 in seq_len(ncontrastes)) {
                combboots[[boot3]][, boot2 + 
                  nbootin * (boot - 1)] <- output[boot3, 
                  ]  # Fills matrix of increase in PSI
            }
        }
    }
    
    # Obtain p-values and table of the
    # minichunks
    table <- get_beta(combboots, incrPSI_originalChunk, 
        ncontrastes)
    # matplot(t(PSI_original[order(pvalues)[1:30],]),
    # type='b') Plot of events with # small
    # p-value
    
    
    return(table)
}

#' @rdname InternalFunctions
get_YB <- function(PSI_arrayS, l, nsamples, 
    I, J, CTEind) {
    K <- rep(sample(dim(PSI_arrayS)[3], nsamples, 
        replace = TRUE), l)
    
    # Formula to obtain the index of the
    # PSI_arrayS
    IndiceIJK <- CTEind + K * prod(dim(PSI_arrayS)[c(1, 
        2)])
    
    Yb <- PSI_arrayS[IndiceIJK]
    attr(Yb, "dim") <- c(l, nsamples)
    return(Yb)
}

#' @rdname InternalFunctions
getInfo <- function(table, ncontrast) {
    if (class(table$LocalFDR) == "list") {
        data <- data.frame(deltaPSI = table$deltaPSI[, 
            ncontrast], Pvalues = table$Pvalues[, 
            ncontrast], LocalFDR = table$LocalFDR[[ncontrast]]$lfdr2)
    } else {
        data <- data.frame(deltaPSI = table$deltaPSI[, 
            ncontrast], Pvalues = table$Pvalues[, 
            ncontrast], LocalFDR = table$LocalFDR$lfdr2)
    }
    
    return(data)
}



########### function for primers design ----

#' @rdname InternalFunctions
PrimerSequenceGeneral <- function(taqman,FinalExons,
                                  generaldata, SG,Dir,
                                  nPrimers,
                                  Primer3Path=Sys.which("primer3_core"),
                                  maxLength,
                                  minsep,wminsep,
                                  valuethreePpenalty,
                                  wnpaths,qualityfilter)
  {
  thermo.param = file.path(Dir, "primer3_config")
  settings = file.path(Dir,"primer3web_v4_0_0_default_settings.txt")
  PrimersFound <- 0
  n <- 0
  Fdata <- data.frame()
  if (class(FinalExons)!="character"){
    while ((nrow(Fdata) < nPrimers )& (n<nrow(FinalExons))){
      n <- n+1
      # 1) Use only two primers
      if(FinalExons[n,5]==0){
        
        Fdata1 <-PrimerSequenceTwo(FinalExons,SG,generaldata,n,thermo.param,
                                   Primer3Path,settings)
        Fdata <- rbind(Fdata, Fdata1)
        
        # 2) Use common FW primer and two Reverse primers
      } else if (FinalExons[n,5]==1){
        Fdata1 <- PrimerSequenceCommonFor(FinalExons,SG,generaldata,n,
                                          thermo.param,Primer3Path,
                                          settings)
        Fdata <- rbind(Fdata, Fdata1)
        
        
        # 3) Use common reverse primer and two FW primers    
      } else if (FinalExons[n,5]==2){
        Fdata1 <- PrimerSequenceCommonRev(FinalExons,SG,generaldata,n,
                                          thermo.param,Primer3Path,
                                          settings)
        Fdata <- rbind(Fdata, Fdata1)
      }
    }
    if (dim(Fdata)[1]!=0){
      # Sorting results taking into account sequences:
      RankedFdata <- getranksequence(taqman, Fdata,maxLength,minsep,
                                     wminsep,valuethreePpenalty,wnpaths,
                                     qualityfilter)
    }else{
      RankedFdata <- "Not possible to place primers due to the structure of the Event."
      
    }
  } else{
    RankedFdata <- "Not possible to place primers due to the structure of the Event."
  }
  return(RankedFdata)
}


#' @rdname InternalFunctions
PrimerSequenceTwo <-function(FinalExons,SG,
                             generaldata,n,
                             thermo.param,
                             Primer3Path,settings)
  {
  Fdata1<- data.frame()
  # Get the ID for each exon where Primers are placed
  ExonF1 <- match(as.character(FinalExons[n,1]), SG$Edges$From)
  ExonR1 <- match(as.character(FinalExons[n,3]), SG$Edges$From)
  # Get de sequence of each exon:
  Hsapienshg38  <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  FExonSeq <- getSeq(Hsapienshg38,SG$Edges$Chr[ExonF1], as.numeric(SG$Edges$Start[ExonF1]), as.numeric(SG$Edges$End[ExonF1]))
  RExonSeq <- getSeq(Hsapienshg38,SG$Edges$Chr[ExonR1], as.numeric(SG$Edges$Start[ExonR1]), as.numeric(SG$Edges$End[ExonR1]))
  
  # Build input sequence for Primer3
  minlength <- max(c(str_length(FExonSeq),str_length(RExonSeq)))+1
  seq <- paste(as.character(FExonSeq), paste(rep("N",minlength),collapse = "") ,as.character(RExonSeq),sep="")
  maxlength <- str_length(seq)
  # Call Primer3 
  p1 <- callPrimer3(seq, size_range = sprintf("%0.0f-%0.0f", minlength, maxlength),
                    name = "Primer1", Primer3Path = Primer3Path,
                    thermo.param = thermo.param,
                    sequence_target = 110+minlength,
                    settings = settings)
  
  # Build Output binding with new Sequence info and knowed Exon info
  if(is.null(dim(p1))==FALSE){
    for (s in 1: dim(p1)[1]) {
      For1Exon <- FinalExons[n,1]
      Rev1Exon <- FinalExons[n,3]
      For1Seq <- p1[s,2]
      Rev1Seq <- p1[s,3]
      LastPosFor1 <-  p1[s,6]+ p1[s,7] - 1
      LastPosFor2 <- NA
      FirstPosRev1 <- p1[s,8]-str_length(FExonSeq)-minlength-p1[s,9]+1
      FirstPosRev2 <- NA
      distinPrimers <- p1$PRIMER_PAIR_PRODUCT_SIZE[s] - minlength
      Info <- getDistanceseachPath(For1Exon,Rev1Exon,generaldata,distinPrimers,SG)
      
      paste(Info[which(Info=="p0")-1][order(as.integer(Info[which(Info=="p0")-1]),decreasing = FALSE)],collapse = " - ")
      
      # Distances from any path:
      DistancesP0 <- paste(Info[which(Info=="p0")-1][order(as.integer(Info[which(Info=="p0")-1]),decreasing = FALSE)],collapse = " - ")
      # Distances from path1:
      DistancesP1 <- paste(Info[which(Info=="p1")-1][order(as.integer(Info[which(Info=="p1")-1]),decreasing = FALSE)],collapse = " - ")
      # Distances from path2
      DistancesP2 <- paste(Info[which(Info=="p2")-1][order(as.integer(Info[which(Info=="p2")-1]),decreasing = FALSE)],collapse = " - ")
      Distances <- data.frame(DistancesP1,DistancesP2,DistancesP0)
      names(Distances) <-c("DistPath1","DistPath2","DistNoPath")
      PrSeq <- data.frame(For1Seq,NA,Rev1Seq,NA,LastPosFor1,LastPosFor2,FirstPosRev1,FirstPosRev2)
      colnames(PrSeq)<-c("For1Seq","For2Seq","Rev1Seq","Rev2Seq","LastPosFor1","LastPosFor2","FirstPosRev1","FirstPosRev2")
      Fdata1<-cbind(PrSeq , FinalExons[n,-6],Distances)
    }
  }   
  return(Fdata1)
}


#' @rdname InternalFunctions
ProbesSequence <- function(SG,FinalSeq,generaldata,Dir
                           ,Primer3Path=Sys.which("primer3_core"),
                           nProbes)
  {
  if (class(FinalSeq)!="character"){
    thermo.param = file.path(Dir, "primer3_config")
    settings = file.path(Dir,"primer3web_v4_0_0_default_settings.txt")
    # The Probe will be in the reference and in one of the two paths:
    ProbesSeq<- data.frame(cbind(rep(NA,nrow(FinalSeq)),rep(NA,nrow(FinalSeq)),rep(NA,nrow(FinalSeq))))
    names(ProbesSeq) <- c("SeqProbeRef","SeqProbeP1","SeqProbeP2")
    nProbescount<- 0
    n <- 0
    while (n<nrow(FinalSeq) && nProbescount<nProbes){
      n <- n + 1
      
      # We put a probe in the reference:
      ExonsRef <- generaldata$exonsPathsandRef$Reference
      seqref <- CreateSequenceforProbe(SG,ExonsRef,FinalSeq,n)
      probesref <- callPrimer3probes(seqref, name = "Primer1", 
                                     Primer3Path = Primer3Path,
                                     thermo.param = thermo.param,
                                     sequence_target = 20,
                                     settings = settings)
      if (length(probesref)>32){
        ProbesSeq[n,1] <-  probesref[16,2]
        
        # We put a probe in the path1:
        
        ExonsPath1 <- generaldata$exonsPathsandRef$Path1
        seqPath1 <- CreateSequenceforProbe(SG,ExonsPath1,FinalSeq,n)
        probesPath1 <- callPrimer3probes(seqPath1,name = "Primer1",
                                         Primer3Path = Primer3Path,
                                         thermo.param = thermo.param,
                                         sequence_target = 20,
                                         settings = settings)
        if (nrow(probesPath1)>24){
          ProbesSeq[n,2] <-  probesPath1[16,2]
          nProbescount <- nProbescount + 1
        }else{
          # We try to put a probe in the path2:
          ExonsPath2 <- generaldata$exonsPathsandRef$Path2
          seqPath2 <- CreateSequenceforProbe(SG,ExonsPath2,FinalSeq,n)
          probesPath2 <- callPrimer3probes(seqPath2,name = "Primer1",
                                           Primer3Path = Primer3Path,
                                           thermo.param = thermo.param,
                                           sequence_target = 20,
                                           settings = settings)
          if (nrow(probesPath2)>24){
            ProbesSeq[n,3] <-  probesPath2[16,2]
            nProbescount <- nProbescount + 1
          }
        }
      }
    }
  }else{
    ProbesSeq <- "Not possible to place probes due to the structure of the Event."
  }
  return(ProbesSeq)
}


#' @rdname InternalFunctions
sort.exons <- function(namesPath, decreasing = FALSE)
  {
  Indices <- order(as.numeric(unlist(strsplit(namesPath,".", fixed=TRUE))[c(TRUE,FALSE)]),
                   decreasing = decreasing)
  return(namesPath[Indices])
}


#' @rdname InternalFunctions
all_simple_paths2 <- function(wg,from,to,...)
  {
  Adj <- as_adjacency_matrix(wg, attr="weight")
  diag(Adj) <- 0.001
  i1 <- match(from, colnames(Adj))
  i2 <- match(to, colnames(Adj))
  Adj <- Adj[i1:i2,i1:i2, drop = FALSE]
  wg <- graph_from_adjacency_matrix(Adj, weighted=TRUE)
  allPaths <- all_simple_paths(wg,from,to,...)
  return(allPaths)
}


#' @rdname InternalFunctions
callPrimer3 <- function (seq,threeprimers = FALSE,
                         pr,reverse=FALSE ,size_range = "150-500",
                         Tm = c(57, 59, 62), name = "Primer1", 
                         Primer3Path = "primer3-2.3.7/bin/primer3_core",
                         thermo.param = "primer3-2.3.7/src/primer3_config/", 
                         sequence_target = NULL,
                         settings = "primer3-2.3.7/primer3web_v4_0_0_default_settings.txt")
{
  if (sum(file.exists(Primer3Path, thermo.param, settings)) != 
      3) {
    message("Please check your Primer3 paths!")
    return(NULL)
  }
  
  #"/" needed to be added at the end of thermo.param so windows could detect it
  thermo.param <- paste(thermo.param,"/",sep="")
  p3.input = tempfile()
  p3.output = tempfile()
  if (threeprimers==TRUE){
    if (reverse==TRUE)
    {
      cmmd <- paste(sprintf("SEQUENCE_ID=%s\n", name), sprintf("SEQUENCE_TEMPLATE=%s\n", 
                                                               as.character(seq)),sprintf("SEQUENCE_PRIMER_REVCOMP=%s\n", 
                                                                                          as.character(pr)), "PRIMER_TASK=pick_detection_primers\n", 
                    "PRIMER_PICK_LEFT_PRIMER=1\n", "PRIMER_PICK_INTERNAL_OLIGO=0\n", 
                    "PRIMER_PICK_RIGHT_PRIMER=1\n", "PRIMER_EXPLAIN_FLAG=1\n", 
                    "PRIMER_PAIR_MAX_DIFF_TM=5\n", sprintf("PRIMER_MIN_TM=%s\n", 
                                                           Tm[1]), sprintf("PRIMER_OPT_TM=%s\n", Tm[2]), sprintf("PRIMER_MAX_TM=%s\n", 
                                                                                                                 Tm[3]), sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n", 
                                                                                                                                 size_range), sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n", 
                                                                                                                                                      thermo.param), "=", sep = "")
    }else{
      cmmd <- paste(sprintf("SEQUENCE_ID=%s\n", name), sprintf("SEQUENCE_TEMPLATE=%s\n", 
                                                               as.character(seq)),sprintf("SEQUENCE_PRIMER=%s\n", 
                                                                                          as.character(pr)), "PRIMER_TASK=pick_detection_primers\n", 
                    "PRIMER_PICK_LEFT_PRIMER=1\n", "PRIMER_PICK_INTERNAL_OLIGO=0\n", 
                    "PRIMER_PICK_RIGHT_PRIMER=1\n", "PRIMER_EXPLAIN_FLAG=1\n", 
                    "PRIMER_PAIR_MAX_DIFF_TM=5\n", sprintf("PRIMER_MIN_TM=%s\n", 
                                                           Tm[1]), sprintf("PRIMER_OPT_TM=%s\n", Tm[2]), sprintf("PRIMER_MAX_TM=%s\n", 
                                                                                                                 Tm[3]), sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n", 
                                                                                                                                 size_range), sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n", 
                                                                                                                                                      thermo.param), "=", sep = "")
    }
    
  }else if(threeprimers == FALSE){
    cmmd <- paste(sprintf("SEQUENCE_ID=%s\n", name), sprintf("SEQUENCE_TEMPLATE=%s\n", 
                                                             as.character(seq)), "PRIMER_TASK=pick_detection_primers\n", 
                  "PRIMER_PICK_LEFT_PRIMER=1\n", "PRIMER_PICK_INTERNAL_OLIGO=0\n", 
                  "PRIMER_PICK_RIGHT_PRIMER=1\n", "PRIMER_EXPLAIN_FLAG=1\n", 
                  "PRIMER_PAIR_MAX_DIFF_TM=5\n", sprintf("PRIMER_MIN_TM=%s\n", 
                                                         Tm[1]), sprintf("PRIMER_OPT_TM=%s\n", Tm[2]), sprintf("PRIMER_MAX_TM=%s\n", 
                                                                                                               Tm[3]), sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n", 
                                                                                                                               size_range), sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n", 
                                                                                                                                                    thermo.param), "=", sep = "")
    if (!is.null(sequence_target)) {
      cmmd <- paste(sprintf("SEQUENCE_TARGET=%s,2\n", sequence_target), 
                    cmmd)
    }
  }
  write(cmmd, p3.input)
  if(Sys.info()[1]=="Windows") {
    #In Windows only shell worked properly at calling the system 
    shell(paste(Primer3Path," ", p3.input, " -p3_settings_file \"", settings, "\" > ", p3.output, sep =""))
  } else {
    #system worked properly on Mac or linux
    system(paste(Primer3Path," ", p3.input, " -p3_settings_file \"", settings, "\" > ", p3.output, sep =""))
  }  
  out <- read.delim(p3.output, sep = "=", header = FALSE)
  unlink(c(p3.input, p3.output))
  returned.primers = as.numeric(as.vector(out[out[, 1] == "PRIMER_PAIR_NUM_RETURNED", 
                                              ][, 2]))
  if (length(returned.primers) == 0) {
    warning("primers not detected for ", name, call. = FALSE)
    return(NA)
  }
  if ((returned.primers) == 0) {
    warning("primers not detected for ", name, call. = FALSE)
    return(NA)
  }
  if (returned.primers > 0) {
    designed.primers = data.frame()
    for (i in seq(0, returned.primers - 1, 1)) {
      id = sprintf("PRIMER_LEFT_%i_SEQUENCE", i)
      PRIMER_LEFT_SEQUENCE = as.character(out[out[, 1] == 
                                                id, ][, 2])
      id = sprintf("PRIMER_RIGHT_%i_SEQUENCE", i)
      PRIMER_RIGHT_SEQUENCE = as.character(out[out[, 1] == 
                                                 id, ][, 2])
      id = sprintf("PRIMER_LEFT_%i", i)
      PRIMER_LEFT = as.numeric(unlist(strsplit(as.vector((out[out[, 
                                                                  1] == id, ][, 2])), ",")))
      id = sprintf("PRIMER_RIGHT_%i", i)
      PRIMER_RIGHT = as.numeric(unlist(strsplit(as.vector((out[out[, 
                                                                   1] == id, ][, 2])), ",")))
      id = sprintf("PRIMER_LEFT_%i_TM", i)
      PRIMER_LEFT_TM = as.numeric(as.vector((out[out[, 
                                                     1] == id, ][, 2])), ",")
      id = sprintf("PRIMER_RIGHT_%i_TM", i)
      PRIMER_RIGHT_TM = as.numeric(as.vector((out[out[, 
                                                      1] == id, ][, 2])), ",")
      res = out[grep(i, out[, 1]), ]
      extra.inf = t(res)[2, , drop = FALSE]
      colnames(extra.inf) = sub(paste("_", i, sep = ""), 
                                "", res[, 1])
      extra.inf = extra.inf[, -c(4:9), drop = FALSE]
      extra.inf = apply(extra.inf, 2, as.numeric)
      primer.info = data.frame(i, PRIMER_LEFT_SEQUENCE, 
                               PRIMER_RIGHT_SEQUENCE, PRIMER_LEFT_TM, PRIMER_RIGHT_TM, 
                               PRIMER_LEFT_pos = PRIMER_LEFT[1], PRIMER_LEFT_len = PRIMER_LEFT[2], 
                               PRIMER_RIGHT_pos = PRIMER_RIGHT[1], PRIMER_RIGHT_len = PRIMER_RIGHT[2], 
                               t(data.frame(extra.inf)), stringsAsFactors = FALSE)
      rownames(primer.info) = NULL
      designed.primers = rbind(designed.primers, primer.info)
    }
  }
  return(designed.primers)
}

#' @rdname InternalFunctions
callPrimer3probes <- function (seq, name = "Primer1",
                               Primer3Path = "primer3-2.3.7/bin/primer3_core",
                               thermo.param = "primer3-2.3.7/src/primer3_config/", 
                               sequence_target = NULL, 
                               settings = "primer3-2.3.7/primer3web_v4_0_0_default_settings.txt")
{
  if (sum(file.exists(Primer3Path, thermo.param, settings)) != 3) {
    message("Please check your Primer3 paths!")
    return(NULL)
  }
  
  #"/" needed to be added at the end of thermo.param so windows could detect it
  thermo.param <- paste(thermo.param,"/",sep="")
  p3.input = tempfile()
  p3.output = tempfile()
  
  cmmd <- paste(sprintf("SEQUENCE_ID=%s\n", name), 
                sprintf("SEQUENCE_TEMPLATE=%s\n", as.character(seq)),  "PRIMER_TASK=generic\n",
                "PRIMER_PICK_LEFT_PRIMER=0\n", "PRIMER_PICK_INTERNAL_OLIGO=1\n", 
                "PRIMER_PICK_RIGHT_PRIMER=0\n" ,"PRIMER_NUM_RETURN=2\n", 
                sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n",thermo.param),"=", sep = "")
  if (!is.null(sequence_target)) {
    cmmd <- paste(sprintf("SEQUENCE_TARGET=%s,2\n", sequence_target), 
                  cmmd)
  }
  
  write(cmmd, p3.input)
  if(Sys.info()[1]=="Windows") {
    #In Windows only shell worked properly at calling the system 
    shell(paste(Primer3Path," ", p3.input, " -p3_settings_file \"", settings, "\" > ", p3.output, sep =""))
  } else {
    #system worked properly on Mac or linux
    system(paste(Primer3Path," ", p3.input, " -p3_settings_file \"", settings, "\" > ", p3.output, sep =""))
  }  
  out <- read.delim(p3.output, sep = "=", header = FALSE)
  unlink(c(p3.input, p3.output))
  
  out <- as.matrix(out)
  return(out)
}


#' @rdname InternalFunctions
CreateSequenceforProbe <- function(SG,Exons,FinalSeq,n)
  {
  if (isEmpty(Exons)){
    Exons=""}
  seq <- ""
  seqinicio <- ""
  seqfin <- ""
  # We need to take into account that the probe has to be between primers:
  # We compre with For1Exon:
  compare <- match(Exons,as.character(FinalSeq[n,9]))
  if (length(which(compare==1))>0){
    if (length(Exons) ==which(as.logical(compare))){
      Exons <- ""
    }else {
      Exons=Exons[(which(as.logical(compare))+1) :length(Exons)]
    }
    # If there is a primer in an exon we need to put only the part of the sequence 
    # after that primer:
    ExonID <- match(FinalSeq[n,9], SG$Edges$From)
    ExonSeq <- as.character(getSeq(Hsapiens,SG$Edges$Chr[ExonID], as.numeric(SG$Edges$Start[ExonID]), as.numeric(SG$Edges$End[ExonID])))
    ExonSeq <- unlist(strsplit(ExonSeq,""))
    seqinicio <- paste(ExonSeq[(FinalSeq[n,5]+1):length(ExonSeq)],collapse = "") 
  }
  # We compre with For2Exon:
  compare <- match(Exons,as.character(FinalSeq[n,10]))
  if (length(which(compare==1))>0){
    if (length(Exons)==which(as.logical(compare))){
      Exons <- ""
    }else {
      Exons=Exons[(which(as.logical(compare))+1) :length(Exons)]
    }
    # If there is a primer in an exon we need to put only the part of the sequence 
    # after that primer:
    ExonID <- match(FinalSeq[n,10], SG$Edges$From)
    ExonSeq <- as.character(getSeq(Hsapiens,SG$Edges$Chr[ExonID], as.numeric(SG$Edges$Start[ExonID]), as.numeric(SG$Edges$End[ExonID])))
    ExonSeq <- unlist(strsplit(ExonSeq,""))
    seqinicio <- paste(ExonSeq[(FinalSeq[n,6]+1):length(ExonSeq)],collapse = "") 
  }
  # We compre with Rev1Exon:
  compare <- match(Exons,FinalSeq[n,11])
  if (length(which(compare==1))>0){
    if (1==which(as.logical(compare))){
      Exons <- ""
    }else {
      Exons=Exons[0:(which(as.logical(compare))-1)]
    }
    # If there is a primer in an exon we need to put only the part of the sequence 
    # after that primer:
    ExonID <- match(FinalSeq[n,11], SG$Edges$From)
    ExonSeq <- as.character(getSeq(Hsapiens,SG$Edges$Chr[ExonID], as.numeric(SG$Edges$Start[ExonID]), as.numeric(SG$Edges$End[ExonID])))
    ExonSeq <- unlist(strsplit(ExonSeq,""))
    seqfin <- paste(ExonSeq[0:(FinalSeq[n,7]-1)],collapse = "") 
  }
  # We compre with Rev2Exon:
  compare <- match(Exons,FinalSeq[n,12])
  if (length(which(compare==1))>0){
    if (1==which(as.logical(compare))){
      Exons <- ""
    }else {
      Exons=Exons[(which(as.logical(compare))+1) :length(Exons)]
    }
    # If there is a primer in an exon we need to put only the part of the sequence 
    # after that primer:
    ExonID <- match(FinalSeq[n,12], SG$Edges$From)
    ExonSeq <- as.character(getSeq(Hsapiens,SG$Edges$Chr[ExonID], as.numeric(SG$Edges$Start[ExonID]), as.numeric(SG$Edges$End[ExonID])))
    ExonSeq <- unlist(strsplit(ExonSeq,""))
    seqfin <- paste(ExonSeq[0:(FinalSeq[n,8]-1)],collapse = "") 
  }
  # Now we can build the sequence seqinicio, seqfinal and the exons in Exons
  # As we dont put probes in juntions we put Ns between sequence of different exons:
  seq=paste(seq,seqinicio,sep="")
  for (i in 1:length(Exons)){
    # Get the ID for each exon:
    ExonID <- match(Exons[i], SG$Edges$From)
    # Get de sequence of that exon
    if (is.na(ExonID)==FALSE){
      ExonSeq <- as.character(getSeq(Hsapiens,SG$Edges$Chr[ExonID], as.numeric(SG$Edges$Start[ExonID]), as.numeric(SG$Edges$End[ExonID])))
      # Build input sequence for Primer3
      seq <- paste(seq,paste(rep("N",21),collapse = ""),ExonSeq, sep="")
    }else {
      ExonSeq <- ""   
    }
  }
  seq=paste(seq,paste(rep("N",21),collapse = ""),seqfin,sep="")
  
}

#' @rdname InternalFunctions
findPotencialExons <- function(D, namesPath,
                               maxLength, SG,
                               minexonlength)
  {
  
  # TODO: using the trick of the QP find out the nodes (in fact, the edges) 
  # that have a greater than or equal flux to the interrogated ones.
  
  # In order to do that, it is necessary to get the incidence matrix of the graph.
  # This code makes the trick
  # library(intergraph)
  # library(network)
  # as.matrix(asNetwork(wg),matrix.type="incidence")
  
  # The library to perform the QP is LowRankQP
  
  exonLengths <- diag(D[c(TRUE,FALSE),c(FALSE,TRUE)])
  names(exonLengths) <- rownames(D[c(TRUE,FALSE),c(FALSE,TRUE)])
  longExons <- names(which(exonLengths > minexonlength))
  fullsignalExons <- getExonsFullSignal(namesPath, SG)
  
  Reverse <- names(which(D[sort.exons(namesPath,decreasing = TRUE)[1],]<maxLength))
  Reverse <- union(fullExons(namesPath), Reverse)
  Reverse <- union(Reverse,includeaexons(Reverse))
  Reverse <- intersect(Reverse, longExons)
  Reverse <- intersect(Reverse, fullsignalExons)
  Reverse <- sort.exons(unique(Reverse))
  
  Forward <- names(which(D[,sort.exons(namesPath)[1]]<maxLength))
  Forward <- union(Forward,includeaexons(Forward))
  Forward <- union(fullExons(namesPath), Forward)
  Forward <- intersect(Forward, longExons)
  Forward <- intersect(Forward, fullsignalExons)
  Forward <- sort.exons(unique(Forward))
  
  return(list(Reverse=Reverse, Forward = Forward))
}
#' @rdname InternalFunctions
fullExons <- function(namesPath)
  {
  exonNumbers <- unlist(strsplit(namesPath,".", fixed=TRUE))[c(TRUE,FALSE)]
  if(length(which(duplicated(exonNumbers)))==0)
    return(NULL)
  exonNames <- paste(exonNumbers[which(duplicated(exonNumbers))],".a",sep="")
  return(exonNames)
}
#' @rdname InternalFunctions
includeaexons <- function(Forward)
  {
  exonNumbers <- unlist(strsplit(Forward,".", fixed=TRUE))[c(TRUE,FALSE)]
  exonNumbers <- unique(exonNumbers)
  exonNames <- paste(exonNumbers,".a",sep="")
  return(exonNames)
}

#' @rdname InternalFunctions
genreverse <- function(FinalInfo)
  {
  names <- names(FinalInfo)
  NewFinal <- FinalInfo
  # We aplly the reverse and complement to probes:
  rev <- function(x) {
    if (is.na(x))
    {return(NA)}
    else
    {return(as.character(reverseComplement(DNAString(x))))}}
  
  for (i in 13:15){
    x <- sapply(FinalInfo[,i],FUN = rev)
    NewFinal[,i] <- x
  }
  # We swap forward and reverse columns keeping names:
  NewFinal[,1:8] <- cbind(NewFinal[,c(3,4,1,2,7,8,5,6)])
  names(NewFinal) <- names
  return(NewFinal)
}


#' @rdname InternalFunctions
getDistanceseachPath<-function(Exon1,Exon2,
                               generaldata,
                               distinPrimers,
                               SG)
  {
  Event <- generaldata$Event
  nt<- generaldata$nt
  ntexons <- nt[SG$Edges$Type == "E"]
  Exon1 <- as.character(Exon1)
  Exon2 <- as.character(Exon2)
  namesP1 <- as.matrix(Event$P1[,1:2])
  namesP2 <- as.matrix(Event$P2[,1:2])
  Exon2mod <- str_replace(Exon2,"a","b")
  ways <- all_simple_paths2(generaldata$wg,Exon1,Exon2mod)
  listmatrixways <- list()
  # Convert each way to matrix of pairs of nodes which form juntions
  for(n in seq(1,length(ways)))
  {
    way <- names(unlist(ways[n]))
    modway <- c(way[1:length(way)-1],way[2:length(way)])
    matrixway<- matrix(modway,ncol = 2,byrow = FALSE)
    colnames(matrixway) <- colnames(as.matrix(Event$P1[,1:2])) 
    listmatrixways[[n]] <- matrixway
  }
  names(ntexons) <- SG$Edges$From[SG$Edges$Type == "E"]
  Info <- vector()
  
  for(i in seq(1,length(ways)))
  {
    # sum of exons that are in the path without primers
    entirexons <- names(unlist(ways[i]))[c(TRUE,FALSE)]
    if (entirexons[1]==Exon1){
      entirexons <- entirexons[-1]
    }
    if (entirexons[length(entirexons)]==Exon2){
      entirexons <- entirexons[-length(entirexons)]
    }
    a <- sum(ntexons[match(entirexons,names(ntexons))])
    Suma <- a + distinPrimers
    
    # There are 3 possibilities:
    # -Path1->1
    # -Path2->2
    # -No interrogate any path->0
    classPath <- "p0"
    
    if(sum(as.numeric(apply(listmatrixways[[i]],1, function(x)apply(namesP1,1,function(y) identical(y,x)))))>0){
      classPath <- "p1"
    } 
    if(sum(as.numeric(apply(listmatrixways[[i]],1, function(x)apply(namesP2,1,function(y) identical(y,x)))))>0){
      classPath <- "p2"
    }
    Info <- c(Info,Suma,classPath)
    
  }
  return(Info)
}

#' @rdname InternalFunctions
getDominants2<- function(PrimersTwo,Primers1,
                         commonForward,commonReverse,
                         namesRef,D,numberOfPaths,
                         nprimerstwo, ED,
                         wNpaths = 1000,
                         wP12inRef =1000)
  {
  
  # Use common references
  Primers1$Forwardfortwo  <- commonForward
  Primers1$Reversefortwo  <- commonReverse
  
  # Different measures
  NP1 <- numberOfPaths[Primers1$Forwardfortwo, Primers1$Reversefortwo,drop=FALSE]
  D1 <- D[Primers1$Forwardfortwo, Primers1$Reversefortwo,drop=FALSE]
  ED1 <- ED[Primers1$Forwardfortwo, Primers1$Reversefortwo,drop=FALSE]
  P1FinRef <- rownames(NP1) %in% namesRef
  names(P1FinRef) <- rownames(NP1)
  P1RinRef <- colnames(NP1) %in% namesRef
  names(P1RinRef) <- colnames(NP1)
  P1inRef <- outer(P1FinRef,P1RinRef,"+")
  # Matrix that measure the quality of the primers
  W1 <-  (NP1-2) * wNpaths + (2-P1inRef) * wP12inRef + D1 +ED1
  W1 <-  as.matrix(W1)
  size<-dim(W1)
  size<-size[1]*size[2]
  
  # Initializing the dataframe
  cols <- min(size,nprimerstwo)
  position<-matrix(ncol=2,nrow=cols)
  names<-matrix(ncol=5,nrow=cols)
  names[,5]<-rep(0,cols)
  value<-rep(0,cols)
  FINALvalue<-rep(0,cols)
  
  c=0
  for (n in 1:nprimerstwo) {
    
    a=which(W1 == min(W1), arr.ind = TRUE)
    if (n + c > nprimerstwo || n + c > size){
      break
    }else if (dim(a)[1]>1){
      for(m in 1:dim(a)[1]){
        position[n+c+m-1,]<-a[m,]
        value[n+c+m-1]=min(W1)*2
        W1[position[n+c+m-1,1],position[n+m+c-1,2]]=Inf
        names[n+c+m-1,1] <- rownames(W1)[position[n+c+m-1,1]]
        names[n+c+m-1,3] <- colnames(W1)[position[n+c+m-1,2]]
        if (n + c + m - 1>= nprimerstwo || n + c + m - 1 >= size){
          break
        }
      }
      
      c = c+dim(a)[1]-1
      
    }else{
      position[n+c,] <- a 
      value[n+c] = min(W1)*2
      W1[position[n+c,1],position[n,2]] = Inf
      names[n+c,1]<-rownames(W1)[position[n+c,1]]
      names[n+c,3]<- colnames(W1)[position[n+c,2]]
    }
    if (n + c >= nprimerstwo || n + c >= size){
      break
    }
    
  }
  position1 <- data.frame(names,value,FINALvalue)
  return(position1)
}

#' @rdname InternalFunctions
getDominantsFor <- function(Primers1,Primers2,
                            commonForward,namesRef,
                            D,numberOfPaths,Event,
                            ncommonForward,ED,
                            wNpaths = 1000, 
                            wP12inRef =1000)
  {
  
  # Use common references in the case of Forward, in Reverse just change the name
  
  Primers1$ForwardforcommForw <- Primers2$ForwardforcommForw <- commonForward
  Primers1$ReverseforcommForw <- Primers1$Reverse
  Primers2$ReverseforcommForw <- Primers2$Reverse
  
  # Adding exons in Paths to the PotencialPrimers in Reverse in the case of commonForward
  # Getting exons in path 1 and path 2
  
  exonsp1 <- unique( as.vector(as.matrix(Event$P1[,c(1,2)]))[grep("a",as.vector(as.matrix(Event$P1[,c(1,2)])))])
  exonsp2 <- unique( as.vector(as.matrix(Event$P2[,c(1,2)]))[grep("a",as.vector(as.matrix(Event$P2[,c(1,2)])))])
  
  # Adding those exons in Potencial Primers for Reverse
  Primers1$ReverseforcommForw <- unique(c(Primers1$ReverseforcommForw,exonsp1))
  Primers2$ReverseforcommForw <- unique(c(Primers2$ReverseforcommForw,exonsp2))
  
  # Different measures
  NP1 <- numberOfPaths[Primers1$ForwardforcommForw, Primers1$ReverseforcommForw,drop=FALSE]
  D1 <- D[Primers1$ForwardforcommForw, Primers1$ReverseforcommForw,drop=FALSE]
  ED1 <- ED[Primers1$ForwardforcommForw, Primers1$ReverseforcommForw,drop=FALSE]
  P1FinRef <- rownames(NP1) %in% namesRef
  names(P1FinRef) <- rownames(NP1)
  P1RinRef <- colnames(NP1) %in% namesRef
  names(P1RinRef) <- colnames(NP1)
  P1inRef <- outer(P1FinRef,P1RinRef,"+")
  
  NP2 <- numberOfPaths[Primers2$ForwardforcommForw, Primers2$ReverseforcommForw,drop=FALSE]
  D2 <- D[Primers2$ForwardforcommForw, Primers2$ReverseforcommForw,drop=FALSE]
  ED2 <- ED[Primers2$ForwardforcommForw, Primers2$ReverseforcommForw,drop=FALSE]
  P2FinRef <- rownames(NP2) %in% namesRef
  names(P2FinRef) <- rownames(NP2)
  P2RinRef <- colnames(NP2) %in% namesRef
  names(P2RinRef) <- colnames(NP2)
  P2inRef <- outer(P2FinRef,P2RinRef,"+")
  
  # There are three matrices that measure the quality of the primers.
  W1 <-  (NP1-1) * wNpaths + (2-P1inRef) * wP12inRef + D1 +ED1
  W2 <-  (NP2-1) * wNpaths + (2-P2inRef) * wP12inRef + D2 +ED2
  
  sizeW1<-dim(W1)
  numberW1<-sizeW1[1]*sizeW1[2]
  sizeW2<-dim(W2)
  numberW2<-sizeW2[1]*sizeW2[2]
  
  Sum <-matrix(nrow=length(commonForward),ncol=sizeW1[2]*sizeW2[2])
  
  positionW1<-matrix(ncol=2,nrow=ncommonForward)
  valueW1=rep(0,ncommonForward)
  namesW1<-matrix(ncol=2,nrow=ncommonForward)
  
  positionW2<-matrix(ncol=2,nrow=ncommonForward)
  valueW2=rep(0,ncommonForward)
  namesW2<-matrix(ncol=2,nrow=ncommonForward)
  
  positioncommonForward<-matrix(ncol=2,nrow=ncommonForward)
  value <- rep(0,ncommonForward)
  FINALvalue <- rep(0,ncommonForward)
  namescommonForward<-matrix(ncol=5,nrow=ncommonForward)
  namescommonForward[,5]<-rep(1,ncommonForward)
  
  c=0
  d=0
  e=0
  
  for (j in 1:sizeW1[2]) 
  {
    Sum[,(j-1)*sizeW2[2]+1:sizeW2[2]]=  as.matrix(as.numeric(W1[commonForward,j ]) + W2[commonForward,,drop=FALSE])
  }
  
  for (n in 1:ncommonForward) 
  {
    c=which(Sum == min(Sum), arr.ind = TRUE)
    if (n+e>ncommonForward )
    {
      break
    }else if (dim(c)[1]>1){
      for(m in 1:dim(c)[1])
      {
        positioncommonForward[n+e+m-1,]<-c[m,]
        value[n+e+m-1]=min(Sum)
        Sum[positioncommonForward[n+e+m-1,1],positioncommonForward[n+e+m-1,2]]=Inf
        namescommonForward[n+e+m-1,1]<-rownames(W1)[positioncommonForward[n+m-1,1]]
        namescommonForward[n+e+m-1,3]<-colnames(W1)[ceiling(positioncommonForward[n+e+m-1,2]/sizeW2[2])]
        namescommonForward[n+e+m-1,4]<-colnames(W2)[positioncommonForward[n+e+m-1,2]-ceiling(positioncommonForward[n+e+m-1,2]/sizeW2[2])*sizeW2[2]+sizeW2[2]]
        if (n+e+m-1>=ncommonForward ){
          break
        }
      }
      e=e+dim(c)[1]-1
      
    }else{
      positioncommonForward[n+e,] <- c
      value[n+e] = min(Sum)
      Sum[positioncommonForward[n+e,1],positioncommonForward[n+e,2]]=Inf
      namescommonForward[n+e,1] <- rownames(W1)[positioncommonForward[n+e,1]]
      namescommonForward[n+e,3] <- colnames(W1)[ceiling(positioncommonForward[n+e,2]/sizeW2[2])]
      namescommonForward[n+e,4] <- colnames(W2)[positioncommonForward[n+e,2]-ceiling(positioncommonForward[n+e,2]/sizeW2[2])*sizeW2[2]+sizeW2[2]]
    }
  }
  position<-data.frame(namescommonForward,value, FINALvalue)
  return(position)
}

#' @rdname InternalFunctions
getDominantsRev<- function(Primers1,Primers2,
                           commonReverse,namesRef,
                           D,numberOfPaths,Event,
                           ncommonReverse,ED,
                           wNpaths = 1000,
                           wP12inRef =1000)
  {
  # Use common references in the case of Reverse, in Forward just change the name
  
  Primers1$ForwardforcommRevw <- Primers1$Forward
  Primers2$ForwardforcommRevw <- Primers2$Forward
  Primers1$ReverseforcommRevw <- Primers2$ReverseforcommRevw<- commonReverse
  
  # Adding exons in Paths to the PotencialPrimers in Forward in the case of commonReverse
  # Getting exons in path 1 and path 2
  exonsp1 <- unique( as.vector(as.matrix(Event$P1[,c(1,2)]))[grep("a",as.vector(as.matrix(Event$P1[,c(1,2)])))])
  exonsp2 <- unique( as.vector(as.matrix(Event$P2[,c(1,2)]))[grep("a",as.vector(as.matrix(Event$P2[,c(1,2)])))])
  # Adding those exons in Potencial Primers for Reverse
  Primers1$ReverseforcommRev <- unique(c(Primers1$ReverseforcommRev,exonsp1))
  Primers2$ReverseforcommRev <- unique(c(Primers2$ReverseforcommRev,exonsp2))
  
  # Different measures
  NP1 <- numberOfPaths[Primers1$ForwardforcommRevw, Primers1$ReverseforcommRevw,drop=FALSE]
  D1 <- D[Primers1$ForwardforcommRevw, Primers1$ReverseforcommRevw,drop=FALSE]
  ED1 <- ED[Primers1$ForwardforcommRevw, Primers1$ReverseforcommRevw,drop=FALSE]
  P1FinRef <- rownames(NP1) %in% namesRef
  names(P1FinRef) <- rownames(NP1)
  P1RinRef <- colnames(NP1) %in% namesRef
  names(P1RinRef) <- colnames(NP1)
  P1inRef <- outer(P1FinRef,P1RinRef,"+")
  
  NP2 <- numberOfPaths[Primers2$ForwardforcommRevw, Primers2$ReverseforcommRevw,drop=FALSE]
  D2 <- D[Primers2$ForwardforcommRevw, Primers2$ReverseforcommRevw,drop=FALSE]
  ED2 <- ED[Primers2$ForwardforcommRevw, Primers2$ReverseforcommRevw,drop=FALSE]
  
  P2FinRef <- rownames(NP2) %in% namesRef
  names(P2FinRef) <- rownames(NP2)
  P2RinRef <- colnames(NP2) %in% namesRef
  names(P2RinRef) <- colnames(NP2)
  P2inRef <- outer(P2FinRef,P2RinRef,"+")
  
  # There are three matrices that measure the quality of the primers.
  W1 <-  (NP1-1) * wNpaths + (2-P1inRef) * wP12inRef + D1 +ED1
  W2 <-  (NP2-1) * wNpaths + (2-P2inRef) * wP12inRef + D2 +ED2
  
  
  W1<-t(W1)
  W2<-t(W2)
  sizeW1<-dim(W1)
  numberW1<-sizeW1[1]*sizeW1[2]
  sizeW2<-dim(W2)
  numberW2<-sizeW2[1]*sizeW2[2]
  
  Sum  <-matrix(nrow=length(commonReverse),ncol=sizeW1[2]*sizeW2[2])
  
  positionW1<-matrix(ncol=2,nrow=ncommonReverse)
  valueW1=rep(0,ncommonReverse)
  namesW1<-matrix(ncol=2,nrow=ncommonReverse)
  
  positionW2<-matrix(ncol=2,nrow=ncommonReverse)
  valueW2=rep(0,ncommonReverse)
  namesW2<-matrix(ncol=2,nrow=ncommonReverse)
  
  positioncommonReverse<-matrix(ncol=2,nrow=ncommonReverse)
  value=rep(0,ncommonReverse)
  FINALvalue <- rep(0,ncommonReverse)
  namescommonReverse<-matrix(ncol=5,nrow=ncommonReverse)
  namescommonReverse[,5]<-rep(2,ncommonReverse)
  
  c=0
  d=0
  e=0
  
  for (j in 1:sizeW1[2]) {
    Sum[,(j-1)*sizeW2[2]+1:sizeW2[2]]= as.numeric(W1[commonReverse,j]) + 
      as.matrix(W2[commonReverse,,drop= FALSE])
  }
  for (n in 1:ncommonReverse) {
    c=which(Sum == min(Sum), arr.ind = TRUE)
    if (n+e>ncommonReverse ){
      break
      
    }else if (dim(c)[1]>1){
      for(m in 1:dim(c)[1]){
        positioncommonReverse[n+e+m-1,]<-c[m,]
        value[n+e+m-1]=min(Sum)
        Sum[positioncommonReverse[n+e+m-1,1],positioncommonReverse[n+e+m-1,2]]=Inf
        namescommonReverse[n+e+m-1,3]<-rownames(W1)[positioncommonReverse[n+m-1,1]]
        namescommonReverse[n+e+m-1,1]<-colnames(W1)[ceiling(positioncommonReverse[n+e+m-1,2]/sizeW2[2])]
        namescommonReverse[n+e+m-1,2]<-colnames(W2)[positioncommonReverse[n+e+m-1,2]-ceiling(positioncommonReverse[n+e+m-1,2]/sizeW2[2])*sizeW2[2]+sizeW2[2]]
        
        if (n+e+m-1>=ncommonReverse ){
          break
        }
      }
      e = e + dim(c)[1]-1
      
    }else{
      positioncommonReverse[n+e,] <- c
      value[n+e] = min(Sum)
      Sum[positioncommonReverse[n + e,1],positioncommonReverse[n + e,2]] = Inf
      namescommonReverse[n + e,3] <- rownames(W1)[positioncommonReverse[n + e,1]]
      namescommonReverse[n + e,1] <- colnames(W1)[ceiling(positioncommonReverse[n+e,2]/sizeW2[2])]
      namescommonReverse[n + e,2] <- colnames(W2)[positioncommonReverse[n + e,2]-ceiling(positioncommonReverse[n+e,2]/sizeW2[2])*sizeW2[2]+sizeW2[2]]
      
    }
    
  }
  position <- data.frame(namescommonReverse,value,FINALvalue)
  
  return(position)
}

#' @rdname InternalFunctions
getExonsFullSignal <- function(namesPath, SG)
  {
  # This function finds out the exons whose expression is larger than the 
  # concentration path under study for any isoform expression distribution.
  
  # It is solved by using quadratic programming:
  # min |v|^2
  # s.t.
  # A v_(-i) = 0
  # v_i = 1
  
  # The solution of this optimization problem will be a distribution of fluxes in which, the fluxes that are 
  # bigger or equal than the flux under study will be close to one.
  
  # Input: SG, splicing graph (igraph object).
  # namesPath: set of nodes under study in the path.
  
  # Get the incidence matrix that corresponds to the splicing graph
  A <- SG$Incidence
  # Remove Start and End nodes.
  A <- A[-match(c("S","E"),rownames(A)),]
  
  # Find the flow that must be one
  # Index <- which(A[namesPath[grep("b",namesPath)],,drop=FALSE]==1, arr.ind = TRUE)[1,2]
  Index <- as.vector(which(colSums(abs(A[namesPath,]))>=2))[1]
  # Note: I look only in the "b" nodes. The corresponding row of the incidence matrix must 
  # have one 1 (and only one) for these nodes. I keep the first one (any of them would 
  # be OK by construction).
  A <- rbind(A,0)
  A[nrow(A),Index] <- 1
  b <- rep(0, nrow(A))
  b[nrow(A)] <- 1
  b<- c(b,rep(0,dim(A)[2]))
  A <- rbind(A,diag(dim(A)[2])*1e-5)
  # ---
  Output <- nnls(A,b)
  flows <- which(abs(Output$x - 1) < 1e-5)
  nodes <- rownames(which(abs(A[1:(nrow(A)-1),flows])==1, arr.ind = TRUE))
  # ---
  # capture.output(Sol <- LowRankQP(diag(ncol(A)),rep(0,ncol(A)),A,b,
  #                  uvec=rep(1e4,ncol(A)),method="CHOL", verbose=FALSE))
  # # Get the flows (and the corresponding nodes)
  # flows <- which(abs(Sol$alpha - 1) < 1e-5)
  # nodes <- rownames(which(abs(A[1:(nrow(A)-1),flows])==1, arr.ind = TRUE))
  # ---
  return(unique(nodes))
}

#' @rdname InternalFunctions
getFinalExons<- function(generaldata,maxLength,
                         nPrimerstwo,ncommonForward,
                         ncommonReverse,nExons,
                         minsep,wminsep,
                         valuethreePpenalty,
                         minexonlength)
  {
  # gen info
  SG<-generaldata[[1]]
  # Splicing Event
  Event<-generaldata[[2]]
  nt<-generaldata[[3]]
  # weighted splicing graph
  wg<-generaldata[[4]]
  # Distance matrix
  D<-generaldata[[5]]
  # Number of paths connecting two nodes in the graph
  numberOfPaths<-generaldata[[6]]
  # exons length penalty
  ED <- generaldata[[7]]
  # Ref nodes
  exonsPathsandRef <- generaldata[[8]]
  namesRef <- exonsPathsandRef$Reference
  
  
  # Potential primers based only on distances to Paths 1 and 2
  
  namesP1 <- unique(as.vector(as.matrix(Event$P1[,1:2])))
  namesP1 <- namesP1[namesP1 %in% colnames(D)]
  Primers1 <- findPotencialExons(D, namesP1, maxLength, SG=SG,minexonlength)
  Primers1$Reverse  <- Primers1$Reverse[grep("a",Primers1$Reverse)]
  Primers1$Forward  <- Primers1$Forward[grep("a",Primers1$Forward)]
  
  
  namesP2 <- unique(as.vector(as.matrix(Event$P2[,1:2])))
  namesP2 <- namesP2[namesP2 %in% colnames(D)]
  Primers2 <- findPotencialExons(D, namesP2, maxLength, SG=SG,minexonlength)
  Primers2$Reverse  <- Primers2$Reverse[grep("a",Primers2$Reverse)]
  Primers2$Forward  <- Primers2$Forward[grep("a",Primers2$Forward)]
  
  # Check for three primers: both forward and/or both reverse must have nodes that overlap
  commonForward <- intersect(Primers1$Forward, Primers2$Forward)
  commonReverse <- intersect(Primers1$Reverse, Primers2$Reverse)
  
  # Initializing for using rbind later
  PrimersTwo<- PrimersFor<- PrimersRev<-matrix(ncol=7,nrow=0)
  
  
  if(length(Primers1$Forward)==0 || length(Primers1$Reverse)==0 || length(Primers2$Forward)==0 || length(Primers2$Reverse)==0  )
  {
    FinalExons<-"Not possible to place primers due to the structure of the Event."
  }else{
    
    # If two primers are sufficient there are three possible cases:
    # 1) Use only two primers
    # 2) Use common FW primer and two Reverse primers
    # 3) Use commong reverse primer and two FW primers
    
    # 1) Use only two primers
    # Check if two primers could be sufficient
    
    if ((length(commonReverse) >0) &  (length(commonForward) >0) & (nPrimerstwo!=0)) 
    {
      PrimersTwo<-getDominants2(PrimersTwo,Primers1,commonForward,commonReverse,namesRef,D,numberOfPaths,nPrimerstwo,ED)
    }
    
    # 2) Using common Forward
    # Check if it is possible three primers in common Forward
    
    if (length(commonForward) >0 & (ncommonForward!=0)) 
    {
      PrimersFor<-getDominantsFor(Primers1,Primers2,commonForward,namesRef,D,numberOfPaths,Event,ncommonForward,ED)
      # We get in order to delete the identical
      # With PrimersFor
      PrimersFor[,3] <- as.character(PrimersFor[,3])
      PrimersFor[,4] <- as.character(PrimersFor[,4])
      PrimersFor[which(t(apply(PrimersFor[,3:4],1,FUN=order))[,1]==2),3:4] <- PrimersFor[which(t(apply(PrimersFor[,3:4],1,FUN=order))[,1]==2),4:3]
      PrimersFor[,3] <- as.factor(PrimersFor[,3])
      PrimersFor[,4] <- as.factor(PrimersFor[,4])
    }
    
    # 3) Using common Reverse
    # Check if it is possible three primers in common Reverse
    if (length(commonReverse) >0 & (ncommonReverse!=0)) 
    {
      PrimersRev<-getDominantsRev(Primers1,Primers2,commonReverse,namesRef,D,numberOfPaths,Event,ncommonReverse,ED)
      # We get in order to delete the identical
      # With PrimersRev
      PrimersRev[,1] <- as.character(PrimersRev[,1])
      PrimersRev[,2] <- as.character(PrimersRev[,2])
      PrimersRev[which(t(apply(PrimersRev[,1:2],1,FUN=order))[,1]==2),1:2] <- PrimersRev[which(t(apply(PrimersRev[,1:2],1,FUN=order))[,1]==2),2:1]
      PrimersRev[,1] <- as.factor(PrimersRev[,1])
      PrimersRev[,2] <- as.factor(PrimersRev[,2])
    }
    
    colnames(PrimersTwo)<-colnames(PrimersFor)<-colnames(PrimersRev)<-c("For1Exon","For2Exon","Rev1Exon","Rev2Exon","3Primers","value","Finalvalue")
    
    # Getting results together
    Dominants<- rbind(PrimersTwo,PrimersFor,PrimersRev)
    Dominants <- Dominants[rownames(unique(Dominants[,1:4])),]
    BadOnes <- which(as.character(Dominants$For1Exon)==as.character(Dominants$For2Exon))
    if(length(BadOnes)>0) Dominants <- Dominants[-BadOnes,]
    BadOnes <- which(as.character(Dominants$Rev1Exon)==as.character(Dominants$Rev2Exon))
    if(length(BadOnes)>0) Dominants <- Dominants[-BadOnes,]
    # Order quality of Dominants with differences of distances and 3 or 2 primers
    
    if (nrow(Dominants)>0){
      FinalExons<-getrankexons(SG,Dominants,nt,wg,nExons,minsep,wminsep,valuethreePpenalty,D)
    }else{
      FinalExons<-"Not possible to place primers due to the structure of the Event."
    }
  }
  return(FinalExons)
}


#' @rdname InternalFunctions
getgeneraldata <- function(SG, Event,shortdistpenalty)
  {
  
  # Get adjacency matrix
  Adj <- SG$Adjacency
  
  # Create a weighted adjacency matrix
  WAdj <- Adj
  nt <- abs(as.numeric(SG$Edges$End) - as.numeric(SG$Edges$Start))+1
  nt[SG$Edges$Type == "J"] <- .001
  WAdj[cbind(match(SG$Edges$From,rownames(Adj)),match(SG$Edges$To,colnames(Adj)))] <- nt
  
  
  # Get weighted splicing graph
  wg <- graph_from_adjacency_matrix(WAdj, weighted = TRUE)
  
  # Get Distance matrix
  D <- distances(wg, mode = "out")
  Keep <- !colnames(D) %in% c("S","E")
  D <- D[Keep, Keep]
  
  # Get Distance of exons
  exonlist <- SG$Edges[which(SG$Edges[["Type"]]=="E"),]
  exonlength <-abs(as.numeric(exonlist$End) - as.numeric(exonlist$Start))+1
  names(exonlength)<- exonlist$From
  # Calculeta Penaltydistance for primers in exons
  shortexonpenalty<- (shortdistpenalty/0.24)*exp(-0.04 *exonlength)
  names(shortexonpenalty)<- exonlist$From
  nexons<-length(shortexonpenalty)
  # Create matrix ED with penalty of dexondistances
  EDcol <- matrix(data = rep(shortexonpenalty,nexons),nrow = nexons,ncol = nexons,byrow = TRUE )
  EDrow <- t(EDcol)
  ED <- EDrow + EDcol
  colnames(ED) <- rownames(ED) <- exonlist$From
  
  # Number of paths connecting two nodes in the graph
  numberOfPaths <- solve(Diagonal(ncol(Adj))-Adj)
  colnames(numberOfPaths) <- rownames(numberOfPaths) <- colnames(Adj)
  Keep <- !colnames(numberOfPaths) %in% c("S","E")
  numberOfPaths <- numberOfPaths[Keep, Keep]
  
  # Calculate vectors with exons of P1, P2 and reference:
  Path1 <- Event$P1[which(Event$P1[,"Type"]=="E"),"From"]
  Path2 <- Event$P2[which(Event$P2[,"Type"]=="E"),"From"]
  Reference <- Event$Ref[which(Event$Ref[,"Type"]=="E"),"From"]
  exonsPathsandRef <- list(Path1,Path2,Reference)
  names(exonsPathsandRef) <- c("Path1","Path2","Reference")
  
  datalist <- list(SG,Event,nt, wg, D, numberOfPaths,ED,exonsPathsandRef)
  names(datalist)<- c("SG","Event","nt", "wg", "D", "numberOfPaths","ED","exonsPathsandRef")
  return (datalist)
}


#' @rdname InternalFunctions
getrankexons<- function(SG,Dominants,nt,
                        wg,items,minsep,
                        wminsep,valuethreePpenalty,
                        D)
  {
  #  puntuation using value,3 or 2 primers and difdistances.
  # neccesary to calculate difdistancespenalty
  nDominants<-dim(Dominants)[1]
  items <- min(items,nDominants)
  Dominants[,7]<-rep(0,nDominants)
  
  names(Dominants)[7]<-"FINALvalue"
  # calculate exons distances
  ntexons <- nt[SG$Edges$Type == "E"]
  names(ntexons) <- SG$Edges$From[SG$Edges$Type == "E"]
  
  for (k in 1:nDominants) {
    #  two possible cases:
    # 1) Use only two primers
    # 2) Use 3 primers
    
    # Calculate penalty with 3 primers
    if (Dominants[k,5]!=0){
      # There are two cases:
      #   2.1)commonForward
      #   2.2)commonReverse
      if (Dominants[k,5]==1){
        commonPrimer <- as.vector(Dominants[k,1]) 
        path1Primer <- as.vector(Dominants[k,3])
        path2Primer <- as.vector(Dominants[k,4])
        path1firstdist <- D[commonPrimer,path1Primer]
        path2firstdist <- D[commonPrimer,path2Primer]
        path1secdist   <- path1firstdist + ntexons[Dominants[k,3]]
        path2secdist   <- path2firstdist + ntexons[Dominants[k,4]]
      }
      
      if (Dominants[k,5]==2){
        commonPrimer <- as.vector(Dominants[k,3])
        path1Primer <- as.vector(Dominants[k,1])
        path2Primer <- as.vector(Dominants[k,2])
        path1firstdist <- D[path1Primer,commonPrimer]
        path2firstdist <- D[path2Primer,commonPrimer]
        path1secdist   <- path1firstdist + ntexons[Dominants[k,1]]
        path2secdist   <- path2firstdist + ntexons[Dominants[k,2]]
      }
      # Estimation with Uniform distribution of the mean of difdistances:
      dar <- runif(1000,path1firstdist,path1secdist)
      dbr <- runif(1000,path2firstdist,path2secdist)
      DeltaABr <- abs(dar - dbr)
      mindist <- mean(DeltaABr)
    }
    
    # Calculate penalty with 2 primers
    if (Dominants[k,5]==0){
      a <- as.vector(Dominants[k,1])
      b <- as.vector(Dominants[k,3])
      allPaths <- all_simple_paths2(wg,from = a , to = b)
      if (length(allPaths)==0)
        distances <- 0
      else
        distances<-sapply(allPaths, FUN = function(x) {return(sum(ntexons[as_ids(x)[c(TRUE,FALSE)]]))})
      distances<-sort(distances)
      difdistances<-diff(distances)
      mindist <- min(difdistances)
    }
    
    # calculate difdistancespenalty with distances if they are lower than minsparation between distance
    difdistancespenalty<-0
    if (mindist < minsep) {
      difdistancespenalty <- wminsep/1000*(mindist-minsep)^2
    }
    
    threePrimerspenalty <- 0
    if (as.numeric(as.vector(Dominants[k,5]))!= 0){
      threePrimerspenalty <- 1
    }
    
    # type of primers either two(0) or three(1) is in col=5
    PreFinvalue <- Dominants[k,6] + threePrimerspenalty * valuethreePpenalty + difdistancespenalty
    Dominants[k,7] <- PreFinvalue
    # PreFinal value is in col=7
  }
  # Order and return of as many items as user wants:
  X <- order(Dominants[,7])
  Dominants <- Dominants[X,]
  FinalExons<-Dominants[1:items,]
  return (FinalExons)
}


#' @rdname InternalFunctions
getranksequence<- function(taqman,Fdata,maxLength,minsep
                           ,wminsep,valuethreePpenalty
                           ,wnpaths,qualityfilter)
  {
  for (i in 1:dim(Fdata)[1]){
    ValueForSequence <- 0
    # Calculate penalty for primers:
    # ExtractDistances:
    DistPath1 <- as.character(Fdata[i,15])
    DistPath1 <- as.integer(unlist(strsplit(DistPath1," - ")))
    DistPath2 <- as.character(Fdata[i,16])
    DistPath2 <- as.integer(unlist(strsplit(DistPath2," - ")))
    DistNoPath <- as.character(Fdata[i,17])
    DistNoPath <- as.integer(unlist(strsplit(DistNoPath," - ")))
    
    # First we apply common penalties:
    # Penalty for number of primers:
    if(Fdata[i,13]!= 0 ){
      ValueForSequence <- ValueForSequence +  valuethreePpenalty
    }
    # Penlty for number of paths:
    NPaths <- length(DistPath1)+length(DistPath2) + length(DistNoPath)
    ValueForSequence <- ValueForSequence + wnpaths * NPaths
    
    # Penlty for long paths: penalty for each long path and if all paths are long we
    # apply another penalty:
    numberlongpaths <- sum( c(DistPath1, DistPath2) > maxLength)
    length <- max(min(DistPath1),min(DistPath2))
    ValueForSequence <- ValueForSequence + length
    ValueForSequence <- ValueForSequence + numberlongpaths * wnpaths* 1.5
    if (length > maxLength){
      ValueForSequence <- ValueForSequence + 10000
    }
    
    
    # We apply specific penalties for taqman and conventional PCR:
    if (taqman == 1){
      longer <- max(c(min(DistPath1),min(DistPath2)))
      ValueForSequence <- ValueForSequence + longer 
      
    }else{
      # Penalty of diference of distances between all paths that are shorter than maxLenght
      distances <- c(DistPath1,DistPath2)
      distances[distances<maxLength]
      distances<-sort(distances)
      difdistances<-diff(distances)
      mindist <- min(difdistances)
      if (mindist<minsep) {
        difdistancespenalty <- wminsep/1000*(mindist-minsep)^2
        ValueForSequence <- ValueForSequence + difdistancespenalty
      }
    }
    Fdata[i,14] <- ValueForSequence
  }
  # We order the data with the penalties
  RankedFdata<- Fdata[order(Fdata[,14]),]
  # We take only values lower than qualityfilter, if there are not good primers we take only first 3 primers:
  RankedFdata <- unique(rbind(RankedFdata[1:3,],RankedFdata[RankedFdata[,14]<qualityfilter,]))
  RankedFdata <- RankedFdata[is.na(RankedFdata)[,1]==FALSE,]
  return(RankedFdata)
}


#' @rdname InternalFunctions
PrimerSequenceCommonFor <-function(FinalExons,SG,
                                   generaldata,n,
                                   thermo.param,
                                   Primer3Path,
                                   settings)
  {
  Fdata1 <- data.frame()
  # Case1: Find sequence in the first Exon and once we have set that sequence we look for the sequence in the second exon.
  For1Exon <- FinalExons[n,1]
  Rev1Exon <- FinalExons[n,3]
  Rev2Exon <- FinalExons[n,4]
  # Get the ID for each exon where Primers are placed
  ExonF1 <- match(as.character(FinalExons[n,1]), SG$Edges$From)
  ExonR1 <- match(as.character(FinalExons[n,3]), SG$Edges$From)
  ExonR2 <- match(as.character(FinalExons[n,4]), SG$Edges$From)
  # Get de sequence of each exon
  Hsapienshg38  <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  FExonSeq  <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonF1], as.numeric(SG$Edges$Start[ExonF1]), as.numeric(SG$Edges$End[ExonF1])))
  RExonSeq1 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonR1], as.numeric(SG$Edges$Start[ExonR1]), as.numeric(SG$Edges$End[ExonR1])))
  RExonSeq2 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonR2], as.numeric(SG$Edges$Start[ExonR2]), as.numeric(SG$Edges$End[ExonR2])))
  # Call primer3 with two the sequence of two of the exons:
  minlength <- max(c(str_length(FExonSeq),str_length(RExonSeq1)))+1
  seq1 <- paste(as.character(FExonSeq), paste(rep("N",minlength),collapse = "") ,as.character(RExonSeq1),sep="")
  maxlength <- str_length(seq1) 
  
  p1 <- callPrimer3(seq1, size_range =sprintf("%0.0f-%0.0f", minlength, maxlength),
                    name = "Primer1", Primer3Path = Primer3Path,
                    thermo.param = thermo.param,
                    sequence_target = 110,
                    settings = settings)
  # When a sequence is finded in first exon we look for de sequence in the other exon:
  if(is.null(dim(p1))==FALSE){
    for (s in 1: dim(p1)[1]) {
      pr=p1[s,2]
      minlength2 <- max(c(str_length(FExonSeq),str_length(RExonSeq2)))+1
      seq2 <- paste(as.character(FExonSeq), paste(rep("N",minlength2),collapse = "") ,as.character(RExonSeq2),sep="")
      maxlength2 <- str_length(seq2)
      p2 <- callPrimer3(seq2,threeprimers = TRUE,pr=pr,reverse = FALSE , size_range = sprintf("%0.0f-%0.0f", minlength2, maxlength2),
                        name = "Primer1", Primer3Path = Primer3Path,
                        thermo.param = thermo.param,
                        sequence_target = 110+ minlength2,
                        settings = settings)
      # If we have all Sequences needed we build Output binding with new Sequence info and knowed Exon info
      if(is.null(dim(p2))==FALSE){
        for (d in 1: dim(p2)[1]) {
          # We calculate primers positions in exons
          For1Seq <- p1[s,2]
          Rev1Seq <- p1[s,3]
          Rev2Seq <- p2[d,3]
          LastPosFor1 <- p1[s,6]+ p1[s,7] - 1
          LastPosFor2 <- NA
          FirstPosRev1 <- p1[s,8]-str_length(FExonSeq)-minlength-p1[s,9]+1
          FirstPosRev2 <- p2[d,8]-str_length(FExonSeq)-minlength2-p2[d,9]+1
          distinPrimers1 <- p1$PRIMER_PAIR_PRODUCT_SIZE[s]- minlength
          InfoC1 <- getDistanceseachPath(For1Exon,Rev1Exon,generaldata,distinPrimers1,SG)
          distinPrimers2 <- p2$PRIMER_PAIR_PRODUCT_SIZE[d]- minlength2
          InfoC2 <- getDistanceseachPath(For1Exon,Rev2Exon,generaldata,distinPrimers2,SG)
          Info <- c(InfoC1,InfoC2)
          # Distances from any path:
          DistancesP0 <- paste(Info[which(Info=="p0")-1][order(as.integer(Info[which(Info=="p0")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path1:
          DistancesP1 <- paste(Info[which(Info=="p1")-1][order(as.integer(Info[which(Info=="p1")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path2
          DistancesP2 <- paste(Info[which(Info=="p2")-1][order(as.integer(Info[which(Info=="p2")-1]),decreasing = FALSE)],collapse = " - ")
          Distances <- data.frame(DistancesP1,DistancesP2,DistancesP0)
          names(Distances) <-c("DistPath1","DistPath2","DistNoPath")
          FExons <- data.frame(FinalExons[n,1],FinalExons[n,2],FinalExons[n,3],FinalExons[n,4],FinalExons[n,5],FinalExons[n,7])
          colnames(FExons) <- colnames(FinalExons)[-6]
          PrSeq<-data.frame(For1Seq,NA,Rev1Seq,Rev2Seq,LastPosFor1,LastPosFor2,FirstPosRev1,FirstPosRev2)
          colnames(PrSeq)<-c("For1Seq","For2Seq","Rev1Seq","Rev2Seq","LastPosFor1","LastPosFor2","FirstPosRev1","FirstPosRev2")
          Fdata1<-rbind(Fdata1,cbind(PrSeq , FExons,Distances))
        }
      }
    }
  }
  # 2ยบ Case: Find sequence in the second Exon and once we have set that sequence we look for the sequence in the first exon.
  # For that reason we swap exons:
  a <-  FinalExons[,3]
  FinalExons[,3] <- FinalExons[,4]
  FinalExons[,4] <- a
  
  For1Exon <- FinalExons[n,1]
  Rev1Exon <- FinalExons[n,3]
  Rev2Exon <- FinalExons[n,4]
  # Get the ID for each exon where Primers are placed
  ExonF1 <- match(as.character(FinalExons[n,1]), SG$Edges$From)
  ExonR1 <- match(as.character(FinalExons[n,3]), SG$Edges$From)
  ExonR2 <- match(as.character(FinalExons[n,4]), SG$Edges$From)
  # Get de sequence of each exon
  FExonSeq  <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonF1], as.numeric(SG$Edges$Start[ExonF1]), as.numeric(SG$Edges$End[ExonF1])))
  RExonSeq1 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonR1], as.numeric(SG$Edges$Start[ExonR1]), as.numeric(SG$Edges$End[ExonR1])))
  RExonSeq2 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonR2], as.numeric(SG$Edges$Start[ExonR2]), as.numeric(SG$Edges$End[ExonR2])))
  # Call primer3 with two the sequence of two of the exons:
  minlength <- max(c(str_length(FExonSeq),str_length(RExonSeq1)))+1
  seq1 <- paste(as.character(FExonSeq), paste(rep("N",minlength),collapse = "") ,as.character(RExonSeq1),sep="")
  maxlength <- str_length(seq1) 
  
  p1 <- callPrimer3(seq1, size_range =sprintf("%0.0f-%0.0f", minlength, maxlength),
                    name = "Primer1", Primer3Path = Primer3Path,
                    thermo.param = thermo.param,
                    sequence_target = 110,
                    settings = settings)
  # When a sequence is finded in first exon we look for de sequence in the other exon:
  if(is.null(dim(p1))==FALSE){
    for (s in 1: dim(p1)[1]) {
      pr=p1[s,2]
      minlength2 <- max(c(str_length(FExonSeq),str_length(RExonSeq2)))+1
      seq2 <- paste(as.character(FExonSeq), paste(rep("N",minlength2),collapse = "") ,as.character(RExonSeq2),sep="")
      maxlength2 <- str_length(seq2)
      p2 <- callPrimer3(seq2,threeprimers = TRUE,pr=pr,reverse = FALSE , size_range = sprintf("%0.0f-%0.0f", minlength2, maxlength2),
                        name = "Primer1", Primer3Path = Primer3Path,
                        thermo.param = thermo.param,
                        sequence_target = 110+ minlength2,
                        settings = settings)
      # If we have all Sequences needed we build Output binding with new Sequence info and knowed Exon info
      if(is.null(dim(p2))==FALSE){
        for (d in 1: dim(p2)[1]) {
          # We calculate primers positions in exons
          For1Seq <- p1[s,2]
          Rev1Seq <- p1[s,3]
          Rev2Seq <- p2[d,3]
          LastPosFor1 <- p1[s,6]+ p1[s,7] - 1
          LastPosFor2 <- NA
          FirstPosRev1 <- p1[s,8]-str_length(FExonSeq)-minlength-p1[s,9]+1
          FirstPosRev2 <- p2[d,8]-str_length(FExonSeq)-minlength2-p2[d,9]+1
          distinPrimers1 <- p1$PRIMER_PAIR_PRODUCT_SIZE[s]- minlength
          InfoC1 <- getDistanceseachPath(For1Exon,Rev1Exon,generaldata,distinPrimers1,SG)
          distinPrimers2 <- p2$PRIMER_PAIR_PRODUCT_SIZE[d]- minlength2
          InfoC2 <- getDistanceseachPath(For1Exon,Rev2Exon,generaldata,distinPrimers2,SG)
          Info <- c(InfoC1,InfoC2)
          # Distances from any path:
          DistancesP0 <- paste(Info[which(Info=="p0")-1][order(as.integer(Info[which(Info=="p0")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path1:
          DistancesP1 <- paste(Info[which(Info=="p1")-1][order(as.integer(Info[which(Info=="p1")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path2
          DistancesP2 <- paste(Info[which(Info=="p2")-1][order(as.integer(Info[which(Info=="p2")-1]),decreasing = FALSE)],collapse = " - ")
          Distances <- data.frame(DistancesP1,DistancesP2,DistancesP0)
          names(Distances) <-c("DistPath1","DistPath2","DistNoPath")
          FExons <- data.frame(FinalExons[n,1],FinalExons[n,2],FinalExons[n,4],FinalExons[n,3],FinalExons[n,5],FinalExons[n,7])
          colnames(FExons) <- colnames(FinalExons)[-6]
          PrSeq<-data.frame(For1Seq,NA,Rev2Seq,Rev1Seq,LastPosFor1,LastPosFor2,FirstPosRev2,FirstPosRev1)
          colnames(PrSeq)<-c("For1Seq","For2Seq","Rev1Seq","Rev2Seq","LastPosFor1","LastPosFor2","FirstPosRev1","FirstPosRev2")
          Fdata1<-rbind(Fdata1,cbind(PrSeq , FExons,Distances))
        }
      }
    }
  }
  return(unique(Fdata1))
}


#' @rdname InternalFunctions
PrimerSequenceCommonRev <-function(FinalExons,SG,
                                   generaldata,n,
                                   thermo.param,
                                   Primer3Path,
                                   settings)
  {
  Fdata1 <- data.frame()
  # Case1: Find sequence in the first Exon and once we have set that sequence we look for the sequence in the second exon.
  For1Exon <- FinalExons[n,1]
  For2Exon <- FinalExons[n,2]
  Rev1Exon <- FinalExons[n,3]
  # Get the ID for each exon where Primers are placed
  ExonF1 <- match(as.character(FinalExons[n,1]), SG$Edges$From)
  ExonF2 <- match(as.character(FinalExons[n,2]), SG$Edges$From)
  ExonR1 <- match(as.character(FinalExons[n,3]), SG$Edges$From)
  # Get de sequence of each exon
  Hsapienshg38  <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  FExonSeq1 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonF1], as.numeric(SG$Edges$Start[ExonF1]), as.numeric(SG$Edges$End[ExonF1])))
  FExonSeq2 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonF2], as.numeric(SG$Edges$Start[ExonF2]), as.numeric(SG$Edges$End[ExonF2])))
  SeqExonR  <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonR1], as.numeric(SG$Edges$Start[ExonR1]), as.numeric(SG$Edges$End[ExonR1])))
  # Call primer3 with two the sequence of two of the exons:
  minlength <- max(c(str_length(FExonSeq1),str_length(SeqExonR)))+1
  seq1 <- paste(as.character(FExonSeq1), paste(rep("N",minlength),collapse = "") ,as.character(SeqExonR),sep="")
  maxlength <- str_length(seq1)
  
  p1 <- callPrimer3(seq1, size_range = sprintf("%0.0f-%0.0f", minlength, maxlength),
                    name = "Primer1", Primer3Path = Primer3Path,
                    thermo.param = thermo.param,
                    sequence_target = 110,
                    settings = settings)
  # When a sequence is finded in first exon we look for de sequence in the other exon:
  if(is.null(dim(p1))==FALSE){
    for (s in 1: dim(p1)[1]) {
      pr=p1[s,3]
      minlength2 <- max(c(str_length(FExonSeq2),str_length(SeqExonR)))+1
      seq2 <- paste(as.character(FExonSeq2), paste(rep("N",minlength2),collapse = "") ,as.character(SeqExonR),sep="")
      maxlength2 <- str_length(seq2)
      p2 <- callPrimer3(seq2,threeprimers = TRUE,pr=pr,reverse = TRUE, size_range = sprintf("%0.0f-%0.0f", minlength2, maxlength2),
                        name = "Primer1", Primer3Path = Primer3Path,
                        thermo.param = thermo.param,
                        sequence_target = 110+minlength2,
                        settings = settings)
      # If we have all Sequences needed we build Output binding with new Sequence info and knowed Exon info
      if(is.null(dim(p2))==FALSE){
        for (d in 1: dim(p2)[1]) {
          # We calculate primers positions in exons
          For1Seq <- p1[s,2]
          For2Seq <- p2[d,2]
          Rev1Seq <- p1[s,3]
          # We calculate primers positions in exons
          LastPosFor1 <- p1[s,6]+ p1[s,7] - 1
          LastPosFor2 <- p2[d,6]+ p2[d,7] - 1
          FirstPosRev1 <- p1[s,8]-str_length(FExonSeq1)-minlength-p1[s,9]+1
          FirstPosRev2 <- NA
          distinPrimers1 <- p1$PRIMER_PAIR_PRODUCT_SIZE[s] - minlength
          InfoC1 <- getDistanceseachPath(For1Exon,Rev1Exon,generaldata,distinPrimers1,SG)
          distinPrimers2 <- p2$PRIMER_PAIR_PRODUCT_SIZE[d] - minlength2
          InfoC2 <- getDistanceseachPath(For2Exon,Rev1Exon,generaldata,distinPrimers2,SG)
          Info <- c(InfoC1,InfoC2)
          # Distances from any path:
          DistancesP0 <- paste(Info[which(Info=="p0")-1][order(as.integer(Info[which(Info=="p0")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path1:
          DistancesP1 <- paste(Info[which(Info=="p1")-1][order(as.integer(Info[which(Info=="p1")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path2
          DistancesP2 <- paste(Info[which(Info=="p2")-1][order(as.integer(Info[which(Info=="p2")-1]),decreasing = FALSE)],collapse = " - ")
          Distances <- data.frame(DistancesP1,DistancesP2,DistancesP0)
          names(Distances) <-c("DistPath1","DistPath2","DistNoPath")
          FExons <- data.frame(FinalExons[n,1],FinalExons[n,2],FinalExons[n,3],FinalExons[n,4],FinalExons[n,5],FinalExons[n,7])
          colnames(FExons) <- colnames(FinalExons)[-6]
          PrSeq<-data.frame(For1Seq,For2Seq,Rev1Seq,NA,LastPosFor1,LastPosFor2,FirstPosRev1,FirstPosRev2)
          colnames(PrSeq)<-c("For1Seq","For2Seq","Rev1Seq","Rev2Seq","LastPosFor1","LastPosFor2","FirstPosRev1","FirstPosRev2")
          Fdata1<-rbind(Fdata1,cbind(PrSeq , FExons,Distances))
        }
      }
    }
  }
  # 2ยบ Case: Find sequence in the second Exon and once we have set that sequence we look for the sequence in the first exon.
  # For that reason we swap exons:
  a <-  FinalExons[,2]
  FinalExons[,2] <- FinalExons[,1]
  FinalExons[,1] <- a
  
  For1Exon <- FinalExons[n,1]
  For2Exon <- FinalExons[n,2]
  Rev1Exon <- FinalExons[n,3]
  # Get the ID for each exon where Primers are placed
  ExonF1 <- match(as.character(For1Exon), SG$Edges$From)
  ExonF2 <- match(as.character(For2Exon), SG$Edges$From)
  ExonR1 <- match(as.character(Rev1Exon), SG$Edges$From)
  # Get de sequence of each exon
  FExonSeq1 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonF1], as.numeric(SG$Edges$Start[ExonF1]), as.numeric(SG$Edges$End[ExonF1])))
  FExonSeq2 <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonF2], as.numeric(SG$Edges$Start[ExonF2]), as.numeric(SG$Edges$End[ExonF2])))
  SeqExonR  <- as.character(getSeq(Hsapienshg38,SG$Edges$Chr[ExonR1], as.numeric(SG$Edges$Start[ExonR1]), as.numeric(SG$Edges$End[ExonR1])))
  # Call primer3 with two the sequence of two of the exons:
  minlength <- max(c(str_length(FExonSeq1),str_length(SeqExonR)))+1
  seq1 <- paste(as.character(FExonSeq1), paste(rep("N",minlength),collapse = "") ,as.character(SeqExonR),sep="")
  maxlength <- str_length(seq1)
  
  p1 <- callPrimer3(seq1, size_range = sprintf("%0.0f-%0.0f", minlength, maxlength),
                    name = "Primer1", Primer3Path = Primer3Path,
                    thermo.param = thermo.param,
                    sequence_target = 110,
                    settings = settings)
  # When a sequence is finded in first exon we look for de sequence in the other exon:
  if(is.null(dim(p1))==FALSE){
    for (s in 1: dim(p1)[1]) {
      pr=p1[s,3]
      minlength2 <- max(c(str_length(FExonSeq2),str_length(SeqExonR)))+1
      seq2 <- paste(as.character(FExonSeq2), paste(rep("N",minlength2),collapse = "") ,as.character(SeqExonR),sep="")
      maxlength2 <- str_length(seq2)
      p2 <- callPrimer3(seq2,threeprimers = TRUE,pr=pr,reverse = TRUE, size_range = sprintf("%0.0f-%0.0f", minlength2, maxlength2),
                        name = "Primer1", Primer3Path = Primer3Path,
                        thermo.param = thermo.param,
                        sequence_target = 110+minlength2,
                        settings = settings)
      # If we have all Sequences needed we build Output binding with new Sequence info and knowed Exon info
      if(is.null(dim(p2))==FALSE){
        for (d in 1: dim(p2)[1]) {
          For1Seq <- p1[s,2]
          For2Seq <- p2[d,2]
          Rev1Seq <- p1[s,3]
          # We calculate primers positions in exons
          LastPosFor1 <- p1[s,6]+ p1[s,7] - 1
          LastPosFor2 <- p2[d,6]+ p2[d,7] - 1
          FirstPosRev1 <- p1[s,8]-str_length(FExonSeq1)-minlength-p1[s,9]+1
          FirstPosRev2 <- NA
          distinPrimers1 <- p1$PRIMER_PAIR_PRODUCT_SIZE[s] - minlength
          InfoC1 <- getDistanceseachPath(For1Exon,Rev1Exon,generaldata,distinPrimers1,SG)
          distinPrimers2 <- p2$PRIMER_PAIR_PRODUCT_SIZE[d] -minlength2
          InfoC2 <- getDistanceseachPath(For2Exon,Rev1Exon,generaldata,distinPrimers2,SG)
          Info <- c(InfoC1,InfoC2)
          # Distances from any path:
          DistancesP0 <- paste(Info[which(Info=="p0")-1][order(as.integer(Info[which(Info=="p0")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path1:
          DistancesP1 <- paste(Info[which(Info=="p1")-1][order(as.integer(Info[which(Info=="p1")-1]),decreasing = FALSE)],collapse = " - ")
          # Distances from path2
          DistancesP2 <- paste(Info[which(Info=="p2")-1][order(as.integer(Info[which(Info=="p2")-1]),decreasing = FALSE)],collapse = " - ")
          Distances <- data.frame(DistancesP1,DistancesP2,DistancesP0)
          names(Distances) <-c("DistPath1","DistPath2","DistNoPath")
          FExons <- data.frame(FinalExons[n,2],FinalExons[n,1],FinalExons[n,3],FinalExons[n,4],FinalExons[n,5],FinalExons[n,7])
          colnames(FExons) <- colnames(FinalExons)[-6]
          PrSeq<-data.frame(For2Seq,For1Seq,Rev1Seq,NA,LastPosFor2,LastPosFor1,FirstPosRev1,FirstPosRev2)
          colnames(PrSeq)<-c("For1Seq","For2Seq","Rev1Seq","Rev2Seq","LastPosFor1","LastPosFor2","FirstPosRev1","FirstPosRev2")
          Fdata1<-rbind(Fdata1,cbind(PrSeq , FExons,Distances))
        }
      }
    }
  }
  return(unique(Fdata1))
}





































