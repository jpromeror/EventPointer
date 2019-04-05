#' CDF file creation for EventPointer (MultiPath)
#'
#' Generates the CDF file to be used under the aroma.affymetrix framework.
#'
#' @param input Reference transcriptome used to build the CDF file. Must be one of Ensembl,
#'              UCSC or GTF.
#' @param inputFile If input is GTF, inputFile should point to the GTF file to be used.
#' @param PSR Path to the Exon probes txt file
#' @param Junc Path to the Junction probes txt file
#' @param PathCDF Directory where the output will be saved
#' @param microarray Microarray used to create the CDF file. Must be one of: HTA-2_0,
#'                  ClariomD, RTA or MTA
#' @param paths Maximum number of paths of the events to find.
#'
#' @return The function displays a progress bar to show the user the progress of the function.
#' However, there is no value returned in R as the function creates three files that are used
#' later by other EventPointer functions. 1) EventsFound.txt : Tab separated file with all
#' the information of all the alternative splcing events found. 2) .flat file : Used to build
#' the corresponding CDF file. 3) .CDF file: Output required for the aroma.affymetrix preprocessing
#' pipeline. Both the .flat and .CDF file take large ammounts of memory in the hard drive, it is
#' recommended to have at least 1.5 GB of free space.
#'
#' @examples
#'    PathFiles<-system.file('extdata',package='EventPointer')
#'    DONSON_GTF<-paste(PathFiles,'/DONSON.gtf',sep='')
#'    PSRProbes<-paste(PathFiles,'/PSR_Probes.txt',sep='')
#'    JunctionProbes<-paste(PathFiles,'/Junction_Probes.txt',sep='')
#'    Directory<-tempdir()
#'    microarray<-'HTA-2_0'
#'
#'   # Run the function
#'
#'    CDFfromGTF_Multipath(input='AffyGTF',inputFile=DONSON_GTF,PSR=PSRProbes,Junc=JunctionProbes,
#'                         PathCDF=Directory,microarray=microarray,paths=3)
#'
#' @export
#' @import Matrix
#' @import SGSeq
#' @import SummarizedExperiment
#' @importFrom affxparser writeCdf
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom GenomicFeatures makeTxDbFromBiomart makeTxDbFromUCSC makeTxDbFromGFF
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar combn write.table read.table
#' @importFrom stringr str_count
#' @importFrom GenomeInfoDb 'seqlevelsStyle<-' seqlevelsStyle seqnames
#' @importFrom igraph graph_from_data_frame as_adj clusters graph_from_adjacency_matrix
#' @importFrom igraph graph.data.frame
#' @importFrom MASS Null ginv
#' @importFrom stats dist qnorm quantile runif rnorm
#' @importFrom nnls nnls
#' @importFrom limma lmFit contrasts.fit eBayes topTable voom
#' @importFrom matrixStats rowMaxs colCumsums
#' @importFrom RBGL connectedComp
#' @importFrom methods as
#' @importFrom graph ftM2graphNEL
#' @importFrom prodlim row.match
#' @importFrom 'methods' is
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors queryHits subjectHits


CDFfromGTF_Multipath <- function(input = "Ensembl", 
    inputFile = NULL, PSR, Junc, PathCDF, 
    microarray = NULL, paths = 2) {
    
    # CDFfromGTF is the corresponding
    # function to create a CDF file based on
    # the GTF file provided that contains the
    # 'reference' trasncriptome to be used to
    # identify alternative splicing events.
    
    # Input Files: input -> GTF File of the
    # transcriptome Dsource -> Information
    # about the CDF (not required) PSR -> Txt
    # file with the information of the PSRs
    # of the array Junc -> Txt file with the
    # information of junction probes in the
    # array
    
    
    # Crate TxDB object that contains all the
    # information contained in the gtf file
    # that was given in the input, the
    # function makeTxDbFromGFF corresponds to
    # the GenomicFeatures R package
    
    cat("Creating SG Information...")
    
    # Possible Arrays: HTA-2_0, ClariomD, RTA
    # and MTA
    
    if (is.null(microarray)) {
        stop("Microarray field empty")
    }
    
    if (microarray == "HTA-2_0") {
        if (input == "Ensembl") {
            TxDb <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl", 
                host = "grch37.ensembl.org")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "UCSC") {
            TxDb <- makeTxDbFromUCSC(genome = "hg19", 
                tablename = "knownGene")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "AffyGTF" | input == 
            "CustomGTF") {
            if (is.null(inputFile)) {
                stop("inputFile parameter is empty")
            }
            
            TxDb <- makeTxDbFromGFF(file = inputFile, 
                format = "gtf", dataSource = "Custom GTF")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else {
            stop("Unknown reference genome")
        }
    } else if (microarray == "ClariomD") {
        if (input == "Ensembl") {
            TxDb <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "UCSC") {
            TxDb <- makeTxDbFromUCSC(genome = "hg38", 
                tablename = "knownGene")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "GTF") {
            if (is.null(inputFile)) {
                stop("inputFile parameter is empty")
            }
            
            TxDb <- makeTxDbFromGFF(file = inputFile, 
                format = "gtf", dataSource = "Custom GTF")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else {
            stop("Unknown reference genome")
        }
    } else if (microarray == "RTA") {
        if (input == "Ensembl") {
            TxDb <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "rnorvegicus_gene_ensembl")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "UCSC") {
            TxDb <- makeTxDbFromUCSC(genome = "rn6", 
                tablename = "refGene")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "GTF") {
            if (is.null(inputFile)) {
                stop("inputFile parameter is empty")
            }
            
            TxDb <- makeTxDbFromGFF(file = inputFile, 
                format = "gtf", dataSource = "Custom GTF")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else {
            stop("Unknown reference genome")
        }
    } else if (microarray == "MTA") {
        if (input == "Ensembl") {
            TxDb <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "mmusculus_gene_ensembl")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "UCSC") {
            TxDb <- makeTxDbFromUCSC(genome = "mm10", 
                tablename = "knownGene")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else if (input == "GTF") {
            if (is.null(inputFile)) {
                stop("inputFile parameter is empty")
            }
            
            TxDb <- makeTxDbFromGFF(file = inputFile, 
                format = "gtf", dataSource = "Custom GTF")
            TranscriptFeatures <- convertToTxFeatures(TxDb)
            SplicingGraphFeatures <- convertToSGFeatures(TranscriptFeatures)
        } else {
            stop("Unknown reference genome")
        }
    } else {
        stop("Microarray should be: HTA-2_0, ClariomD, MTA or RTA")
    }
    
    # Read Information for the probes in the
    # array
    
    cat("\nReading Information On Probes...")
    
    stopifnot(!is.null(PSR) & !is.null(Junc))
    # Read ProbeSets TXT
    ProbeSets <- read.delim(file = PSR, sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE)
    
    # Arrange the information in the txt file
    # to have the information in the way the
    # algorithm needs it
    
    ProbeSets <- PrepareProbes(ProbeSets, 
        "PSR")
    
    # Read Junctions TXT
    Junctions <- read.delim(file = Junc, 
        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
    # Arrange the information in the txt file
    # to have the information in the way the
    # algorithm needs it
    
    Junctions <- PrepareProbes(Junctions, 
        "Junction")
    
    
    if (input == "Ensembl" | input == "UCSC") {
        Junc <- makeGRangesFromDataFrame(Junctions)
        PSRs <- makeGRangesFromDataFrame(ProbeSets)
        
        seqlevelsStyle(Junc) <- seqlevelsStyle(SplicingGraphFeatures)
        seqlevelsStyle(PSRs) <- seqlevelsStyle(SplicingGraphFeatures)
        
        
        Overlap_Junc <- findOverlaps(Junc, 
            SplicingGraphFeatures)
        
        GeN <- geneName(SplicingGraphFeatures[subjectHits(Overlap_Junc)])
        GeN <- sapply(GeN, function(X) {
            A <- X[1]
            return(A)
        })
        
        GeN_iix <- cbind(queryHits(Overlap_Junc), 
            GeN)
        GeN_iix <- unique(GeN_iix)
        Junctions[as.numeric(GeN_iix[, 1]), 
            "Gene"] <- GeN_iix[, 2]
        
        
        Overlap_PSR <- findOverlaps(PSRs, 
            SplicingGraphFeatures)
        
        GeN <- geneName(SplicingGraphFeatures[subjectHits(Overlap_PSR)])
        GeN <- sapply(GeN, function(X) {
            A <- X[1]
            return(A)
        })
        
        GeN_iix <- cbind(queryHits(Overlap_PSR), 
            GeN)
        GeN_iix <- unique(GeN_iix)
        ProbeSets[as.numeric(GeN_iix[, 1]), 
            "Gene"] <- GeN_iix[, 2]
    }
    
    cat("Done")
    
    # Get genes with probes and index them to
    # get locations in Junctions and Probes
    
    cat("\nIndexing Genes and Probes...")
    
    # Obtain the genenames('TCXXXX.hg') and
    # geneIDs (1,2...) for all the genes in
    # the transcriptome used
    
    geneIds <- geneID(SplicingGraphFeatures)
    Genes <- as.list(geneName(SplicingGraphFeatures))
    
    # As certain genes share geneID we repeat
    # each geneID as many times it is used by
    # a gene
    LLs <- unlist(lapply(Genes, length))
    geneIds <- rep(geneIds, LLs)
    Genes <- unlist(Genes)
    
    # Create a matrix with both Gene Names
    # and Gene IDs
    GeneIndex <- cbind(Genes, geneIds)
    
    # Get Genes with Junctions probes, PSR
    # probes and already located in the
    # GeneIndex matrix
    KeptGenes <- intersect(GeneIndex[, 1], 
        intersect(Junctions[, "Gene"], ProbeSets[, 
            "Gene"]))
    
    # Order GeneIndex as the genes in
    # KeptGenes
    iix <- match(KeptGenes, GeneIndex[, 1])
    GeneIndex <- GeneIndex[iix, , drop = FALSE]
    
    # For the probe variables, we only keep
    # the genes that appear in KeptGenes
    ii.P <- match(ProbeSets[, "Gene"], KeptGenes)
    ii.P <- which(!is.na(ii.P))
    ProbeSets <- ProbeSets[ii.P, , drop = FALSE]
    
    ii.J <- match(Junctions[, "Gene"], KeptGenes)
    ii.J <- which(!is.na(ii.J))
    Junctions <- Junctions[ii.J, , drop = FALSE]
    
    # Order all the variables in the same
    # way, so the genes in different
    # variables have the same 'index'
    Ord <- order(KeptGenes)
    GeneIndex <- GeneIndex[Ord, , drop = FALSE]
    KeptGenes <- KeptGenes[Ord]
    ProbeSets <- ProbeSets[order(ProbeSets[, 
        "Gene"]), , drop = FALSE]
    Junctions <- Junctions[order(Junctions[, 
        "Gene"]), , drop = FALSE]
    
    # With rle we get how many times a
    # certain value appears in the variable
    # this way we can get the indices that
    # correspond to the probes for a
    # particular gene.
    
    # For example, with the indices we can
    # now that a gene X has PSR probes
    # between lines 235 and 312 of the
    # variable.
    
    Junc_rle <- rle(Junctions[, "Gene"])
    PS_rle <- rle(ProbeSets[, "Gene"])
    
    indexE_Junc <- cumsum(Junc_rle$lengths)
    indexS_Junc <- c(1, indexE_Junc[seq_len((length(indexE_Junc) - 
        1))] + 1)
    
    indexE_PS <- cumsum(PS_rle$lengths)
    indexS_PS <- c(1, indexE_PS[seq_len((length(indexE_PS) - 
        1))] + 1)
    
    cat("Done")
    
    # Main loop to apply functions to all the
    # genes in the GTF
    
    # The list Result will contain all the
    # information of the Events found
    
    Result <- vector("list", length = length(KeptGenes))
    Flat <- vector("list", length = length(KeptGenes))
    
    # Set a progress bar to track the overall
    # progress
    gc()
    
    cat("\nFinding Events...")
    cat("\n")
    
    pb <- txtProgressBar(min = 0, max = nrow(GeneIndex), 
        style = 3)
    
    for (jj in seq_len(nrow(GeneIndex))) {
        # cat(jj,' of 30476 \n')
        setTxtProgressBar(pb, jj)
        
        Gene <- GeneIndex[jj, 1]
        
        # Obtain the required information to
        # build the SG, the Exons correspond to
        # Nodes while Junctions correspond to
        # edges. It is important to take into
        # account that the Exons are already
        # 'divided' into subexons, non
        # overlapping regions between transcripts
        
        
        SG_Gene <- SplicingGraphFeatures[unlist(geneID(SplicingGraphFeatures)) == 
            GeneIndex[jj, 2]]
        
        
        # Get the number of transcripts that
        # appear for the gene that is being
        # analyzed
        
        Txx <- unique(unlist(as.data.frame(SG_Gene)[, 
            "txName"]))
        
        # In order to have alternative splicing a
        # gene must have at least more than one
        # exon and more than one transcript
        
        if (nrow(as.data.frame(SG_Gene)) != 
            1 & length(Txx) != 1) {
            
            # Function to obtain the information of
            # the SG, the information returned
            # Corresponds to the Edges, Adjacency
            # matrix and Incidence Matrix of the SG
            # that is being analyzed
            
            # SG <- SG_Info(SG_Gene)
            
            # SG <- SG_creation(SG_Gene)
            
            SG <- SG_creation_RNASeq(SG_Gene)
            
            # Using Ax=b with A the incidence matrix
            # and b=0, we solve to obtain the null
            # space and obtain the flux through each
            # edge if we relate the SG with an
            # hydraulic installation
            
            # The new parameter (ncol) is used to
            # generate 10 fluxes to ensure that the
            # eventdetection is exaclty the same
            # every step.
            
            randSol <- getRandomFlow(SG$Incidence, 
                ncol = 10)
            
            # To find events, we need to have fluxes
            # different from 1 in different edges
            
            if (any(round(randSol) != 1)) {
                # Identify the alternative splicing
                # events using the fluxes obtained
                # previously.  The number of fluxes
                # correspond to the number of edges so we
                # can relate events with edges
                Events <- findTriplets2(SG$Incidence, 
                  paths = paths, randSol)
                
                
                if (nrow(Events$multipaths) > 
                  0) {
                  twopaths <- which(rowSums(Events$multipaths != 
                    0) == 4)
                  Events <- getEventMultiPaths(Events, 
                    SG, twopaths, paths)
                  
                  # Events<-TrimMultiEvents(Events,paths)
                  if (length(Events) > 0) {
                    Events <- ClassifyEvents(SG, 
                      Events, twopaths)
                    
                    # Get PSRs that map to the gene that is
                    # being used
                    
                    PSR_Gene <- ProbeSets[indexS_PS[jj]:indexE_PS[jj], 
                      ]
                    
                    # Get junctions that map to the gene that
                    # is being used
                    
                    Junc_Gene <- Junctions[indexS_Junc[jj]:indexE_Junc[jj], 
                      ]
                    
                    
                    # Obtain the corresponding
                    # PSRs/JunctionProbes for each of the
                    # paths of every event
                    
                    Events <- annotateEventsMultipath(Events, 
                      PSR_Gene, Junc_Gene, 
                      GeneIndex[jj, 1], paths)
                    
                    Result[[jj]] <- Events[[1]]
                    Flat[[jj]] <- Events[[2]]
                  }
                }
            }
        }
    }
    
    
    close(pb)
    
    # Create .flat file used in flat2CDF
    # function
    
    cat("\nCreating .flat ...")
    Result <- Filter(Negate(is.null), Result)
    Result <- do.call(rbind, Result)
    colnames <- "colnames(Result) <- c('Affy Gene Id', 'Gene Name', 'Event Number', 'Event Type',
                                      'Genomic Position','Num of Paths',"
    for (kk in seq_len(paths)) {
        colnames <- paste0(colnames, "'Path ", 
            kk, "', ")
    }
    colnames <- paste0(colnames, "'Path Reference',")
    for (kk in seq_len(paths)) {
        colnames <- paste0(colnames, "'Probes P", 
            kk, "', ")
    }
    colnames <- paste0(colnames, "'Probes Ref')")
    
    eval(parse(text = colnames))
    
    # if(input=='GTF') {
    # GeneNames<-read.delim(file=paste(GeneNamesFile,sep=''),sep='\t',header=T)
    # iix<-match(Result[,1],GeneNames[,1])
    # Nmb<-as.vector(GeneNames[iix,2])
    # Result[,2]<-Nmb }else
    # if(input=='Ensembl' | input=='UCSC') {
    # Result[,2]<-Result[,1] }
    
    
    Flat <- Filter(Negate(is.null), Flat)
    Flat <- do.call(rbind, Flat)
    colnames(Flat) <- c("Probe_ID", "X", 
        "Y", "Probe_Sequence", "Group_ID", 
        "Unit_ID")
    
    
    # tag<-'r'
    
    ROWS <- 
72
    COLS <- 2680
    
    cat("Done")
    
    write.table(Result, file = paste(PathCDF, 
        "/EventsFound_", microarray, ".txt", 
        sep = ""), sep = "\t", quote = FALSE, 
        col.names = TRUE, row.names = FALSE)
    
    write.table(Flat, file = paste(PathCDF, 
        "/", microarray, ".flat", sep = ""), 
        sep = "\t", quote = FALSE, col.names = TRUE, 
        row.names = FALSE)
    
    # Cant change wd during BioCCheck
    # setwd(PathCDF)
    
    flat2Cdf(file = paste(PathCDF, "/", microarray, 
        ".flat", sep = ""), chipType = microarray, 
        rows = ROWS, cols = COLS, Directory = PathCDF)
}
