#' @noRd
NULL

# ep_st_wf.R ################
#' @rdname InternalFunctions
voomEventPointerST <- function(PSI,Design,Contrast){
  # Compute the abundance as the minimum value of the three paths
  # Using the minimum
  PSI_boots <- PSI$PSI
  Events <- PSI$ExpEvs
  abundance <- unlist(lapply(Events,
                             function(y) min(colMeans2(y))))
  
  averageP1 <- unlist(lapply(Events,
                             function(y) colMeans2(y)[1]))
  averageP2 <- unlist(lapply(Events,
                             function(y) colMeans2(y)[2]))
  averageRef <- unlist(lapply(Events,
                              function(y) colMeans2(y)[3]))
  
  # Remove some values
  dummy <- (rowSds(PSI_boots[,1,],useNames =T)<1e-6)
  dummy[is.na(dummy)] <- TRUE
  Quitar <- dummy
  
  PSI <- PSI_boots[!Quitar,1,]
  abundance <- abundance[!Quitar]
  averageP1 <- averageP1[!Quitar]
  averageP2 <- averageP2[!Quitar]
  averageRef <- averageRef[!Quitar]
  betaest <- solve(base::crossprod(Design),t(Design)%*%t(PSI))
  error <- t(PSI) - Design %*% betaest
  errorSE <- sqrt(colSums(error^2)/-diff(dim(Design)))
  
  # If the value of PSI is close to 1 or close to zero, the variance is small
  p <- rowMeans(PSI)
  varEst <- p*(1-p)
  # Let's combine both sources of information
  modelo <- lm((errorSE)^.5 ~ I(log2(abundance)) +
                 I(log2(averageRef)) + I(log2(averageP1)) + I(log2(averageP2)))
  
  overallVar <- predict(modelo)
  
  overallVar <- 1 + max(overallVar) - overallVar
  
  logabd <- overallVar
  PSI_changed <- (5+PSI) * logabd/rowMeans(5+PSI)
  
  # We will try to apply voom to these data.
  # Remove NAs for the time being
  
  modelo <- voom2(PSI_changed, Design, plot=F, span=.1)
  fit <- suppressWarnings(lmFit(modelo, Design, method="robust"))
  fit <- contrasts.fit(fit, Contrast)
  fit <- suppressWarnings(eBayes(fit, robust=TRUE,proportion = .1))
  deltaPSI <- t(Contrast) %*% (solve(base::crossprod(Design), t(Design))) %*% t(PSI)
  colnames(deltaPSI) <- rownames(PSI_changed)
  res <- list()
  for(coef in c(1:ncol(fit$coefficients))){
    result <- topTable(fit, coef = coef, number = Inf)
    result <- data.frame(deltaPSI = deltaPSI[coef,rownames(result)], result)
    res[[coef]] <- result
  }
  return(res)
}


#EventPointerStats_BAM.R ################

#' @rdname InternalFunctions
voomEventPointerBAM <- function(PSI_boots,Events,Design,Contrast){
  # Compute the abundance as the minimum value of the three paths
  # Using the minimum
  abundance <- unlist(sapply(Events,
                             function(y) sapply(y, function(x) min(rowMeans2(x$Counts)))))
  averageRef <- unlist(sapply(Events,
                              function(y) sapply(y, function(x) rowMeans2(x$Counts)[3])))
  averageP1 <- unlist(sapply(Events,
                             function(y) sapply(y, function(x) rowMeans2(x$Counts)[1])))
  averageP2 <- unlist(sapply(Events,
                             function(y) sapply(y, function(x) rowMeans2(x$Counts)[2])))
  
  # Remove some values
  dummy <- (rowSds(PSI_boots[,1,],useNames =T )<1e-6)
  dummy[is.na(dummy)] <- TRUE
  Quitar <- dummy
  
  PSI <- PSI_boots[!Quitar,1,]
  abundance <- abundance[!Quitar]
  averageP1 <- averageP1[!Quitar]
  averageP2 <- averageP2[!Quitar]
  averageRef <- averageRef[!Quitar]
  betaest <- solve(base::crossprod(Design),t(Design)%*%t(PSI))
  error <- t(PSI) - Design %*% betaest
  errorSE <- sqrt(colSums(error^2)/-diff(dim(Design)))
  
  # If the value of PSI is close to 1 or close to zero, the variance is small
  p <- rowMeans(PSI)
  varEst <- p*(1-p)
  # Let's combine both sources of information
  modelo <- lm((errorSE)^.5 ~ I(log2(abundance)) +
                 I(log2(averageRef)) + I(log2(averageP1)) + I(log2(averageP2)))
  
  overallVar <- predict(modelo)
  
  overallVar <- 1 + max(overallVar) - overallVar
  
  logabd <- overallVar
  PSI_changed <- (5+PSI) * logabd/rowMeans(5+PSI)
  
  # We will try to apply voom to these data.
  # Remove NAs for the time being
  
  modelo <- voom2(PSI_changed, Design, plot=F, span=.1)
  fit <- suppressWarnings(lmFit(modelo, Design, method="robust"))
  fit <- contrasts.fit(fit, Contrast)
  fit <- suppressWarnings(eBayes(fit, robust=TRUE,proportion = .1))
  deltaPSI <- t(Contrast) %*% (solve(base::crossprod(Design), t(Design))) %*% t(PSI)
  colnames(deltaPSI) <- rownames(PSI_changed)
  res <- list()
  for(coef in c(1:ncol(fit$coefficients))){
    result <- topTable(fit, coef = coef, number = Inf)
    result <- data.frame(deltaPSI = deltaPSI[coef,rownames(result)], result)
    res[[coef]] <- result
  }
  return(res)
}

#' @rdname InternalFunctions
voom2 <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none",
                   block = NULL, correlation = NULL, weights = NULL, span = 0.5,
                   plot = FALSE, save.plot = FALSE, keepMax=T){
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >0)
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size))
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts)))
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts)))
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  m <- min(counts)
  if (is.na(m))
    stop("NA counts not allowed")
  if (m < 0)
    stop("Negative counts not allowed")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size))
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- limma:::normalizeBetweenArrays(y, method = normalize.method)
  y <- counts
  fit <- lmFit(y, design, block = block, correlation = correlation,
               weights = weights)
  if (is.null(fit$Amean))
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  NWithReps <- sum(fit$df.residual > 0L)
  if (NWithReps < 2L) {
    if (NWithReps == 0L)
      warning("The experimental design has no replication. Setting weights to 1.")
    if (NWithReps == 1L)
      warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    out$design <- design
    if (is.null(out$targets))
      out$targets <- data.frame(lib.size = lib.size)
    else out$targets$lib.size <- lib.size
    return(new("EList", out))
  }
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if(keepMax) {
    I <- which.max(l$y)
    ymax <- l$y[I]
    l$y[1:I] <- ymax
  }
  if (plot) {
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
         pch = 16, cex = 0.25)
    title("voom: Mean-variance trend")
    lines(l, col = "red")
  }
  f <- approxfun(l, rule = 2, ties = list("ordered", mean))
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coefficients[, j, drop = FALSE] %*%
      t(fit$design[, j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coefficients %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets))
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )",
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  new("EList", out)
}


#getSGFeatureCounts_SGSeqMod.R #################

#' @rdname InternalFunctions
sampleWhichCount <- function(sample_info, list_features) {
  list_which <- unlist(lapply(list_features, function(chrFeature) {
    if (length(chrFeature) != 0) {
      which <- range(chrFeature)
      seq <- as.character(unique(seqnames(which)))
      newWhich <- GRanges(
        seqnames = Rle(c(seq), c(1)),
        ranges = IRanges(min(start(which)), end = max(end(which))),
        strand = Rle(strand(c("-")), c(1))
      )
      newWhich
    }
  }))
  
  res <- c()
  for (sample in c(1:length(sample_info$file_bam))) {
    res <- c(
      res,
      lapply(
        X = list_which,
        FUN = listCHRFileCount,
        valSample = sample_info[sample, ],
        list_features = list_features
      )
    )
    
  }
  return(res)
}

#' @rdname InternalFunctions
listCHRFileCount <- function(which, valSample, list_features) {
  seq <- as.character(seqnames(which))
  features <- list_features[[seq]]
  file_bam <- valSample$file_bam
  paired_end <- valSample$paired_end
  sample_name <- valSample$sample_name
  addWhich <-list(file_bam, paired_end, sample_name, which, features)
  names(addWhich) <-c("file_bam", "paired_end", "sample_name", "which", "features")
  return(addWhich)
}

#' @rdname InternalFunctions
getSGFeatureCounts <- function(sample_info,
                               features,
                               min_anchor = 6,
                               counts_only = FALSE,
                               retain_coverage = FALSE,
                               verbose = FALSE,
                               cores = 1)
{
  
  if (!is(features, "SGFeatures")) {
    stop("features must be an SGFeatures object")
  }
  #Get which from features
  
  mergeFeatures <- features
  list_range <- range(mergeFeatures)
  hits <- findOverlaps(mergeFeatures, list_range)
  list_index <- split(queryHits(hits), subjectHits(hits))
  mergeFeatures <- mergeFeatures[queryHits(hits)]
  list_features <- split(mergeFeatures, seqnames(mergeFeatures))
  
  list_element_count <- sampleWhichCount(sample_info, list_features)
  clC <- makePSOCKcluster(cores, outfile = "reportCount.txt")
  list_counts <- clusterApplyLB(
    cl = clC,
    x = list_element_count,
    fun = getSGFeatureCountsTotal,
    min_anchor = min_anchor,
    retain_coverage = retain_coverage,
    verbose = verbose
  )
  list_counts <- unlist(list_counts, recursive = FALSE)
  
  on.exit(stopCluster(cl = clC))
  CountPerSample <- list()
  for (sample in sample_info$sample_name) {
    listCountSample <- list()
    for (chr in c(1:length(list_counts))) {
      nameCHR <- names(list_counts[chr])
      chr <- list_counts[[chr]]
      listStrand <- list()
      for (chrStrand in chr) {
        if (sample == chrStrand$sample_name) {
          listStrand <- c(listStrand, list(chrStrand$counts))
        }
      }
      if (length(listStrand) != 0) {
        listCountSample[[nameCHR]] <- listStrand
      }
    }
    
    listCountSample <-
      listCountSample[unique(as.character(seqnames(mergeFeatures)))]
    counts <- c()
    if (retain_coverage) {
      counts <- do.call(rbind, listCountSample)
      counts <- counts[order(unlist(list_index)),]
      
    } else {
      counts <- unlist(listCountSample, use.names = FALSE)
      counts <- counts[order(unlist(list_index))]
      
    }
    if (length(CountPerSample) == 0) {
      CountPerSample <- list(counts)
    } else{
      CountPerSample <- c(CountPerSample, list(counts))
    }
  }
  
  counts <- do.call(cbind, CountPerSample)
  
  if (counts_only)
    return(counts)
  
  sgfc <- makeSGFeatureCounts(
    rowRanges = features,
    colData = sample_info,
    counts = counts,
    min_anchor = min_anchor
  )
  
  return(sgfc)
  
}

#' @rdname InternalFunctions
getSGFeatureCountsTotal <- function(mainCount,min_anchor, retain_coverage, verbose)
{
  features <- mainCount$features
  sample_name <- mainCount$sample_name
  file_bam <- mainCount$file_bam
  paired_end <- mainCount$paired_end
  sample_name <- mainCount$sample_name
  which <- mainCount$which
  whichMod<- which
  chrName <- as.character(which@seqnames)
  seqlevel <- as.character(seqnames(which))
  strand <- unique(as.character(strand(features)))
  print(paste("READ CHR",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
  pairGap <- readGapPair(file_bam, paired_end, which, sample_name, verbose)
  listCount <- list()
  if (length(strand)==2) {
    featuresPlus <- features[strand(features) == "+"]
    featuresMinus <- features[strand(features) == "-"]
    print(paste("Start + strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    listCount[[chrName]] <- list("+" = processCounts(pairGap$gapPlus, featuresPlus, strand = "+",sample_name, min_anchor, retain_coverage, verbose))
    print(paste("End + strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    
    gc()
    print(paste("Start - strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    listCount[[chrName]] <- c(listCount[[chrName]],list("-" = processCounts(pairGap$gapMinus, featuresMinus, strand = "-",sample_name, min_anchor, retain_coverage, verbose)))
    print(paste("End - strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    
    gc()
  }else{
    if (strand == "+"){
      featuresPlus <- features[strand(features) == "+"]
      listCount[[chrName]] <- list("+" = processCounts(pairGap$gapPlus, featuresPlus, strand = "+",sample_name, min_anchor, retain_coverage, verbose))
      gc()
    }else{
      featuresMinus <- features[strand(features) == "-"]
      listCount[[chrName]] <- list("-" = processCounts(pairGap$gapMinus, featuresMinus, strand = "-",sample_name, min_anchor, retain_coverage, verbose))
      gc()
    }
  }
  
  return(listCount)
  
}

#' @rdname InternalFunctions
processCounts <- function(gap, features, strand,
                          sample_name, min_anchor, retain_coverage, verbose){
  
  which <- range(features)
  chrName <- as.character(which@seqnames)
  if (isEmpty(gap$frag_exonic)) {
    res_count <- list(counts = as.integer(rep.int(0,length(features))), sample_name=sample_name)
    rm(gap)
    gc()
  }else{
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron
    
    rm(gap)
    gc()
    
    rangeExon <- unlist(range(frag_exonic))
    
    selectionRangeExon <- range(which(end(rangeExon) >= min(start(features)) & start(rangeExon) <= max(end(features))))
    print(c(selectionRangeExon, chrName, strand))
    if (is.infinite(selectionRangeExon[1]) | is.infinite(selectionRangeExon[2]) ) {
      rm(frag_exonic)
      rm(frag_intron)
      gc()
      res_count <- list(counts = as.integer(rep.int(0,length(features))), sample_name=sample_name)
    }else{
      frag_exonic <- frag_exonic[c(selectionRangeExon[1]:selectionRangeExon[2])]
      frag_intron <- frag_intron[c(selectionRangeExon[1]:selectionRangeExon[2])]
      
      ir <- SGSeq:::extractRangesFromFeatures(features)
      
      ## extract feature type and spliced boundaries
      type <- mcols(ir)$type
      spliceL <- mcols(ir)$spliceL
      spliceR <- mcols(ir)$spliceR
      
      i_J <- which(type == "J")
      i_E <- which(type == "E")
      i_S <- which(type %in% c("spliceL", "spliceR"))
      
      N <- rep(NA_integer_, length(ir))
      if (length(i_J) > 0) {
        N[i_J] <- junctionCompatible(ir[i_J], frag_exonic, frag_intron,min_anchor)
      }
      if (length(i_E) > 0) {
        E_index <- exonCompatible(ir[i_E], spliceL[i_E], spliceR[i_E],frag_exonic, frag_intron, FALSE)
        N[i_E] <- elementNROWS(E_index)
      }
      if (length(i_S) > 0) {
        N[i_S] <- splicesiteOverlap(ir[i_S],sub("splice", "", type[i_S], fixed = TRUE),frag_exonic, frag_intron, min_anchor, "unspliced")
      }
      if (retain_coverage) {
        
        counts <- DataFrame(N = N)
        counts$N_splicesite <- IntegerList(vector("list", nrow(counts)))
        counts$coverage <- RleList(IntegerList(vector("list", nrow(counts))))
        if (length(i_J) > 0) {
          counts$N_splicesite[i_J] <- splicesiteCounts(ir[i_J],frag_exonic, frag_intron, min_anchor, "junction", "all")
        }
        if (length(i_E) > 0) {
          counts$N_splicesite[i_E] <- splicesiteCounts(ir[i_E],frag_exonic, frag_intron, min_anchor, "exon", "spliced")
          counts$coverage[i_E] <- exonCoverage(ir[i_E], E_index,frag_exonic)
        }
        if (strand == "-") {
          counts$N_splicesite <- endoapply(counts$N_splicesite, rev)
          counts$coverage <- endoapply(counts$coverage, rev)
        }
      } else {
        counts <- N
      }
      if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))
      rm(N)
      rm(E_index)
      res_count <- list(counts = counts, sample_name=sample_name)
    }
  }
  return(res_count)
  
}



# OverlapModule.R #################



#' @rdname InternalFunctions
junctionCompatible <- function(junctions, frag_exonic, frag_intron,
                               min_anchor, counts = TRUE)
{
  
  frag_intron <- filterIntrons(frag_intron, frag_exonic, min_anchor)
  
  introns <- junctions - 1
  
  hits <- findOverlapsRanges(introns, frag_intron, "equal")
  junction_index <- as.list(hits)
  rm(frag_exonic)
  rm(frag_intron)
  rm(junctions)
  gc(full=TRUE)
  if (counts) elementNROWS(junction_index)
  else junction_index
  
}

#' @rdname InternalFunctions
filterIntrons <- function(frag_intron, frag_exonic, min_anchor)
{
  unlisted_intron <- unlist(frag_intron)
  
  f_L <- as(flank(unlisted_intron, min_anchor, TRUE), "IRangesList")
  f_L <- intersect(f_L, frag_exonic[togroup0(frag_intron)])
  w_L <- sum(width(f_L))
  f_R <- as(flank(unlisted_intron, min_anchor, FALSE), "IRangesList")
  f_R <- intersect(f_R, frag_exonic[togroup0(frag_intron)])
  w_R <- sum(width(f_R))
  i <- which(w_L == min_anchor & w_R == min_anchor)
  
  filtered <- setNames(split(unlisted_intron[i],factor(togroup0(frag_intron)[i], seq_along(frag_intron))), NULL)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(filtered)
  
}

##' Identify fragments compatible with exons.


#' @rdname InternalFunctions
dimMatrixFix <- function(matrix1, matrix2){
  if(dim(matrix1)[1] != dim(matrix2)[1]){
    if(nrow(matrix1) > nrow(matrix2)){
      rows_to_add <- nrow(matrix1) - nrow(matrix2)
      matrix2 <- rbind(matrix2,
                       sparseMatrix(i = NULL,j=NULL,dims = c(rows_to_add,ncol(matrix2)))
      )
    }
  }
  return(matrix2)
}



#' @rdname InternalFunctions
splicesiteOverlap <- function(splicesites, side, frag_exonic, frag_intron,min_anchor, include = c("all", "spliced", "unspliced"), counts = TRUE){
  
  include <- match.arg(include)
  
  if (length(side) == 1) side <- rep(side, length(splicesites))
  
  start <- setNames(c(L = TRUE, R = FALSE)[side], NULL)
  flanking_intron <- flank(splicesites, 1, start = start)
  hits <- findOverlapsRanges(splicesites, frag_exonic, out="hits")
  
  if (include == "spliced" || include == "all") {
    
    frag_intron <- filterIntrons(frag_intron, frag_exonic, min_anchor)
    tmp <- findOverlapsRanges(flanking_intron, frag_intron, out="hits")
    hits_spliced <- intersect(hits, tmp)
    
  }
  
  if (include == "unspliced" || include == "all") {
    
    tmp <- findOverlapsRanges(flanking_intron, frag_exonic, out="hits")
    hits_unspliced <- intersect(hits, tmp)
    rm(tmp)
    qH <- queryHits(hits_unspliced)
    sH <- subjectHits(hits_unspliced)
    intron_anchor <- as(flank(splicesites, min_anchor, start = start),"IRangesList")
    exonic_anchor <- as(flank(splicesites, -min_anchor, start = start),"IRangesList")
    w_E <- sum(width(intersect(exonic_anchor[qH], frag_exonic[sH])))
    w_I <- sum(width(intersect(intron_anchor[qH], frag_exonic[sH])))
    hits_unspliced <- hits_unspliced[w_E == min_anchor & w_I == min_anchor]
    
  }
  
  if (include == "spliced") {
    hits <- hits_spliced
    rm(hits_spliced)
  } else if (include == "unspliced") {
    hits <- hits_unspliced
    rm(hits_unspliced)
    rm()
  } else if (include == "all") {
    hits <- union(hits_spliced, hits_unspliced)
    rm(hits_spliced)
    rm(hits_unspliced)
  }
  
  splicesite_index <- as.list(hits)
  rm(hits)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  
  if (counts) elementNROWS(splicesite_index)
  else splicesite_index
}


#' @rdname InternalFunctions
exonCompatible<- function(exons, spliceL, spliceR, frag_exonic,frag_intron, counts = TRUE){
  
  if (length(spliceL) == 1) spliceL <- rep(spliceL, length(exons))
  if (length(spliceR) == 1) spliceR <- rep(spliceR, length(exons))
  res<- c()
  
  lenSubject <- length(frag_exonic)
  lenQuery <- length(exons)
  
  subject_unlisted <- unlist(frag_exonic)
  subject_togroup <- togroup0(frag_exonic)
  query_unlisted <- exons
  query_togroup <- seq_along(exons)
  
  df <- cbind(start(subject_unlisted), end(subject_unlisted))
  options(digits=12)
  uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2),nrow = 2))
  
  posiciones_repetidas <- split(seq_along(uniqueVector), uniqueVector)[as.character(unique(uniqueVector))]
  
  new_subject_unlist <- df[!duplicated(uniqueVector),]
  rm(df)
  rm(uniqueVector)
  if (length(nrow(new_subject_unlist)) == 0){
    if (length(new_subject_unlist) == 2) {
      new_subject_unlist <- IRanges(start = new_subject_unlist[1],end = new_subject_unlist[2])
    }else{
      new_subject_unlist <- IRanges(start = NULL,end = NULL)
    }
  }else{
    new_subject_unlist <- IRanges(start = new_subject_unlist[,1],end = new_subject_unlist[,2])
  }
  
  hits_exonic <- list()
  positionsVector <- split(c(1:lenQuery), ceiling(seq_along(query_unlisted)/3000))
  for(val in c(1:length(positionsVector))){
    pos <- positionsVector[val]
    hits_listed <- as.list(findOverlaps(query_unlisted[unlist(pos)], new_subject_unlist, type = "any"))
    if (length(hits_listed) == 0) {
      hits_exonic <- hits_listed
      rm(hits_listed)
    }else{
      hits_exonic <- c(hits_exonic,hits_listed)
      rm(hits_listed)
    }
    if(val%%10 == 0){
      gc()
    }
  }
  
  gc()
  
  hits_introns <- findOverlapsRanges(exons, frag_intron)
  
  excl1 <- findOverlapsRanges(flank(exons, 1, TRUE), frag_exonic)
  excl1[which(!spliceL)] <- 0L
  
  excl2 <- findOverlapsRanges(flank(exons, 1, FALSE), frag_exonic)
  excl2[which(!spliceR)] <- 0L
  lastErase <- 1
  if (counts){
    res <- lapply(seq_along(excl1),function(i){
      eraseGroup <- c(hits_introns[[i]],excl1[[i]],excl2[[i]])
      lgroup <- length(setdiff(subject_togroup[unlist(posiciones_repetidas[hits_exonic[[i]]])],eraseGroup))
      if(i%%10000 == 0){
        excl1[c(lastErase:i)] <<- 0L
        excl2[c(lastErase:i)] <<- 0L
        hits_introns[c(lastErase:i)] <<- 0L
        lastErase <- i
        gc()
      }
      
      return(lgroup)
    })
    
  }else{
    res <- lapply(seq_along(excl1),function(i){
      eraseGroup <- c(hits_introns[[i]],excl1[[i]],excl2[[i]])
      lgroup <- setdiff(subject_togroup[unlist(posiciones_repetidas[hits_exonic[[i]]])],eraseGroup)
      res <- c(res, list(lgroup))
      if(i != 1 & i%%10000 == 0){
        excl1[c(lastErase:i)] <<- 0L
        excl2[c(lastErase:i)] <<- 0L
        hits_introns[c(lastErase:i)] <<- 0L
        lastErase <<- i
        gc()
      }
      return(lgroup)
    })
    
  }
  
  rm(excl1)
  rm(excl2)
  rm(hits_introns)
  gc()
  
  rm(new_subject_unlist)
  rm(subject_togroup)
  rm(query_unlisted)
  rm(subject_unlisted)
  rm(posiciones_repetidas)
  
  gc()
  
  return(res)
}

#' @rdname InternalFunctions
findOverlapsRanges <- function(query, subject, type = "any", out = "list")
{
  
  subject_hits <- c()
  query_hits <- c()
  lenSubject <- length(subject)
  lenQuery <- length(query)
  
  if (is(query, "IRangesList")) {
    query_unlisted <- unlist(query)
    query_togroup <- togroup0(query)
  }else {
    query_unlisted <- query
    query_togroup <- seq_along(query)
  }
  if (is(subject, "IRangesList")) {
    subject_unlisted <- unlist(subject)
    subject_togroup <- togroup0(subject)
    rm(subject)
  }else {
    subject_unlisted <- subject
    subject_togroup <- seq_along(subject)
    rm(subject)
  }
  hits <- c()
  gc()
  if (type == "equal") {
    hits_unlisted <- findMatches(query_unlisted, subject_unlisted)
    subject_hits <- subjectHits(hits_unlisted)
    query_hits <- queryHits(hits_unlisted)
    qH <- query_togroup[query_hits]
    sH <- subject_togroup[subject_hits]
    df <- cbind(qH, sH)
    uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2),nrow = 2))
    df <- as.data.frame(df)
    hits <- df[!duplicated(uniqueVector),]
    
    rm(uniqueVector)
    rm(df)
    rm(sH)
    rm(qH)
    rm(query_hits)
    rm(subject_hits)
    
    rm(hits_unlisted)
    gc(full=TRUE)
    
    hits <- Hits(as.integer(hits[, 1]), as.integer(hits[, 2]),lenQuery, lenSubject, sort.by.query = TRUE)
    
  }else {
    if (!is(query, "IRangesList")) {
      hits <- modFindOverlap(query_unlisted,subject_unlisted,subject_togroup,lenQuery, lenSubject, out)
    }else{
      hits_unlisted <- findOverlaps(query_unlisted, subject_unlisted, type = "any")
      subject_hits <- subjectHits(hits_unlisted)
      query_hits <- queryHits(hits_unlisted)
      qH <- query_togroup[query_hits]
      sH <- subject_togroup[subject_hits]
      df <- cbind(qH, sH)
      uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2)))
      df <- as.data.frame(df)
      hits <- df[!duplicated(uniqueVector),]
      rm(uniqueVector)
      rm(df)
      rm(sH)
      rm(qH)
      rm(query_hits)
      rm(subject_hits)
      
      rm(hits_unlisted)
      gc(full=TRUE)
      
      hits <- Hits(as.integer(hits[, 1]), as.integer(hits[, 2]),lenQuery, lenSubject, sort.by.query = TRUE)
    }
  }
  return(hits)
}

#' @rdname InternalFunctions
modFindOverlap <- function(query_unlisted,subject_unlisted,subject_togroup,lenQuery, lenSubject, out){
  
  
  df <- cbind(start(subject_unlisted), end(subject_unlisted))
  options(digits=12)
  uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2),nrow = 2))
  posiciones_repetidas <- split(seq_along(uniqueVector), uniqueVector)[as.character(unique(uniqueVector))]
  new_subject_unlist <- df[!duplicated(uniqueVector),]
  rm(df)
  rm(uniqueVector)
  if (length(nrow(new_subject_unlist)) == 0){
    if (length(new_subject_unlist) == 2) {
      new_subject_unlist <- IRanges(start = new_subject_unlist[1], end = new_subject_unlist[2])
    }else{
      new_subject_unlist <- IRanges(start = NULL, end = NULL)
    } 
  }else{
    new_subject_unlist <- IRanges(start = new_subject_unlist[,1], end = new_subject_unlist[,2])
  }
  
  hits_listed_res <- list()
  positionsVector <- split(c(1:lenQuery), ceiling(seq_along(query_unlisted)/3000))
  for(val in c(1:length(positionsVector))){
    pos <- positionsVector[val]
    hits_listed <- as.list(findOverlaps(query_unlisted[unlist(pos)], new_subject_unlist, type = "any"))
    for (overlap in c(1:length(hits_listed))) {
      hits_listed[[overlap]] <- unique(subject_togroup[unlist(posiciones_repetidas[hits_listed[[overlap]]])])
    }
    if (length(hits_listed) == 0) {
      hits_listed_res <- hits_listed
      rm(hits_listed)
    }else{
      hits_listed_res <- c(hits_listed_res,hits_listed)
      rm(hits_listed)  
    }
    if(val%%10 == 0){
      gc()
    }
  }
  rm(posiciones_repetidas)
  gc()
  if (out != "list") {
    hits_listed_res <- Hits(rep(c(1:length(query_unlisted)),lengths(hits_listed_res)),unlist(hits_listed_res),lenQuery, lenSubject, sort.by.query = FALSE)
  }
  
  rm(new_subject_unlist)
  rm(subject_togroup)
  rm(query_unlisted)
  rm(subject_unlisted)
  
  gc()
  return(hits_listed_res)
}






#' @rdname InternalFunctions
splicesiteCounts <- function(x, frag_exonic, frag_intron, min_anchor,
                             option = c("junction", "exon"), include)
{
  
  option <- match.arg(option)
  
  N_L <- splicesiteOverlap(flank(x, -1, TRUE),
                           switch(option, junction = "R", exon = "L"),
                           frag_exonic, frag_intron, min_anchor, include)
  N_R <- splicesiteOverlap(flank(x, -1, FALSE),
                           switch(option, junction = "L", exon = "R"),
                           frag_exonic, frag_intron, min_anchor, include)
  N <- IntegerList(mapply(c, N_L, N_R, SIMPLIFY = FALSE))
  rm(frag_exonic)
  rm(frag_intron)
  rm(N_L)
  rm(N_R)
  gc(full=TRUE)
  return(N)
  
}

#' @rdname InternalFunctions
exonCoverage <- function(exons, exons_i_frag, frag_exonic)
{
  
  expanded_exon <- factor(togroup0(exons_i_frag), seq_along(exons))
  expanded_frag_exonic <- frag_exonic[unlist(exons_i_frag)]
  rm(exons_i_frag)
  rm(frag_exonic)
  
  expanded_exon <- expanded_exon[togroup0(expanded_frag_exonic)]
  expanded_frag_exonic <- unlist(expanded_frag_exonic)
  
  irl <- split(expanded_frag_exonic, expanded_exon)
  rm(expanded_exon)
  rm(expanded_frag_exonic)
  gc(full=TRUE)
  coverage <- coverage(irl, shift = -start(exons) + 1, width = width(exons))
  ## coverage() returns a SimpleRleList
  ## need RleList() to obtain a CompressedRleList
  coverage <- RleList(coverage)
  
  rm(exons)
  gc(full=TRUE)
  return(coverage)
  
}



# predictTxFeatures_SGSeqMod.R #################
#' @rdname InternalFunctions
predictTxFeatures <- function(sample_info,
                              which = NULL,
                              alpha = 2,
                              psi = 0,
                              beta = 0.2,
                              gamma = 0.2,
                              include_counts = FALSE,
                              retain_coverage = FALSE,
                              junctions_only = FALSE,
                              min_junction_count = NULL,
                              min_anchor = 6,
                              max_complexity = 20,
                              min_n_sample = 1,
                              min_overhang = NA,
                              verbose = FALSE,
                              cores = 1)
{
  SGSeq:::checkSampleInfo(sample_info)
  cat("Amo a ve")
  list_whichSamples <- c()
  if (is.null(which)) {
    which <- globalWhich(sample_info)
    list_whichSamples <- sampleWhichPredict(sample_info, alpha, min_junction_count, which, TRUE)
  } else {
    which <- expandUnstrandedRanges(which)
    list_whichSamples <- sampleWhichPredict(sample_info, alpha, min_junction_count, which, FALSE)
  }
  clF <- makeCluster(cores, outfile = "reportPrediction.txt")
  
  list_features <- parLapplyLB(
    cl = clF,
    X = list_whichSamples,
    fun = predictTxFeaturesTotal,
    psi = psi,
    beta = beta,
    gamma = gamma,
    min_anchor = min_anchor,
    include_counts = include_counts,
    retain_coverage = retain_coverage,
    junctions_only = junctions_only,
    max_complexity = max_complexity,
    verbose = verbose
  )
  on.exit(stopCluster(cl = clF))
  featuresPerSample <- c()
  for (sample in sample_info$file_bam) {
    lFeaturesSample <- lapply(list_features, function(x) {
      if (sample == x$fileBamName) {
        x$gr
      }
    })
    
    lFeaturesSample <- unlist(lFeaturesSample)
    lFeaturesSample <-lFeaturesSample[!vapply(lFeaturesSample, is.null, logical(1))]
    
    if (length(lFeaturesSample) == 0) {
      features <- TxFeatures()
      
    } else {
      features <- do.call(c, setNames(lFeaturesSample, NULL))
      features <- sort(features)
      features <- TxFeatures(features)
    }
    featuresPerSample <- c(featuresPerSample, features)
    
  }
  
  features <-
    mergeTxFeatures(featuresPerSample, min_n_sample = min_n_sample)
  
  if (!is.null(min_overhang)) {
    features <- processTerminalExons(features, min_overhang)
    
  }
  
  return(features)
  
}

#' @rdname InternalFunctions
expandUnstrandedRanges <- function(x)
{
  i <- which(strand(x) == "*")
  if (length(i) > 0) {
    additional <- x[i]
    strand(additional) <- "-"
    strand(x)[i] <- "+"
    x <- c(x, additional)
  }
  return(x)
}

#' @rdname InternalFunctions
validationParameters <-function(sample_info,alpha,min_junction_count,which) {
  file_bam = sample_info$file_bam
  paired_end = sample_info$paired_end
  read_length = sample_info$read_length
  frag_length = sample_info$frag_length
  lib_size = sample_info$lib_size
  sample_name = sample_info$sample_name
  
  if (is.null(min_junction_count) && is.null(alpha)) {
    stop("Need to provide min_junction_count or alpha")
  }
  if (is.null(min_junction_count) &&(is.null(read_length) || is.null(frag_length) || is.null(lib_size))) {
    stop("For use with alpha, need to provide read_length,frag_length and lib_size")
  }
  if (!is.null(which) && !is(which, "GRanges")) {
    stop("argument which must be NULL or a GRanges object")
  }
  if (is.null(min_junction_count)) {
    min_junction_count <- SGSeq:::convertFpkmToCount(alpha, paired_end,read_length, frag_length, lib_size)
  }
  
  return(list(min_junction_count = min_junction_count,
              file_bam = file_bam,
              paired_end = paired_end,
              sample_name = sample_name
  )
  )
}


#' @rdname InternalFunctions
globalWhich <- function(sample_info) {
  res <- c()
  for (i in c(1:length(sample_info$file_bam))) {
    file_bam <- seqinfo(BamFile(sample_info$file_bam[i]))
    sl <- rep(seqlevels(file_bam), 2)
    st <- rep(c("+", "-"), rep(length(file_bam), 2))
    which <- GRanges(sl, IRanges(1, seqlengths(file_bam)[sl]), st)
    if (length(res) == 0) {
      res <- which
    } else{
      res <- c(res, which)
    }
  }
  res <- reduce(res)
  return(res)
  
}

#' @rdname InternalFunctions
sampleWhichPredict <- function(sample_info,alpha,min_junction_count,which,novo) {
  res <- c()
  for (sample in c(1:length(sample_info$file_bam))) {
    valSample <- validationParameters(sample_info[sample, ],alpha, min_junction_count, which)
    paired_end = valSample$paired_end
    searchWhich <- which
    if (paired_end & novo) {
      searchWhich <- searchWhich[searchWhich@strand == "+"]
    }
    bam_index <- idxstatsBam(file = sample_info[sample, ]$file_bam)
    res <-c(res,lapply(X = split(searchWhich, seq_along(searchWhich)),
                       FUN = listCHRFilePredict,
                       valSample = valSample,
                       bam_index = bam_index
    )
    )
  }
  return(res)
}

#' @rdname InternalFunctions
listCHRFilePredict <- function(range, valSample, bam_index) {
  file_bam = valSample$file_bam
  paired_end = valSample$paired_end
  sample_name = valSample$sample_name
  as.character(range@seqnames)
  index_chr <- bam_index[which(bam_index$seqnames == as.character(range@seqnames)), ]
  coefComplex <- index_chr$mapped / index_chr$seqlength
  min_junction_count = valSample$min_junction_count
  addWhich <-list(file_bam,paired_end,sample_name,min_junction_count,range,coefComplex)
  names(addWhich) <-c("file_bam","paired_end","sample_name","min_junction_count","which","coefComplex")
  return(addWhich)
}

#' @rdname InternalFunctions
predictTxFeaturesTotal <- function(sWhich,file_bam, paired_end, which,
                                   min_junction_count, psi, beta, gamma, min_anchor, include_counts,
                                   retain_coverage, junctions_only, max_complexity, sample_name, verbose){
  which <- sWhich$which
  fileBamName <- sWhich$file_bam
  sample_name <- sWhich$sample_name
  seqlevel <- as.character(seqnames(which))
  
  file_bam <- sWhich$file_bam
  paired_end <- sWhich$paired_end
  
  min_junction_count <- sWhich$min_junction_count
  
  if (is(file_bam, "BamFile")) {
    
    si <- seqinfo(file_bam)
    
  } else {
    
    si <- seqinfo(BamFile(file_bam))
    
  }
  seqlevel <- as.character(seqnames(which))
  strand <- as.character(strand(which))
  
  if (paired_end) {
    pairGap <- readGapPair(file_bam, paired_end, which, sample_name, verbose)
    print(paste(sample_name,"read CHR ", seqlevel, Sys.time(), sep = " "))
    
    gc(full = TRUE)
    
    frag_exonic <- pairGap$gapPlus$frag_exonic
    frag_intron <- pairGap$gapPlus$frag_intron
    junctions_df <- pairGap$gapPlus$junctions_df
    grPlus <- constructPaired(frag_exonic, frag_intron, junctions_df, min_junction_count,
                              psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                              junctions_only, max_complexity, sample_name, seqlevel, "+", si)
    print(paste("    ",sample_name,seqlevel, "+",Sys.time(), sep = " "))
    
    frag_exonic <- pairGap$gapMinus$frag_exonic
    frag_intron <- pairGap$gapMinus$frag_intron
    junctions_df <- pairGap$gapMinus$junctions_df
    grMinus <- constructPaired(frag_exonic, frag_intron, junctions_df, min_junction_count,
                               psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                               junctions_only, max_complexity, sample_name, seqlevel, "-", si)
    print(paste("    ",sample_name,seqlevel, "-",Sys.time(), sep = " "))
    rm(frag_exonic)
    rm(frag_intron)
    rm(pairGap)
    gc()
    gr <- c(grPlus,grMinus)
  }else{
    gap <- readGAlignments(file_bam, paired_end, which, sample_name, verbose)
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron
    gr <- constructPaired(frag_exonic, frag_intron, min_junction_count,
                          psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                          junctions_only, max_complexity, sample_name, seqlevel, strand, si)
    rm(frag_exonic)
    rm(frag_intron)
    rm(gap)
  }
  
  gc(full=TRUE)
  if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))
  res <- list(fileBamName=fileBamName, gr=gr)
  
  return(res)
  
}
#' @rdname InternalFunctions
predictJunctions <- function(frag_exonic, frag_intron,
                             min_junction_count, psi, min_anchor, retain_coverage)
{
  junctions <- unique(unlist(frag_intron)) + 1
  if (length(junctions) == 0) { return() }
  mcols(junctions) <- DataFrame("type" = rep("J", length(junctions)))
  mcols(junctions)$N <- junctionCompatible(junctions, frag_exonic,frag_intron, min_anchor)
  
  junctions <- junctions[which(mcols(junctions)$N >= min_junction_count)]
  
  if (length(junctions) == 0) { return() }
  
  ## consider splice junctions and splicesites with
  ## counts >= psi * max(splicesite counts)
  
  ## Note left/right (L/R) nomenclature: for the LHS splice site,
  ## the spliced boundary is situated on the right, for the RHS
  ## splice site, the spliced boundary is situated on the left
  
  if (psi > 0 || retain_coverage) {
    mcols(junctions)$N_splicesite <- splicesiteCounts(junctions,frag_exonic, frag_intron, min_anchor, "junction", "all")
    junctions <- junctions[which(mcols(junctions)$N >=psi * max(mcols(junctions)$N_splicesite))]
    if (length(junctions) == 0) { return() }
  }
  if (!retain_coverage) {
    mcols(junctions)$N_splicesite <- NULL
  }
  junctions <- SGSeq:::completeMcols(junctions, retain_coverage)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(junctions)
}

#' @rdname InternalFunctions
readGapPair <- function(file, paired_end, which, sample_name, verbose)
{
  
  if (length(which) != 1) {
    stop("which argument must have length 1")
  }
  
  flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                      isPaired = T, isProperPair = T, 
                      isUnmappedQuery = F,
                      isNotPassingQualityControls=F)
  param <- ScanBamParam(flag = flag, tag = "XS", which = which)
  if (paired_end) {
    gap <- suppressWarnings(readGAlignmentPairs(file = file,param = param)) 
    gap <- SGSeq:::propagateXS(gap)
  }else{
    gap <- suppressWarnings(readGAlignments(file = file, param = param))
  }
  
  gap <- SGSeq:::filterGap(gap)
  
  mcols(gap)$strand <- SGSeq:::XS2strand(mcols(gap)$XS)
  
  gapPlus <- fragExonIntron(gap = gap, strand="+",verbose)
  gapMinus <- fragExonIntron(gap = gap, strand="-",verbose)
  rm(gap)
  gc()
  return(list(gapPlus=gapPlus,gapMinus=gapMinus))
  
}
#' @rdname InternalFunctions
generateWarningMessage <- function (fun_name, item, msg) 
{
  message(makeWarningMessage(fun_name, item, msg))
}

#' @rdname InternalFunctions
fragExonIntron <- function(gap,strand,verbose){
  gap <- gap[mcols(gap)$strand %in% c(strand, "*")]
  frag_exonic <- reduce(ranges(grglist(gap, drop.D.ranges = TRUE)))
  frag_intron <- IRangesList(rep(list(IRanges()), length(gap)))
  pos_junctions <- which((njunc(gap@first)<=1 | njunc(gap@last)<=1))
  junctionsRanges <- ranges(junctions(gap[pos_junctions]))
  frag_intron[pos_junctions] <- junctionsRanges
  junctions_df <- NULL
  diff <- setdiff(frag_exonic, frag_intron)
  excl <- which(sum(width(frag_exonic)) > sum(width(diff)))
  
  if (length(excl) > 0) {
    
    if (verbose) {
      
      msg <- paste("filtered",length(excl),"inconsistent paired alignments in",gr2co(which))
      SGSeq:::generateWarningMessage("readGap",sample_name,msg)
      
    }
    
    frag_exonic <- frag_exonic[-excl]
    frag_intron <- frag_intron[-excl]
    
  }
  
  gap <- list(frag_exonic = frag_exonic, frag_intron = frag_intron, junctions_df = junctions_df)
  rm(frag_exonic)
  rm(frag_intron)
  return(gap)
}

#' @rdname InternalFunctions
togroup0 <- S4Vectors:::quick_togroup

#' @rdname InternalFunctions
constructPaired <- function(frag_exonic, frag_intron, junctions_df, min_junction_count,
                            psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                            junctions_only, max_complexity, sample_name, seqlevel, strand, si){
  
  if (length(frag_exonic) == 0) {
    
    gr <- NULL
    
  } else {
    ir <- predictSpliced(frag_exonic, frag_intron, junctions_df, min_junction_count,
                         psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                         junctions_only, max_complexity, sample_name, seqlevel, strand)
    if (is.null(ir)) {
      gr <- NULL
      
    } else {
      gr <- constructGRangesFromRanges(ir, seqlevel, strand, si)
    }
    
  }
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(gr)
}




#' @rdname InternalFunctions
co2str <- function (seqlevel, start, end, strand) 
{
  paste0(seqlevel, ":", start, "-", end, ":", strand)
}

#' @rdname InternalFunctions
predictSpliced <- function(frag_exonic, frag_intron, junctions_df, min_junction_count,
                           psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                           junctions_only, max_complexity, sample_name, seqlevel, strand)
{
  junctions <- predictJunctions(frag_exonic, frag_intron,min_junction_count, psi, min_anchor, retain_coverage)
  
  if (is.null(junctions)) { return() }
  
  if (!junctions_only) {
    lower <- max(min_junction_count * beta, 1)
    frag_coverage <- coverage(unlist(frag_exonic))
    islands <- slice(frag_coverage, lower, rangesOnly = TRUE)
    
    ## skip problematic regions
    if (!is.na(max_complexity)) {
      ir <- as(slice(coverage(junctions), max_complexity), "IRanges")
      if (length(ir) > 0) {
        
        junctions_stripped <- junctions
        mcols(junctions_stripped) <- NULL
        loci <- reduce(c(junctions_stripped, islands))
        excl <- loci[loci %over% ir]
        
        junctions <- junctions[!junctions %over% excl]
        islands <- islands[!islands %over% excl]
        
        excl_str <- co2str(seqlevel, start(excl), end(excl), strand)
        
        SGSeq:::generateWarningMessage("predictSpliced",sample_name,paste("skipping", excl_str))
        
        if (length(junctions) == 0) { return() }
        
      }
      
    }
    
    features <- junctions
    splicesites_L <- extractSplicesitesFromJunctions(junctions, "L")
    splicesites_R <- extractSplicesitesFromJunctions(junctions, "R")
    splicesites <- c(splicesites_L, splicesites_R)
    candidates <- predictCandidatesInternal(islands, splicesites,frag_coverage, beta)
    exons_I <- predictExonsInternal(candidates, frag_exonic,frag_intron, beta, min_anchor, include_counts, retain_coverage)
    if (!is.null(exons_I)) { features <- c(features, exons_I) }
    candidates <- predictCandidatesTerminal(islands, splicesites, "exon_L")
    exons_L <- predictExonsTerminal(candidates, frag_exonic, frag_intron,gamma, min_anchor, "exon_L", include_counts, retain_coverage)
    if (!is.null(exons_L)) { features <- c(features, exons_L) }
    candidates <- predictCandidatesTerminal(islands, splicesites, "exon_R")
    exons_R <- predictExonsTerminal(candidates, frag_exonic, frag_intron,gamma, min_anchor, "exon_R", include_counts, retain_coverage)
    if (!is.null(exons_R)) { features <- c(features, exons_R) }
  } else {
    features <- junctions
  }
  if (!include_counts) {
    mcols(features)$N <- NULL
  }
  rm(frag_exonic)
  rm(frag_intron)
  rm(junctions)
  gc(full=TRUE)
  features <- sort(features)
  
  return(features)
  
}

#' @rdname InternalFunctions
extractSplicesitesFromJunctions <- function(junctions, type = c("L", "R")){
  
  type <- match.arg(type)
  S <- flank(junctions, -1, switch(type, "L" = TRUE, "R" = FALSE))
  S_pos <- as.character(start(S))
  pos_N <- tapply(mcols(junctions)$N, S_pos, sum)
  pos_N <- setNames(as.integer(pos_N), names(pos_N))
  i <- which(!duplicated(S_pos))
  S <- S[i]
  S_pos <- S_pos[i]
  mcols(S) <- DataFrame(type = rep(type, length(S)),N = pos_N[match(S_pos, names(pos_N))])
  return(S)
  
}

#' @rdname InternalFunctions
predictExonsTerminal <- function(candidates, frag_exonic, frag_intron, relCov,
                                 min_anchor, type = c("exon_L", "exon_R"), include_counts, retain_coverage){
  type <- match.arg(type)
  if (length(candidates) == 0) { return() }
  
  spliceL <- switch(type, "exon_L" = FALSE, "exon_R" = TRUE)
  spliceR <- switch(type, "exon_L" = TRUE, "exon_R" = FALSE)
  
  index <- exonCompatible(candidates, spliceL, spliceR,frag_exonic, frag_intron, FALSE)
  coverage <- exonCoverage(candidates, index, frag_exonic)
  
  splicesite <- flank(candidates, -1,start = switch(type, "exon_L" = FALSE, "exon_R" = TRUE))
  N_splicesite <- splicesiteOverlap(splicesite,
                                    switch(type, "exon_L" = "R", "exon_R" = "L"),
                                    frag_exonic, frag_intron, min_anchor, "spliced")
  
  ranges <- (coverage >= relCov * N_splicesite)
  el <- elementNROWS(ranges)
  rl <- runLength(ranges)
  
  if (type == "exon_L") {
    ir <- IRanges(end = el, width = SGSeq:::plast(rl))
  }
  if (type == "exon_R") {
    ir <- IRanges(start = 1, width = SGSeq:::pfirst(rl))
  }
  if (length(ir) ==    0) { return() }
  exons <- SummarizedExperiment::shift(ir, start(candidates) - 1)
  mcols(exons) <- DataFrame(type = rep(type, length(exons)))
  
  if (include_counts) {
    mcols(exons)$N <- exonCompatible(exons, spliceL, spliceR,frag_exonic, frag_intron)
  }
  if (retain_coverage) {
    mcols(exons)$N_splicesite <- as(N_splicesite, "IntegerList")
    mcols(exons)$coverage <- coverage[setNames(split(ir, seq_along(ir)), NULL)]
  }
  
  exons <- SGSeq:::completeMcols(exons, retain_coverage)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(exons)
  
}
#' @rdname InternalFunctions
predictCandidatesTerminal <- function(islands, splicesites,type = c("exon_L", "exon_R"))
{
  
  type <- match.arg(type)
  splicesites <- splicesites[mcols(splicesites)$type ==switch(type, "exon_L" = "L", "exon_R" = "R")]
  hits <- findOverlaps(splicesites, islands)
  spliced_boundary <- splicesites[queryHits(hits)]
  island <- islands[subjectHits(hits)]
  
  if (type == "exon_L") {
    candidates <- IRanges(start(island), end(spliced_boundary))
  }
  if (type == "exon_R") {
    candidates <- IRanges(start(spliced_boundary), end(island))
  }
  mcols(candidates) <- DataFrame(N = mcols(splicesites)$N[queryHits(hits)])
  return(candidates)
}
#' @rdname InternalFunctions
predictExonsInternal <- function(candidates, frag_exonic, frag_intron, relCov,min_anchor, include_counts, retain_coverage)
{
  
  if (length(candidates) == 0) { return() }
  candidate_index <- exonCompatible(candidates, TRUE, TRUE,frag_exonic, frag_intron, FALSE)
  candidate_coverage <- exonCoverage(candidates, candidate_index,frag_exonic)
  candidate_N_splicesite <- splicesiteCounts(candidates, frag_exonic,frag_intron, min_anchor, "exon", "spliced")
  index <- which(min(candidate_coverage) >=relCov * min(candidate_N_splicesite))
  if (length(index) == 0) { return() }
  exons <- candidates[index]
  mcols(exons) <- DataFrame("type" = rep("I", length(exons)))
  
  if (include_counts) {
    mcols(exons)$N <- exonCompatible(exons, TRUE, TRUE, frag_exonic,frag_intron)
  }
  if (retain_coverage) {
    mcols(exons)$N_splicesite <- candidate_N_splicesite[index]
    mcols(exons)$coverage <- candidate_coverage[index]
  }
  
  exons <- SGSeq:::completeMcols(exons, retain_coverage)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(exons)
  
}
#' @rdname InternalFunctions
predictCandidatesInternal <- function(islands, splicesites, frag_coverage,
                                      relCov)
{
  
  ## for each island, identify overlapping splice sites
  island_splicesite <- as.list(findOverlaps(islands, splicesites))
  
  ## for each island, obtain all pairs of overlapping splice sites
  island_splicesite_pairs <- mapply(expand.grid,island_splicesite, island_splicesite, SIMPLIFY = FALSE)
  splicesite_pairs <- unique(do.call(SGSeq:::rbindDfsWithoutRowNames,island_splicesite_pairs))
  
  ## retain pairs of splice sites that are consistent with
  ## flanking an internal exon
  N_1 <- mcols(splicesites)$N[splicesite_pairs[, 1]]
  N_2 <- mcols(splicesites)$N[splicesite_pairs[, 2]]
  type_1 <- mcols(splicesites)$type[splicesite_pairs[, 1]]
  type_2 <- mcols(splicesites)$type[splicesite_pairs[, 2]]
  pos_1 <- start(splicesites)[splicesite_pairs[, 1]]
  pos_2 <- start(splicesites)[splicesite_pairs[, 2]]
  i <- which(type_1 == "R" & type_2 == "L" & pos_1 <= pos_2)
  candidates <- IRanges(pos_1[i], pos_2[i])
  mcols(candidates) <- DataFrame(N = IntegerList(mapply(c, N_1[i], N_2[i], SIMPLIFY = FALSE)))
  
  if (length(candidates) > 0) {
    
    ## retain candidate internal exons with sufficient read coverage
    candidates_frag_coverage <- split(frag_coverage[candidates],togroup0(candidates))
    i <- which(min(candidates_frag_coverage) >=relCov * min(mcols(candidates)$N))
    candidates <- candidates[i]
    
  }
  gc(full=TRUE)
  return(candidates)
  
}

#' @rdname InternalFunctions
constructGRangesFromRanges <- function(x, seqname, strand, seqinfo)
{
  
  x_mcols <- mcols(x)
  mcols(x) <- NULL
  
  if (strand == "+") {
    x_mcols_type <- as.character(x_mcols$type)
    x_mcols_type <- sub("exon_L", "F", x_mcols_type, fixed = TRUE)
    x_mcols_type <- sub("exon_R", "L", x_mcols_type, fixed = TRUE)
    x_mcols$type <- factor(x_mcols_type,levels = c("J", "I", "F", "L", "U"))
  }
  if (strand == "-") {
    x_mcols_type <- as.character(x_mcols$type)
    x_mcols_type <- sub("exon_L", "L", x_mcols_type, fixed = TRUE)
    x_mcols_type <- sub("exon_R", "F", x_mcols_type, fixed = TRUE)
    x_mcols$type <- factor(x_mcols_type,levels = c("J", "I", "F", "L", "U"))
    if ("N_splicesite" %in% names(x_mcols)) {
      x_mcols$N_splicesite <- endoapply(x_mcols$N_splicesite, rev)
    }
    if ("coverage" %in% names(x_mcols)) {
      x_mcols$coverage <- endoapply(x_mcols$coverage, rev)
    }
  }
  
  gr <- GRanges(rep(seqname, length(x)), x, rep(strand, length(x)),x_mcols, seqinfo = seqinfo)
  
  rm(x)
  gc(full=TRUE)
  
  return(gr)
  
}



# PrepareBam_EP_mod.R #################


#' @rdname InternalFunctions
AnnEventsFunc <- function(EventsDetection_pred, EventsDetection_ann, cores){
  
  registerDoParallel(cores=cores)
  listEventsPred <- foreach(event=unlist(EventsDetection_pred, recursive = F), .packages = 'GenomicRanges') %dopar% {
    P1 <- GRanges(event$P1)
    P2 <- GRanges(event$P2)
    Ref <- GRanges(event$Ref)
    GRangesList(c(P1,P2,Ref))[[1]]
  }
  listEventsPred <- GRangesList(listEventsPred)
  listEventsAnn <- foreach(event=unlist(EventsDetection_ann, recursive = F), .packages = 'GenomicRanges') %dopar% {
    P1 <- GRanges(event$P1)
    P2 <- GRanges(event$P2)
    Ref <- GRanges(event$Ref)
    GRangesList(c(P1,P2,Ref))[[1]]
  }
  stopImplicitCluster()
  listEventsAnn <- GRangesList(listEventsAnn)
  fso <- GenomicAlignments:::findSpliceOverlaps(listEventsPred,listEventsAnn, type="equal")
  fso <- fso[which(fso@elementMetadata$compatible)]
  fso <- fso[which(lengths(listEventsAnn[fso@to])==lengths(listEventsPred[fso@from]))]
  percentIdentity <- 1-sum(abs(width(listEventsAnn[fso@to])-width(listEventsPred[fso@from])))/max(sum(width(listEventsAnn[fso@to])), sum(width(listEventsPred[fso@from])))
  fso <- fso[which(percentIdentity >=.98)]
  
  countEvent <- 1
  for (gene in c(1:length(EventsDetection_pred))) {
    for(posEvent in c(1:length(EventsDetection_pred[[gene]]))){
      event <- EventsDetection_pred[[gene]][[posEvent]]
      if(countEvent %in% fso@from){
        event$Info$Ann <- TRUE
      } else{
        event$Info$Ann <- FALSE
      }
      EventsDetection_pred[[gene]][[posEvent]] <- event
      countEvent <- countEvent+1 
    }
    
  }
  return(EventsDetection_pred)
}

#' @rdname InternalFunctions
getBamInfo <- function(PathSamplesAbundance, region, cores = 1)
{
  sample_name <- dir(PathSamplesAbundance,pattern = "*.bam$")
  if(identical(sample_name,character(0))){
    file_bam <- c()
    for(folder in list.dirs(PathSamplesAbundance)){
      files_bam_folder <- list.files(folder,pattern = "*.bam$")
      sample_name <- c(sample_name,files_bam_folder)
      file_bam <- c(file_bam, rep(folder, length(files_bam_folder)))
      
    }
    sample_info <- data.frame(sample_name = sample_name,
                              file_bam = file_bam, stringsAsFactors = FALSE)
  }else{
    sample_info <- data.frame(sample_name = sample_name,
                              file_bam = rep(PathSamplesAbundance, length(sample_name)), stringsAsFactors = FALSE)
  }
  
  SGSeq:::checkSampleInfo(sample_info, FALSE)
  listSamples <- split(sample_info, c(1:nrow(sample_info)))
  clF <- makePSOCKcluster(cores, outfile = "bamInfoSamples.txt")
  list_bamInfo <- clusterApplyLB(cl = clF,
                                 x = listSamples,
                                 fun = getBamInfoPerSample,
                                 region=region)
  on.exit(stopCluster(cl = clF))
  SGSeq:::checkApplyResultsForErrors(list_bamInfo,
                                     "getBamInfoPerSample",
                                     sample_info$sample_name,
                                     "character")
  
  bamInfo <- do.call(SGSeq:::rbindDfsWithoutRowNames, list_bamInfo)
  
  SGSeq:::checkBamInfo(bamInfo)
  
  cols <- c("paired_end", "read_length", "frag_length", "lib_size")
  cols <- cols[cols %in% names(bamInfo)]
  
  for (col in cols) {
    sample_info[[col]] <- bamInfo[[col]]
    
  }
  sample_info$file_bam <- paste0(sample_info$file_bam,"/",sample_info$sample_name)
  return(sample_info)
  
}


#' @rdname InternalFunctions
getBamInfoPerSample <- function(sample_info, region)
{
  file_bam <- paste0(sample_info$file_bam,"/",sample_info$sample_name)
  sample_name <- sample_info$sample_name
  if (is(file_bam, "BamFile")) {
    file_tmp <- file_bam
    
  } else {
    file_tmp <- BamFile(file_bam)
    
  }
  
  si <- seqinfo(file_tmp)
  sl <- rep(seqlevels(si),2)
  st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
  which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
  
  flag <-
    scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      hasUnmappedMate = FALSE,
      isDuplicate = FALSE
      
    )
  what <- c("qname", "flag", "qwidth", "isize")
  param <-
    ScanBamParam(
      flag = flag,
      what = what,
      which = which[region[1]],
      tag = "XS"
    )
  
  bam <- scanBam(file = file_tmp, param = param)[[1]]
  gc()
  XS <- !is.null(bam$tag$XS)
  paired_end <- any(bamFlagTest(bam$flag, "isPaired"))
  read_length <- median(bam$qwidth, na.rm = TRUE)
  x <- data.frame(
    XS = XS,
    paired_end = paired_end,
    read_length = read_length,
    stringsAsFactors = FALSE
  )
  if (paired_end) {
    isize <- bam$isize
    frag_length <- median(isize[which(isize > 0)], na.rm = TRUE)
    infoBam <- idxstatsBam(file = file_bam, index = file_bam)
    x$lib_size <- sum(infoBam$mapped)/2
    
  } else {
    infoBam <- idxstatsBam(file = file_bam, index = file_bam)
    desv <- (infoBam$mapped[1]/length(unique(bam$qname))-1)
    x$lib_size <- sum(infoBam$mapped) - (sum(infoBam$mapped) * desv)
    frag_length <- NA_real_
    
  }
  rm(bam)
  gc()
  x$frag_length <- frag_length
  
  
  message(sample_name)
  
  return(x)
  
}
