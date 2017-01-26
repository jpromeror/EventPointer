# Creation of CDF file from .flat file
#
#' Internal function called by EvenPointer
#'
#' Function to create the CDF file given the corresponding .flat file
#'
#' @keywords internal
#'
#' @param file .flat file
#' @param chipType name of the array
#' @param tags tags to identify cdf file
#' @param rows rows in the array
#' @param cols columns in the array
#' @param verbose output to user
#' @param xynames identifiers for chromosome X and Y
#' @param gcol Internal parameter
#' @param ucol Internal parameter
#' @param splitn Internal parameter
#' @param col.class Internal parameter
#'
#'
#' @return CDF file
#'

####################################################################
# flat2Cdf
#---------
# this function takes a "flat" file and converts it to a binary CDF file
# it was downloaded from: http://www.aroma-project.org/howtos/create_CDF_from_scratch/
# for further details see that link. 
#
#example: flat2Cdf(file="hjay.r1.flat",chipType="hjay",tag="r1,TC")
#file: assumes header...better perhaps to have ... that passes to read.table?; requires header X, Y
#ucol: unit column
#gcol: group column
#col.class: column classes of file (see read.table); NOTE: needs check that right number?
#splitn: parameter that controls the number of initial chunks that are unwrapped (number of characters of unit names used to keep units together for initial chunks)
#rows: 
#cols:
####################################################################
flat2Cdf<-function(file,chipType,tags=NULL,rows=2560,cols=2560,verbose=10,xynames=c("X","Y"),
                   gcol=5,ucol=6,splitn=4,col.class=c("integer","character")[c(1,1,1,2,2,2)],...) {
  split.quick<- 
    function(r,ucol,splitn=3,verbose=TRUE) {
      rn3<-substr(r[,ucol],1,splitn)
      split.matrix<-split.data.frame
      rr<-split(r,factor(rn3))
      if (verbose) cat(" split into",length(rr),"initial chunks ...")
      rr<-unlist(lapply(rr,FUN=function(u) split(u,u[,ucol])),recursive=FALSE)
      if (verbose) cat(" unwrapped into",length(rr),"chunks ...")
      names(rr)<-substr(names(rr),splitn+2,nchar(rr))
      rr
      ##    rr<-unlist(lapply(rr,FUN=function(u) split(u,u[,ucol])),recursive=FALSE,use.names=FALSE)
      ##    namrr<-sapply(rr,function(u){nam<-unique(u[,ucol]); if(length(nam)>1) stop("Programming Error (units", nam,"). Please report") else return(nam)},USE.NAMES=FALSE)
      ##    names(rr)<-namrr
      ####    rr<-unlist(lapply(rr,FUN=function(u) split(u[,-ucol],u[,ucol])),recursive=FALSE)
      ####    names(rr)<-sapply(strsplit(names(rr),"\\."),.subset2,2)
      
    }
  
  if (verbose) cat("Reading TXT file ...")
  file<-read.table(file,header=TRUE,colClasses=col.class,stringsAsFactors=FALSE,comment.char="",...)
  if (verbose) cat(" Done.\n")
  
  if (verbose) cat("Splitting TXT file into units ...")
  gxys<-split.quick(file,ucol,splitn)
  rm(file); gc()
  if (verbose) cat(" Done.\n")
  
  l<-vector("list",length(gxys))
  if (verbose) cat("Creating structure for",length(gxys),"units (dot=250):\n")
  for(i in  1:length(gxys)) {
    sp<-split(gxys[[i]],factor(gxys[[i]][,gcol]))
    e<-vector("list",length(sp))
    for(j in 1:length(sp)) {
      np<-nrow(sp[[j]])
      e[[j]]<-list(x=sp[[j]][,xynames[1]],y=sp[[j]][,xynames[2]],pbase=rep("A",np),tbase=rep("T",np),atom=0:(np-1),indexpos=0:(np-1),
                   groupdirection="sense",natoms=np,ncellsperatom=1)
    }
    names(e)<-names(sp)
    l[[i]]<-list(unittype=1,unitdirection=1,groups=e,natoms=nrow(gxys[[i]]),ncells=nrow(gxys[[i]]),ncellsperatom=1,unitnumber=i)
    if (verbose) { if(i %% 250==0) cat("."); if(i %% 5000==0) cat("(",i,")\n",sep="") }
  }
  cat("\n")
  names(l)<-names(gxys)
  if(!is.null(tags) && tags!="") filename<-paste(chipType,tags,sep=",")
  else filename<-chipType
  filename<-paste(filename,"cdf",sep=".")
  hdr<-list(probesets=length(l),qcprobesets=0,reference="",chiptype=chipType,filename=filename,
            nqcunits=0,nunits=length(l),rows=rows,cols=cols,refseq="",nrows=rows,ncols=cols)
  writeCdf(hdr$filename, cdfheader=hdr, cdf=l, cdfqc=NULL, overwrite=TRUE, verbose=verbose)
  invisible(list(cdfList=l,cdfHeader=hdr))
}
