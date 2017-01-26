# GTF file creation
#
#' Internal function called by EvenPointer
#'
#' Function to write the GTF file in the requried format to be loaded into IGV
#'
#' @keywords internal
#'
#' @param PATH Directory where to write the .gtf file
#' @param Data Basic information of the event (Gene, Type,...,etc)
#' @param Probes Probes mapped to any of the paths of an alternative splicing event
#' @param Paths Paths that represent an alternative splicing event
#'
#'
#' @return Adds lines to the corresponding .GTF file
#'


##############################################################################
# function to create Event plots in IGV
##############################################################################
WriteGTF <- function(PATH,Data,Probes,Paths){
  
  
  
  STRAND <- unique(Paths[,5])
  FILE.probes <- paste(PATH,"/probes.gtf",sep="")
  ### Error en los grep.. si no encuentra regresa 0 y no NA
  PATHS <- as.vector(unique(Probes[,6]))
  
  
  for(i in 1:nrow(Probes))
  {
    
    if (STRAND=="+"){
      START<-Probes[i,3]
      END <- Probes[i,3]+Probes[i,4]
    }else if (STRAND=="-"){
      START<-Probes[i,3]-Probes[i,4]
      END <- Probes[i,3]
    }
    
    # browser()
    if(Probes[i,6]=="Ref")
    {
      COL <- "#B0B0B0"
    }else if(Probes[i,6]=="Path1")
    {
      COL <- "#D00000"
    }else if(Probes[i,6]=="Path2")
    {
      COL <- "#00CC33"
    }
    PROBES <- paste(Probes[i,2],"\t","microarray","\t","probe","\t",Probes[i,3],
                    "\t",END,"\t","0","\t","*","\t","0","\t",
                    "event=",Data[1,4],"; path=",Probes[i,6],"; color=",COL,";","probeID=",Probes[i,1],";",sep="")
    cat(file=FILE.probes,PROBES,sep="\n",append=TRUE)
    
  }
  
  
  
  # Paths GTF
  
  II <- order(Paths[,7])
  Paths <- Paths[II,]
  PATHS <- unique(Paths[,7])
  FILE.paths <- paste(PATH,"/paths.gtf",sep="")
  
  
  GENESSSSS <- paste(Paths[,1],"\t","EventPointer","\t","gene","\t",min(Paths[,2]),
                     "\t",max(Paths[,3]),"\t","0","\t",Paths[,5],"\t","0","\t",
                     paste("gene_id ",Data[1,2],"_",Data[1,3],"; ",
                           "type ",shQuote(as.vector(Data[1,4]),type="cmd"),
                           "; color=","#000000",";",sep=""),sep="")
  GENESSSSS <- unique(GENESSSSS)
  
  cat(file=FILE.paths,GENESSSSS,sep="\n",append=TRUE)
  
  
  # browser()
  for (i in 1:length(PATHS)){
    ii <- which(Paths[,7]==PATHS[i])
    # if (all(!is.na(grep("Ref",PATHS[i])))){
    if (length(!is.na(grep("Ref",PATHS[i])))!=0){
      COL <- "#B0B0B0"
      aaaaaa <- 3
    }
    # if (all(!is.na(match("A",PATHS[i])))){
    if (length(!is.na(grep("A",PATHS[i])))!=0){
      COL <- "#D00000"
      aaaaaa <- 2
    }
    # if (all(!is.na(match("B",PATHS[i])))){
    if (length(!is.na(grep("B",PATHS[i])))!=0){
      COL <- "#00CC33"
      aaaaaa <- 1
    }
    # if (all(!is.na(match("Empty",PATHS[i])))){
    if (length(!is.na(grep("Empty",PATHS[i])))!=0){
      COL <- "#FFFFFF"
    }
    
    
    # browser()
    #TRANS <- paste(Paths[ii,1],"\t","EventPointer","\t","transcript","\t",min(Paths[,2])-10*as.numeric(as.matrix(Data[2]))-aaaaaa
    TRANS <- paste(Paths[ii,1],"\t","EventPointer","\t","transcript","\t",min(Paths[ii,2]),#MINGENE-as.numeric(as.matrix(Data[2])),
                   "\t",max(Paths[ii,3]),"\t","0","\t",Paths[ii,5],"\t","0","\t",
                   paste("gene_id ",Data[1,2],"_",Data[1,3],"; ",
                         "transcript_id ",shQuote(Paths[ii,7],type="cmd"),"_",gsub(" ","_",Data[1,4]),"_",unique(Paths[ii,6]),"_",Data[1,3],"; ",
                         "type ",shQuote(Data[1,4],type="cmd"),
                         "; color=",COL,";",sep=""),sep="")
    TRANS <- unique(TRANS)
    
    GTF <- paste(Paths[ii,1],"\t","EventPointer","\t","exon","\t",Paths[ii,2],
                 "\t",Paths[ii,3],"\t","0","\t",Paths[ii,5],"\t","0","\t",
                 paste("gene_id ",Data[1,2],"_",Data[1,3],"; ",
                       "transcript_id ",shQuote(Paths[ii,7],type="cmd"),"_",gsub(" ","_",Data[1,4]),"_",unique(Paths[ii,6]),"_",Data[1,3],"; ",
                       "type ",shQuote(Data[1,4],type="cmd"),
                       "; exon_number ",1:length(ii),"; color=",COL,";",sep=""),sep="")
    #if (i == 1){
    #  cat(file=FILE.paths,TRANS,sep="\n")
    #}else{
    cat(file=FILE.paths,TRANS,sep="\n",append=TRUE)
    #}
    cat(file=FILE.paths,GTF,sep="\n",append=TRUE)
  }
  
  
}
