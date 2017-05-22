library("seqinr") ## required for read.fasta()
library("plyr") ## required for ddply()

discardDups <- function(str){
  q <- 1
  while (q <= length(str)){
    str <- str[which(upsidedown(str) != str[q])]
    q <- q+1
  }
  return(str)
}

upsidedown <- function(bp){
  # "a/b > c/d" becomes "b/a > d/c"
  paste0(  substr(bp,3,3),"/",
           substr(bp,1,1)," > ",
           substr(bp,9,9),"/",
           substr(bp,7,7)
  )
  
}

cutoff <- 2 ## varient is called parental is it shows up in 'cutoff' or more strains 

## read in data, with mild processing
VCFfolder <- "Z:/R/ENU analysis/VCF/"
VCFfiles <- list.files(VCFfolder)
nfiles <- length(VCFfiles)
allCPRC = NULL
for(i in 1:nfiles){
  print( paste0( c("Reading... ",VCFfiles[i]),collapse="") )
  temp <- read.table(paste0(c(VCFfolder,VCFfiles[i]),collapse=""), stringsAsFactors=FALSE)    ## read in
  temp$CPRC <- apply(temp,1,function(x) paste(x[1],x[2],x[4],x[5]))                           ## add CPRC col
  allCPRC <- c(allCPRC,levels(factor(temp$CPRC)))                                             ## add these CPRCs to growing list of all CPRCs
  name <- paste0(c("VCF",i),collapse="")                                                      ## build variable name
  assign(name, temp)                                                                          ## assign var
  }
print( paste0( c("Read ",nfiles," files."),collapse="") )


## subtract parental
CPRC <- data.frame("CPRC" = allCPRC)                          ## move CPRC into a df, so ddply() can work on it
CPRCcount <- ddply(CPRC,.(CPRC),nrow)                         ## counts times each CPRC shows up, and puts that in a new column
colnames(CPRCcount)[2] <- "numWorms"                          ## renames the new column to something sensible

parental <- CPRCcount[which(CPRCcount$numWorm>=cutoff),]      ## "parental" variants must show up >= cutoff number of times

for(i in 1:nfiles){
  name1 <- paste0(c("VCF",i),collapse="")                              ## build var name to call
  temp1 <- get(name1)                                                  ## call var
  noParents <- temp1[which(is.na(match(temp1$CPRC,parental$CPRC))),]   ## remove variants that match a parental
  assign(name, noParents)                                              ## reassign 
}
  

## concatenate
vcf <- NULL

for(i in 1:nfiles){
  name2 <- paste0(c("VCF",i),collapse="")
  temp2 <- get(name2)
  vcf <- rbind(vcf,temp2)
}

vcf <- vcf[which((nchar(vcf$V4)==1) & nchar(vcf$V5)==1),]  ## drop indels
colnames(vcf)[c(1,2,4,5)] <- c("chr","pos","wt_dna","mut_dna")
vcf$chr <- apply(vcf,1,function(x) substring(x[1],4))

## cleanup the worksapce
rm(list=setdiff(ls(), c("vcf","discardDups","upsidedown"))) 
gc()



## define some basic stuff
purines <- c("A","G")
pyrimidines <- c("T","C")
nt <- c(purines,pyrimidines)
ntnum <- 1:length(nt)
GC <- c("G/C","C/G")
AT <- c("A/T","T/A")
allPairs <- c(GC,AT)
changeRows <- NULL
numPairs <- 1:length(allPairs)
for(i in numPairs){
  for(j in numPairs[-i]){    
    changeRows <- c(changeRows,paste0(allPairs[i]," > ",allPairs[j]))
  }
}
mutRows <- discardDups(changeRows)

#### run Figure_1 now.

#### run Figure 2 now