library("plyr")
library("data.table")


############## import, process VCF##############
VCFfolder <- paste0(workPath,"VCF/")
VCFfiles <- list.files(VCFfolder)
nVCFfiles <- length(VCFfiles)
allCPRC <- NULL
nRead <- 0
for(i in 1:nVCFfiles){
   print( paste0("Reading... ",VCFfiles[i]))
   Ncomments <- sum(substr(readLines(paste0(VCFfolder,VCFfiles[i])),1,1)=="#")  ## Count comment lines
   temp <- read.table(paste0(VCFfolder,VCFfiles[i]), stringsAsFactors=FALSE, fill=TRUE, header=TRUE,comment.char="", skip = Ncomments-1, sep="\t")    ## read in
   if(exists("temp")){if(length(temp)>0){nRead<-nRead+1}}
   temp$worm <- paste0("worm",i)
   temp$CPRC <- apply(temp,1,function(x) paste(x[1],x[2],x[4],x[5]))                           ## add CPRC col
   allCPRC <- c(allCPRC,levels(factor(temp$CPRC)))                                             ## add these CPRCs to growing list of all CPRCs
   name <- paste0("VCF",i)                                                      ## build variable name
   assign(name, temp)   
   rm(temp)
}
print(paste0("Read ",nRead," of ",nVCFfiles," VCF files."))

## subtract parental from VCFs
CPRC <- data.frame("CPRC" = allCPRC)                          ## move CPRC into a df, so ddply() can work on it
CPRCcount <- ddply(CPRC,.(CPRC),nrow)                         ## counts times each CPRC shows up, and puts that in a new column
colnames(CPRCcount)[2] <- "numWorms"                          ## renames the new column to something sensible

parental <- CPRCcount[which(CPRCcount$numWorm>=cutoff),]      ## "parental" variants must show up >= cutoff number of times
for(i in 1:nVCFfiles){
   name <- paste0("VCF",i)                              ## build var name to call
   temp <- get(name)                                                  ## call var
   noParents <- temp[which(is.na(match(temp$CPRC,parental$CPRC))),]   ## remove variants that match a parental
   assign(name, noParents)                                              ## reassign 
}

## concatenate
fullVCF <- NULL
for(i in 1:nVCFfiles){
   name <- paste0("VCF",i)
   temp <- get(name)
   fullVCF <- rbind(fullVCF,temp)
}

############## import, process snpEff##############
snpeffFolder <- paste0(workPath,"snpEff/")
snpeffFiles <- list.files(snpeffFolder)
nSNPfiles <- length(snpeffFiles)
allCPRC = NULL
nRead <- 0
for(i in 1:nSNPfiles){
   print( paste0("Reading... ",snpeffFiles[i]) )
   Ncomments <- sum(substr(readLines(paste0(snpeffFolder,snpeffFiles[i])),1,1)=="#")  ## Count comment lines
   temp <- read.table(paste0(snpeffFolder,snpeffFiles[i]), stringsAsFactors=FALSE, fill=TRUE, header=TRUE,comment.char="", skip = Ncomments-1, sep="\t")    ## read in
   if(exists("temp")){if(length(temp)>0){nRead<-nRead+1}}
   temp$worm <- paste0("worm",i)
   temp$CPRC <- apply(temp,1,function(x) paste(x[1],x[2],x[3],x[4]))                           ## add CPRC col
   allCPRC <- c(allCPRC,levels(factor(temp$CPRC))) 
   name <- paste0("snpEff",i)                                                      ## build variable name
   assign(name, temp)                                                                          ## assign var
   rm(temp)
}
print(paste0("Read ",nRead," of ",nSNPfiles," snpEff files."))

## subtract parental from SnpEffs
CPRC <- data.frame("CPRC" = allCPRC)                          ## move CPRC into a df, so ddply() can work on it
CPRCcount <- ddply(CPRC,.(CPRC),nrow)                         ## counts times each CPRC shows up, and puts that in a new column
colnames(CPRCcount)[2] <- "numWorms"                          ## renames the new column to something sensible

parental <- CPRCcount[which(CPRCcount$numWorm>=cutoff),]      ## "parental" variants must show up >= cutoff number of times
for(i in 1:nSNPfiles){
   name <- paste0("snpEff",i)                              ## build var name to call
   temp <- get(name)                                                  ## call var
   noParents <- temp[which(is.na(match(temp$CPRC,parental$CPRC))),]   ## remove variants that match a parental
   assign(name, noParents)                                              ## reassign 
}

## concatenate
fullsnpeff <- NULL
for(i in 1:nSNPfiles){
   name <- paste0("snpEff",i)
   temp <- get(name)
   fullsnpeff <- rbind(fullsnpeff,temp)
}
fullsnpeff$CPRC<-paste0("chr",fullsnpeff$CPRC)

####### Integrate codon and AA columns from snpeff to VCF ####
print("Merging...")
snpEffcodonHitsOnly <- fullsnpeff[which(nchar(fullsnpeff$old_AA.new_AA)!=0|nchar(fullsnpeff$Old_codon.New_codon)!=0),]
pb<-txtProgressBar(0,nrow(fullVCF),style = 3,width=100)
fullIntegratedList<-apply(fullVCF,1,function(x){
   setTxtProgressBar(pb,which(fullVCF['CPRC']==x['CPRC'])[1])
   snpeffRowIdx<-which(snpEffcodonHitsOnly['CPRC']==x['CPRC'])
   if(length(snpeffRowIdx)==0){ ## if no snpEff, just return the VCF row with empty codon and AA cols
      return(cbind(rbind(x),data.frame("Old_codon.New_codon"="","old_AA.new_AA"="",stringsAsFactors = FALSE),row.names=NULL))
   }else{
      return(cbind(rbind(x),data.frame("Old_codon.New_codon"=snpEffcodonHitsOnly$Old_codon.New_codon[snpeffRowIdx],"old_AA.new_AA"=snpEffcodonHitsOnly$old_AA.new_AA[snpeffRowIdx],stringsAsFactors = FALSE),row.names=NULL))
   }
   
})
close(pb)
fullIntegratedDF<- rbindlist(fullIntegratedList)

noIndel <- fullIntegratedDF[which((nchar(as.character(fullIntegratedDF$REF))==1) & nchar(as.character(fullIntegratedDF$ALT))==1),]  ## drop indels
dat <- data.frame("ID"=noIndel$worm ,"chr"=substring(as.character(noIndel$X.CHROM),4),"pos"=noIndel$POS,"wt_dna"=noIndel$REF,"mut_dna"=noIndel$ALT,"Old_codon.New_codon"=noIndel$Old_codon.New_codon,"old_AA.new_AA"=noIndel$old_AA.new_AA)

rm(list=setdiff(ls(), c("dat",VARS2KEEP))) 
gc()