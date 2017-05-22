library("readr")
library("seqinr")
library("plyr")


## read in data
subfolder<- paste0("F1_",getDat,"_snpeffs/deduped/") ## getDat set in organizer
datfolder <- paste0(workPath,subfolder)
datfiles <- list.files(datfolder,pattern = "[.]")
ndatfiles <- length(datfiles)
allCPRC <- NULL
nRead <- 0
for(i in 1:ndatfiles){
   print( paste0("Reading... ",datfiles[i]))
   pathToFile<-paste0(datfolder,datfiles[i])
   Ncomments <- max(which(substr(readLines(pathToFile,n=100),1,1)=="#"))
   comment <- readLines(pathToFile, n=Ncomments)
   temp <- read_delim(pathToFile, delim="\t", comment="", skip = Ncomments-1)
   if(length(which(nchar(temp$ALT)>1)) != 0){ temp <- temp[-which(nchar(temp$ALT)>1),]}
   if(length(which(nchar(temp$REF)>1)) != 0){ temp <- temp[-which(nchar(temp$REF)>1),]}
   colnames(temp)[1]<-substr(colnames(temp)[1],2,nchar(colnames(temp)[1])) ## gets rid of the # at begin of the first col name
   
   if(exists("temp")){if(length(temp)>0){nRead<-nRead+1}}
   temp$wormID <- strsplit(datfiles[i],"[.]")[[1]][1]
   temp$CPRC <- apply(temp,1,function(x) paste(x[1],x[2],x[4],x[5]))                           ## add CPRC col
   allCPRC <- c(allCPRC,levels(factor(temp$CPRC)))                                             ## add these CPRCs to growing list of all CPRCs
   name <- paste0("dat",i)                                                      ## build variable name
   assign(name, temp)   
   rm(temp)
}
print(paste0("Read ",nRead," of ",ndatfiles," snpEff files."))

## subtract parental from VCFs
CPRC <- data.frame("CPRC" = allCPRC)                          ## move CPRC into a df, so ddply() can work on it
CPRCcount <- ddply(CPRC,.(CPRC),nrow)                         ## counts times each CPRC shows up, and puts that in a new column
colnames(CPRCcount)[2] <- "numWorms"                          ## renames the new column to something sensible

parental <- CPRCcount[which(CPRCcount$numWorm>=cutoff),]      ## "parental" variants must show up >= cutoff number of times
for(i in 1:ndatfiles){
   name <- paste0("dat",i)                              ## build var name to call
   temp <- get(name)                                                  ## call var
   noParents <- temp[which(is.na(match(temp$CPRC,parental$CPRC))),]   ## remove variants that match a parental
   assign(name, noParents)                                              ## reassign 
}

## concatenate
fulldat <- NULL
for(i in 1:ndatfiles){
   name <- paste0("dat",i)
   temp <- get(name)
   if(sum(colnames(temp)=="WORMSAMPLE") == 1){
      colnames(temp)[which(colnames(temp)=="WORMSAMPLE")]<-"worm"
      }
   fulldat <- rbind(fulldat,temp)
}

## fix the old_AA.new_AA column to jive better with Fig3
print("Fixing columns...")
AAchange <-apply(fulldat,1,function(s){
   if(!grepl("/",s['EFF_Codon_Change'])| (nchar(s['EFF_Codon_Change'])!=7)){NA}  ## if EFF_Codon_Change isnt of the format '***/***', put NA in the AA change col
   else{
      codonChange <-strsplit(s['EFF_Codon_Change'],"/")[[1]]
      paste(sapply(codonChange,function(x){translate(s2c(x))}),collapse="/")     ## Otherwise, translate the codons to AA, return the */* format.
   }
})

splitRGSM<-strsplit(fulldat$worm,"[:]")
VAF<-sapply(splitRGSM,function(x){
   a<-strsplit(x[2],",")
   as.numeric(a[[1]][2])/(as.numeric(a[[1]][1])+as.numeric(a[[1]][2]))
}
)
HR<-sapply(splitRGSM,function(x){
   a<-strsplit(x[2],",")
   as.numeric(a[[1]][1])/as.numeric(a[[1]][2])
}
)
DP<-as.numeric(sapply(splitRGSM,function(x){
   x[3]
}
))
altcount<-as.numeric(sapply(splitRGSM,function(x){
   a<-strsplit(x[2],",")
   as.numeric(a[[1]][2])
}
))

fulldat<-cbind(fulldat,HR,VAF,DP,altcount,AAchange)

filtereddat <- with(fulldat, fulldat[-which(
   as.numeric(QD) < 2.5 |
   FS > 30 |
   SOR > 3 |
   MQ < 55 |
   as.numeric(MQRankSum) < -2.5 |
   as.numeric(MQRankSum) > 2.5 |
   as.numeric(ReadPosRankSum) < -2 |
   as.numeric(ReadPosRankSum) > 2 |
   altcount < 2.5
),])


## get things into uniform structure, with only the necessary columns
dat<-with(filtereddat,data.frame("altCount"=altcount,"DP"=DP,"HR"=HR,"VAF"=VAF, "ID"=wormID,"CPRC"=CPRC,"chr"=CHROM,"pos"=POS,"wt_dna"=REF,"mut_dna"=ALT, "Old_codon.New_codon"=EFF_Codon_Change, "old_AA.new_AA"=AAchange))

altCountCutoff <- 8 ## higher = more conservative
HRcutoff<- 3.5 # lower = more conservative


dat<-dat[which( dat$altCount>altCountCutoff & dat$HR<HRcutoff ),]
rm(list=setdiff(ls(), c("dat",VARS2KEEP))) 
gc()








